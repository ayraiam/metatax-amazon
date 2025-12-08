#!/usr/bin/env python3
import argparse, sys
from pathlib import Path
import gzip

def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def parse_args():
    ap = argparse.ArgumentParser(
        description="Build Emu LSU DB (Asco+Basidio) from SILVA LSURef NR99 taxonomy FASTA."
    )
    ap.add_argument("--fasta", required=True,
                    help="SILVA LSURef taxonomy FASTA (e.g. SILVA_138.2_LSURef_NR99_tax_silva.fasta[.gz])")
    ap.add_argument("--out-prefix", required=True,
                    help="Output prefix, e.g. refdb/silva_lsu_138/lsu_silva_asco_basi")
    return ap.parse_args()

def parse_silva_header(hdr: str):
    """
    Example SILVA header:
    >AC152122.8318.11813 Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea;Fungi;Dikarya;Ascomycota;Saccharomycotina;Saccharomycetes;Saccharomycetales;Debaryomycetaceae;Scheffersomyces;Scheffersomyces stipitis
    """
    hdr = hdr.lstrip(">").strip()
    if not hdr:
        raise ValueError("Empty SILVA header")
    parts = hdr.split(maxsplit=1)
    seq_id = parts[0]
    tax_str = parts[1] if len(parts) > 1 else ""

    if not tax_str:
        return seq_id, {}, tax_str

    tokens = [t.strip() for t in tax_str.split(";") if t.strip()]

    # ---------------------------------------------------------
    # find 'Fungi' in the chain and use offsets from there
    # ---------------------------------------------------------
    ranks = {}
    try:
        k_idx = tokens.index("Fungi")
    except ValueError:
        # not fungal at all
        return seq_id, ranks, tax_str

    fungal = tokens[k_idx:]                    # ['Fungi','Dikarya','Ascomycota',...]
    if not fungal:
        return seq_id, ranks, tax_str

    # kingdom = Fungi
    ranks["kingdom"] = fungal[0]

    # phylum, class, order, family, genus, species with safe indexing
    if len(fungal) > 2:
        ranks["phylum"] = fungal[2]            # e.g. Ascomycota
    if len(fungal) > 4:
        ranks["class"]  = fungal[4]            # e.g. Saccharomycetes
    if len(fungal) > 5:
        ranks["order"]  = fungal[5]            # e.g. Saccharomycetales
    if len(fungal) > 6:
        ranks["family"] = fungal[6]            # e.g. Debaryomycetaceae
    if len(fungal) > 7:
        ranks["genus"]  = fungal[7]            # e.g. Scheffersomyces
    if len(fungal) > 8:
        ranks["species"]= fungal[8]            # e.g. Scheffersomyces stipitis

    return seq_id, ranks, tax_str

def main():
    args = parse_args()
    fasta_path = Path(args.fasta)
    out_prefix = Path(args.out_prefix)

    if not fasta_path.exists():
        sys.stderr.write(f"FASTA not found: {fasta_path}\n")
        sys.exit(1)

    kept = []
    # ---------------------------------------------------------
    # NEW: deduplicated taxonomy using dict
    # ---------------------------------------------------------
    taxon2id = {}
    tax_rows = []
    seq2tax_rows = []
    tax_id_counter = 1

    with open_maybe_gz(fasta_path) as fh:
        cur_h = None
        cur_seq = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_h is not None:
                    # process previous
                    seq_id, ranks, _ = parse_silva_header(cur_h)
                    kingdom = ranks.get("kingdom", "")
                    phylum  = ranks.get("phylum", "")

                    # keep only fungal Ascomycota / Basidiomycota
                    if kingdom == "Fungi" and phylum in ("Ascomycota", "Basidiomycota"):
                        seq = "".join(cur_seq)
                        kept.append((seq_id, seq))

                        species = ranks.get("species", "")
                        genus   = ranks.get("genus", "")
                        family  = ranks.get("family", "")
                        order   = ranks.get("order", "")
                        clazz   = ranks.get("class", "")
                        # phylum, kingdom already defined

                        tax_key = (
                            species,
                            genus,
                            family,
                            order,
                            clazz,
                            phylum,
                            kingdom,
                        )

                        tid = taxon2id.get(tax_key)
                        if tid is None:
                            tid = tax_id_counter
                            tax_id_counter += 1
                            taxon2id[tax_key] = tid
                            tax_rows.append([
                                str(tid),
                                species,
                                genus,
                                family,
                                order,
                                clazz,
                                phylum,
                                kingdom
                            ])

                        seq2tax_rows.append([seq_id, str(tid)])

                cur_h = line
                cur_seq = []
            else:
                cur_seq.append(line)
        # final record
        if cur_h is not None:
            seq_id, ranks, _ = parse_silva_header(cur_h)
            kingdom = ranks.get("kingdom", "")
            phylum  = ranks.get("phylum", "")
            if kingdom == "Fungi" and phylum in ("Ascomycota", "Basidiomycota"):
                seq = "".join(cur_seq)
                kept.append((seq_id, seq))

                species = ranks.get("species", "")
                genus   = ranks.get("genus", "")
                family  = ranks.get("family", "")
                order   = ranks.get("order", "")
                clazz   = ranks.get("class", "")

                tax_key = (
                    species,
                    genus,
                    family,
                    order,
                    clazz,
                    phylum,
                    kingdom,
                )

                tid = taxon2id.get(tax_key)
                if tid is None:
                    tid = tax_id_counter
                    tax_id_counter += 1
                    taxon2id[tax_key] = tid
                    tax_rows.append([
                        str(tid),
                        species,
                        genus,
                        family,
                        order,
                        clazz,
                        phylum,
                        kingdom
                    ])

                seq2tax_rows.append([seq_id, str(tid)])

    if not kept:
        sys.stderr.write("No fungal Ascomycota/Basidiomycota sequences found in SILVA LSU FASTA.\n")
        sys.exit(1)

    out_fasta   = out_prefix.with_suffix(".fasta")
    out_seq2tax = out_prefix.with_suffix(".seq2tax.map.tsv")
    out_tax     = out_prefix.with_suffix(".taxonomy.tsv")

    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    # write FASTA
    with out_fasta.open("w") as fh:
        for sid, seq in kept:
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")

    # write seq2tax
    with out_seq2tax.open("w") as fh:
        for sid, tid in seq2tax_rows:
            fh.write(f"{sid}\t{tid}\n")

    # write taxonomy
    with out_tax.open("w") as fh:
        fh.write("tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\n")
        for row in tax_rows:
            fh.write("\t".join(row) + "\n")

    sys.stderr.write("Wrote:\n")
    sys.stderr.write(f"  {out_fasta}\n")
    sys.stderr.write(f"  {out_seq2tax}\n")
    sys.stderr.write(f"  {out_tax}\n")

if __name__ == "__main__":
    main()
