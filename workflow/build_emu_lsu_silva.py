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
    Example SILVA header (one possible style):
    >AF177760.1.1752 Eukaryota;Fungi;Ascomycota;Saccharomycetes;...

    We do NOT assume fixed positions; we return:
      - seq_id
      - ranks dict (best-effort mapping to kingdom/phylum/...)
      - tokens: the full list of taxonomy tokens
    """
    hdr = hdr.lstrip(">").strip()
    if not hdr:
        raise ValueError("Empty SILVA header")
    parts = hdr.split(maxsplit=1)
    seq_id = parts[0]
    tax_str = parts[1] if len(parts) > 1 else ""

    if not tax_str:
        return seq_id, {}, []

    tokens = [t.strip() for t in tax_str.split(";") if t.strip()]

    # Best-effort mapping: we keep the first 7 tokens as
    # kingdom, phylum, class, order, family, genus, species
    ranks = {}
    if len(tokens) >= 1:
        ranks["kingdom"] = tokens[0]
    if len(tokens) >= 2:
        ranks["phylum"]  = tokens[1]
    if len(tokens) >= 3:
        ranks["class"]   = tokens[2]
    if len(tokens) >= 4:
        ranks["order"]   = tokens[3]
    if len(tokens) >= 5:
        ranks["family"]  = tokens[4]
    if len(tokens) >= 6:
        ranks["genus"]   = tokens[5]
    if len(tokens) >= 7:
        ranks["species"] = tokens[6]

    return seq_id, ranks, tokens

def is_fungal_asco_basidio(tokens):
    """
    more flexible detection of fungal Asco/Basidio

    We keep the record if:
      - One of the tokens is 'Fungi' (or ends with 'Fungi')
      - AND one of the tokens is 'Ascomycota' or 'Basidiomycota'
    """
    if not tokens:
        return False, "", ""

    # detect 'Fungi' anywhere
    fungi_like = [t for t in tokens if t == "Fungi" or t.endswith("Fungi")]
    if not fungi_like:
        return False, "", ""

    # detect phylum token
    phylum = None
    for t in tokens:
        if t in ("Ascomycota", "Basidiomycota"):
            phylum = t
            break

    if phylum is None:
        return False, "", ""

    kingdom = fungi_like[0]
    return True, kingdom, phylum

def main():
    args = parse_args()
    fasta_path = Path(args.fasta)
    out_prefix = Path(args.out_prefix)

    if not fasta_path.exists():
        sys.stderr.write(f"FASTA not found: {fasta_path}\n")
        sys.exit(1)

    kept = []
    tax_rows = []
    seq2tax_rows = []

    tax_id = 1

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
                    seq_id, ranks, tokens = parse_silva_header(cur_h)

                    keep, kingdom, phylum = is_fungal_asco_basidio(tokens)
                    if keep:
                        seq = "".join(cur_seq)
                        kept.append((seq_id, seq))

                        # override or fill ranks using detected kingdom/phylum
                        ranks["kingdom"] = kingdom
                        ranks["phylum"]  = phylum

                        tax_rows.append([
                            str(tax_id),
                            ranks.get("species", ""),
                            ranks.get("genus", ""),
                            ranks.get("family", ""),
                            ranks.get("order", ""),
                            ranks.get("class", ""),
                            ranks.get("phylum", ""),
                            ranks.get("kingdom", "")
                        ])
                        seq2tax_rows.append([seq_id, str(tax_id)])
                        tax_id += 1

                cur_h = line
                cur_seq = []
            else:
                cur_seq.append(line)

        # final record
        if cur_h is not None:
            seq_id, ranks, tokens = parse_silva_header(cur_h)
            keep, kingdom, phylum = is_fungal_asco_basidio(tokens)
            if keep:
                seq = "".join(cur_seq)
                kept.append((seq_id, seq))

                ranks["kingdom"] = kingdom
                ranks["phylum"]  = phylum

                tax_rows.append([
                    str(tax_id),
                    ranks.get("species", ""),
                    ranks.get("genus", ""),
                    ranks.get("family", ""),
                    ranks.get("order", ""),
                    ranks.get("class", ""),
                    ranks.get("phylum", ""),
                    ranks.get("kingdom", "")
                ])
                seq2tax_rows.append([seq_id, str(tax_id)])
                tax_id += 1

    if not kept:
        sys.stderr.write(
            "No fungal Ascomycota/Basidiomycota sequences found in SILVA LSU FASTA "
            "(after relaxed parsing). Check the FASTA headers / version.\n"
        )
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
