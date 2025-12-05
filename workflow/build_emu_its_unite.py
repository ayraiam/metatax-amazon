#!/usr/bin/env python3
import sys, argparse
from pathlib import Path

def parse_args():
    ap = argparse.ArgumentParser(
        description="Build Emu ITS DB (Asco+Basidio) from UNITE sh_general_release_dynamic FASTA."
    )
    ap.add_argument("--fasta", required=True,
                    help="UNITE sh_general_release_dynamic_*.fasta")
    ap.add_argument("--out-prefix", required=True,
                    help="Output prefix, e.g. refdb/unite_its_2025/its_unite_asco_basi")
    return ap.parse_args()

def parse_unite_header(header: str):
    """
    Example:
    >Abrothallus_subhalei|MT153946|SH1227328.10FU|refs|k__Fungi;p__Ascomycota;...
    """
    header = header.lstrip(">")
    parts = header.split("|")

    if len(parts) < 5:
        raise ValueError(f"Unexpected UNITE header format: {header}")

    name_part = parts[0]
    acc       = parts[1]        # using this as seq_id
    tax_str   = parts[4]        # k__Fungi;p__Ascomycota;...

    tokens = [t for t in tax_str.split(";") if t]
    ranks = {}
    for t in tokens:
        if "__" in t:
            r, val = t.split("__", 1)
            ranks[r] = val
    return {
        "raw_header": header,
        "seq_id": acc,
        "name": name_part,
        "ranks": ranks
    }

# ==========================================================
# More flexible detection of Ascomycota/Basidiomycota
# ==========================================================
def is_asco_or_basidio(phylum: str) -> bool:
    """Return True if phylum belongs to Ascomycota or Basidiomycota,
    allowing variants like Ascomycota_incertae_sedis or BasidiomycotaXYZ."""
    if not phylum:
        return False
    return (
        phylum.startswith("Ascomycota") or
        phylum.startswith("Basidiomycota")
    )

def main():
    args = parse_args()
    fasta_path = Path(args.fasta)
    out_prefix = Path(args.out_prefix)

    if not fasta_path.exists():
        sys.stderr.write(f"FASTA not found: {fasta_path}\n")
        sys.exit(1)

    # read FASTA
    seqs = []
    cur_h = None
    cur_seq = []
    with fasta_path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_h is not None:
                    seqs.append((cur_h, "".join(cur_seq)))
                cur_h = line
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_h is not None:
            seqs.append((cur_h, "".join(cur_seq)))

    kept_seqs = []
    tax_rows = []
    seq2tax_rows = []

    tax_id = 1
    for hdr, seq in seqs:
        info = parse_unite_header(hdr)
        ranks = info["ranks"]

        phylum = ranks.get("p", "")

        if not is_asco_or_basidio(phylum):
            continue

        kingdom = ranks.get("k", "Fungi")
        clazz   = ranks.get("c", "")
        order   = ranks.get("o", "")
        family  = ranks.get("f", "")
        genus   = ranks.get("g", "")
        species = ranks.get("s", "")

        seq_id = info["seq_id"]

        kept_seqs.append((seq_id, seq))

        tax_rows.append([
            str(tax_id),
            species,
            genus,
            family,
            order,
            clazz,
            phylum,
            kingdom
        ])
        seq2tax_rows.append([seq_id, str(tax_id)])

        tax_id += 1

    if not kept_seqs:
        sys.stderr.write("No Ascomycota/Basidiomycota sequences found; check input.\n")
        sys.exit(1)

    out_fasta   = out_prefix.with_suffix(".fasta")
    out_seq2tax = out_prefix.with_suffix(".seq2tax.map.tsv")
    out_tax     = out_prefix.with_suffix(".taxonomy.tsv")

    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    # write FASTA
    with out_fasta.open("w") as fh:
        for sid, seq in kept_seqs:
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
