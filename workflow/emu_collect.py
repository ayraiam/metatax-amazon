#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Tuple, Optional
import pandas as pd
import numpy as np

def find_first(path: Path, patterns: List[str]) -> Optional[Path]:
    """Return first matching path under 'path' (maxdepth=2), else None."""
    for pat in patterns:
        hits = sorted(path.glob(pat))
        if hits:
            return hits[0]
        # also search one level deeper
        hits = sorted((p for p in path.glob(f"*/{pat}")))
        if hits:
            return hits[0]
    return None

def read_table_guess_sep(p: Path, compression="infer") -> pd.DataFrame:
    """Robust reader: try tab; if single column, try comma; fallback to whitespace."""
    try:
        df = pd.read_csv(p, sep="\t", compression=compression, engine="python")
        if df.shape[1] == 1:
            df = pd.read_csv(p, sep=",", compression=compression, engine="python")
        return df
    except Exception:
        df = pd.read_csv(p, delim_whitespace=True, compression=compression, engine="python")
        return df

def load_abundance_raw(sample_dir: Path) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Load abundance table WITHOUT standardizing columns.
    Keep every original column (e.g., tax_id, abundance, species, genus, ...).
    Add a 'file' column (sample dir name).
    Return (DataFrame, dict-of-lists with original column names).
    """
    p = None
    if (sample_dir / "abundance.tsv").exists():
        p = sample_dir / "abundance.tsv"
    else:
        p = find_first(sample_dir, [
            "*abundance*.tsv", "*abundance*.csv",
            "*abundance*.tsv.gz", "*abundance*.csv.gz"
        ])
    if p is None:
        return pd.DataFrame(), {}

    df = read_table_guess_sep(p, compression="infer")

    # Add sample name as 'file'
    df.insert(0, "file", sample_dir.name)

    # Dict: inner keys are EXACTLY the column names in this abundance file
    abund_dict = {col: df[col].tolist() for col in df.columns}

    return df, abund_dict

def load_read_assignment_distributions(sample_dir: Path) -> Tuple[pd.DataFrame, Dict[str, Dict[str, Dict[str, float]]]]:
    """
    Load the read-assignment distributions table:
      - First column is the read ID (header may be blank/Unnamed)
      - Remaining columns are taxon IDs
      - Cell values are probabilities
    Return (DataFrame, nested dict { sample: { read_id: { taxon_id: prob }}}).
    """
    p = find_first(sample_dir, [
        "*read-assignment-distributions.tsv",
        "*read-assignment-distributions.tsv.gz"
    ])
    if p is None:
        return pd.DataFrame(), {}

    df = read_table_guess_sep(p, compression="infer")
    if df.shape[1] == 0:
        return pd.DataFrame(), {}

    # Force first column to be called 'read_id'
    first = df.columns[0]
    df.rename(columns={first: "read_id"}, inplace=True)

    # Prob columns = all others (taxon IDs already in headers)
    taxon_cols = [c for c in df.columns if c != "read_id"]

    # Convert to numeric probabilities (NaN → 0)
    for c in taxon_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    # Nested dict
    sample = sample_dir.name
    nested: Dict[str, Dict[str, Dict[str, float]]] = {sample: {}}
    for _, row in df.iterrows():
        rid = str(row["read_id"]).strip()
        if not rid or rid.lower() == "nan":
            continue
        nested[sample][rid] = {str(tid): float(row[tid]) for tid in taxon_cols}

    return df, nested

# ------------------- FRACTIONAL mapping stats -------------------
def mapping_stats_fractional(assign_df: pd.DataFrame, total_reads_sidecar: int) -> Dict[str, Any]:
    """
    Compute mapping stats *fractionally* from the read-assignment posterior matrix:
      - Sum probabilities across taxa per read (row-sum), clip to [0,1]
      - Sum across reads -> fractional assigned reads
      - Fractions relative to total_reads (prefer sidecar; else #rows)
    """
    if assign_df.empty:
        total = int(total_reads_sidecar) if total_reads_sidecar > 0 else 0
        return dict(
            total_reads=total,
            assigned_reads=0.0,
            assigned_frac=0.0,
            unassigned_reads=float(total),
            unassigned_frac=0.0 if total == 0 else 1.0
        )

    taxon_cols = [c for c in assign_df.columns if c != "read_id"]
    if not taxon_cols:
        total = int(total_reads_sidecar) if total_reads_sidecar > 0 else int(assign_df.shape[0])
        return dict(
            total_reads=total,
            assigned_reads=0.0,
            assigned_frac=0.0,
            unassigned_reads=float(total),
            unassigned_frac=0.0 if total == 0 else 1.0
        )

    # Row sums of probabilities; defensive clipping & NaN handling
    mat = assign_df[taxon_cols].to_numpy(dtype=float, copy=False)
    row_sums = np.nansum(mat, axis=1)
    row_sums = np.clip(row_sums, 0.0, 1.0)

    fractional_assigned = float(row_sums.sum())

    total_rows = assign_df.shape[0]
    total = int(total_reads_sidecar) if total_reads_sidecar > 0 else int(total_rows)

    assigned_frac = (fractional_assigned / total) if total > 0 else 0.0
    unassigned_reads = max(0.0, float(total) - fractional_assigned)
    unassigned_frac = max(0.0, 1.0 - assigned_frac)

    # Round display values (preserve transparency)
    return dict(
        total_reads=int(total),
        assigned_reads=round(fractional_assigned, 6),
        assigned_frac=round(assigned_frac, 6),
        unassigned_reads=round(unassigned_reads, 6),
        unassigned_frac=round(unassigned_frac, 6),
    )

def main():
    ap = argparse.ArgumentParser(description="Collect Emu abundance + read-assignment distributions.")
    ap.add_argument("--runs-dir", required=True, help="Path to results/emu_runs")
    ap.add_argument("--outdir", required=True, help="Where to write tables (e.g., results/tables)")
    ap.add_argument("--dictdir", default="metadata", help="Where to dump JSON dicts")
    ap.add_argument("--min-prob", type=float, default=0.5, help="(Ignored) Kept for backward-compatibility.")
    ap.add_argument("--no-json", action="store_true", help="Do not write JSON dict files.")
    ap.add_argument("--skip-assign", action="store_true", help="Skip loading read-assignment matrices.")
    args = ap.parse_args()

    if "--min-prob" in sys.argv:
        sys.stderr.write("Note: --min-prob is ignored; using fractional mapping stats.\n")

    runs_dir = Path(args.runs_dir)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    dictdir = Path(args.dictdir); dictdir.mkdir(parents=True, exist_ok=True)

    abundance_dict: Dict[str, Dict[str, List[Any]]] = {}
    read_assign_dict: Dict[str, Dict[str, Dict[str, float]]] = {}

    all_abund, map_rows = [], []

    for sdir in sorted([p for p in runs_dir.glob("*") if p.is_dir()]):
        sname = sdir.name

        # 1) abundance (raw)
        abund_df, abund_d = load_abundance_raw(sdir)
        if not abund_df.empty:
            all_abund.append(abund_df)
        if abund_d:
            abundance_dict[sname] = abund_d

        # 2) (optional) read-assignment distributions
        if args.skip_assign:
            assign_df, assign_d = pd.DataFrame(), {}
        else:
            assign_df, assign_d = load_read_assignment_distributions(sdir)

        if assign_d:
            read_assign_dict[sname] = assign_d[sname]

        # 3) total reads from sidecar (if present)
        total_reads_sidecar = 0
        sidecar = sdir / "input_reads.tsv"
        if sidecar.exists():
            try:
                tmp = pd.read_csv(sidecar, sep="\t")
                if tmp.shape[0] >= 1 and tmp.shape[1] >= 2:
                    total_reads_sidecar = int(tmp.iloc[0, 1])
            except Exception:
                pass

        # 4) fractional mapping stats (if skip-assign: becomes 0 assigned)
        stats = mapping_stats_fractional(assign_df, total_reads_sidecar)
        stats_row = {"file": sname, **stats}
        map_rows.append(stats_row)

    # Write abundance_combined.tsv (unchanged)
    abund_out = outdir / "abundance_combined.tsv"
    if all_abund:
        comb = pd.concat(all_abund, ignore_index=True)
        comb.to_csv(abund_out, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["file", "tax_id", "abundance", "species"]).to_csv(abund_out, sep="\t", index=False)

    # Write mapping_stats.tsv (unchanged)
    map_out = outdir / "mapping_stats.tsv"
    pd.DataFrame(map_rows, columns=[
        "file", "total_reads", "assigned_reads", "assigned_frac", "unassigned_reads", "unassigned_frac"
    ]).to_csv(map_out, sep="\t", index=False)

    # Dump dicts (optional)
    if not args.no_json:
        (dictdir / "abundance_dict.json").write_text(json.dumps(abundance_dict))
        (dictdir / "read_assignment_dict.json").write_text(json.dumps(read_assign_dict))
        sys.stderr.write(f">>> Wrote dicts           -> {dictdir/'abundance_dict.json'}, {dictdir/'read_assignment_dict.json'}\n")
    else:
        sys.stderr.write(">>> Skipped writing dict JSONs (--no-json)\n")

    sys.stderr.write(f">>> Wrote abundance table -> {abund_out}\n")
    sys.stderr.write(f">>> Wrote mapping stats   -> {map_out}\n")

if __name__ == "__main__":
    main()
