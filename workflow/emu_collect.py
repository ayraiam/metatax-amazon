#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, json, sys
from pathlib import Path
from typing import Dict, List, Any, Tuple, Optional
import pandas as pd
import numpy as np

def find_first(path: Path, patterns: List[str]) -> Optional[Path]:
    for pat in patterns:
        hits = sorted(path.glob(pat))
        if hits: return hits[0]
        hits = sorted((p for p in path.glob(f"*/{pat}")))
        if hits: return hits[0]
    return None

def read_table_guess_sep(p: Path, compression="infer") -> pd.DataFrame:
    try:
        df = pd.read_csv(p, sep="\t", compression=compression, engine="python")
        if df.shape[1] == 1:
            df = pd.read_csv(p, sep=",", compression=compression, engine="python")
        return df
    except Exception:
        return pd.read_csv(p, delim_whitespace=True, compression=compression, engine="python")

def load_abundance_raw(sample_dir: Path) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    p = sample_dir / "abundance.tsv"
    if not p.exists():
        p = find_first(sample_dir, ["*abundance*.tsv","*abundance*.csv","*abundance*.tsv.gz","*abundance*.csv.gz"])
    if p is None:
        return pd.DataFrame(), {}
    df = read_table_guess_sep(p, compression="infer")
    df.insert(0, "file", sample_dir.name)
    abund_dict = {col: df[col].tolist() for col in df.columns}
    return df, abund_dict

def load_read_assignment_distributions(sample_dir: Path) -> Tuple[pd.DataFrame, Dict[str, Dict[str, Dict[str, float]]]]:
    p = find_first(sample_dir, ["*read-assignment-distributions.tsv","*read-assignment-distributions.tsv.gz"])
    if p is None:
        return pd.DataFrame(), {}
    df = read_table_guess_sep(p, compression="infer")
    if df.shape[1] == 0:
        return pd.DataFrame(), {}
    first = df.columns[0]
    df.rename(columns={first: "read_id"}, inplace=True)
    taxon_cols = [c for c in df.columns if c != "read_id"]
    for c in taxon_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    sample = sample_dir.name
    nested: Dict[str, Dict[str, Dict[str, float]]] = {sample: {}}
    for _, row in df.iterrows():
        rid = str(row["read_id"]).strip()
        if not rid or rid.lower() == "nan": continue
        nested[sample][rid] = {str(tid): float(row[tid]) for tid in taxon_cols}
    return df, nested

def mapping_stats_fractional(assign_df: pd.DataFrame, total_reads_sidecar: int) -> Dict[str, Any]:
    if assign_df.empty:
        total = int(total_reads_sidecar) if total_reads_sidecar > 0 else 0
        return dict(total_reads=total, assigned_reads=0.0, assigned_frac=0.0,
                    unassigned_reads=float(total), unassigned_frac=(0.0 if total == 0 else 1.0))
    taxon_cols = [c for c in assign_df.columns if c != "read_id"]
    if not taxon_cols:
        total = int(total_reads_sidecar) if total_reads_sidecar > 0 else int(assign_df.shape[0])
        return dict(total_reads=total, assigned_reads=0.0, assigned_frac=0.0,
                    unassigned_reads=float(total), unassigned_frac=(0.0 if total == 0 else 1.0))
    mat = assign_df[taxon_cols].to_numpy(dtype=float, copy=False)
    row_sums = np.nansum(mat, axis=1); row_sums = np.clip(row_sums, 0.0, 1.0)
    fractional_assigned = float(row_sums.sum())
    total_rows = assign_df.shape[0]
    total = int(total_reads_sidecar) if total_reads_sidecar > 0 else int(total_rows)
    assigned_frac = (fractional_assigned / total) if total > 0 else 0.0
    unassigned_reads = max(0.0, float(total) - fractional_assigned)
    unassigned_frac = max(0.0, 1.0 - assigned_frac)
    return dict(
        total_reads=int(total),
        assigned_reads=round(fractional_assigned, 6),
        assigned_frac=round(assigned_frac, 6),
        unassigned_reads=round(unassigned_reads, 6),
        unassigned_frac=round(unassigned_frac, 6),
    )

def main():
    ap = argparse.ArgumentParser(description="Collect Emu abundance + (optional) read-assignment distributions.")
    ap.add_argument("--runs-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--dictdir", default="metadata")
    ap.add_argument("--min-prob", type=float, default=0.5, help="(Ignored; fractional mapping stats.)")
    ap.add_argument("--no-json", action="store_true")
    ap.add_argument("--skip-assign", action="store_true", help="Skip reading read-assignment distributions (saves RAM/IO).")
    args = ap.parse_args()

    if "--min-prob" in sys.argv:
        sys.stderr.write("Note: --min-prob is ignored; using fractional mapping stats.\n")

    runs_dir = Path(args.runs_dir)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    dictdir = Path(args.dictdir); dictdir.mkdir(parents=True, exist_ok=True)

    abundance_dict: Dict[str, Dict[str, List[Any]]] = {}
    read_assign_dict: Dict[str, Dict[str, Dict[str, float]]] = {}

    all_abund = []
    map_rows = []

    for sdir in sorted([p for p in runs_dir.glob("*") if p.is_dir()]):
        sname = sdir.name

        abund_df, abund_d = load_abundance_raw(sdir)
        if not abund_df.empty: all_abund.append(abund_df)
        if abund_d: abundance_dict[sname] = abund_d

        # total reads from sidecar (if present)
        total_reads_sidecar = 0
        sidecar = sdir / "input_reads.tsv"
        if sidecar.exists():
            try:
                tmp = pd.read_csv(sidecar, sep="\t")
                if tmp.shape[0] >= 1 and tmp.shape[1] >= 2:
                    total_reads_sidecar = int(tmp.iloc[0, 1])
            except Exception:
                pass

        if args.skip-assign:
            # If skipping, we can't compute assigned/unassigned; keep totals only
            stats = dict(
                total_reads=int(total_reads_sidecar),
                assigned_reads=np.nan, assigned_frac=np.nan,
                unassigned_reads=np.nan, unassigned_frac=np.nan
            )
        else:
            assign_df, assign_d = load_read_assignment_distributions(sdir)
            if assign_d:
                read_assign_dict[sname] = assign_d[sname]
            stats = mapping_stats_fractional(assign_df, total_reads_sidecar)

        stats_row = {"file": sname, **stats}
        map_rows.append(stats_row)

    # Write batch-scoped outputs
    abund_out = outdir / "abundance_combined.tsv"
    if all_abund:
        pd.concat(all_abund, ignore_index=True).to_csv(abund_out, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["file","tax_id","abundance","species"]).to_csv(abund_out, sep="\t", index=False)

    map_out = outdir / "mapping_stats.tsv"
    pd.DataFrame(map_rows, columns=["file","total_reads","assigned_reads","assigned_frac","unassigned_reads","unassigned_frac"]).to_csv(map_out, sep="\t", index=False)

    if not args.no_json:
        (dictdir / "abundance_dict.json").write_text(json.dumps(abundance_dict))
        if not args.skip_assign:
            (dictdir / "read_assignment_dict.json").write_text(json.dumps(read_assign_dict))
        sys.stderr.write(">>> Wrote JSON dicts.\n")

    sys.stderr.write(f">>> Wrote abundance table -> {abund_out}\n")
    sys.stderr.write(f">>> Wrote mapping stats   -> {map_out}\n")

if __name__ == "__main__":
    main()
