from svmu2.orchestration.parse import parse
from svmu2.orchestration.synteny import resolve_synteny
import argparse
from pathlib import Path
import pandas as pd
import gc

def parse_args():
    parser = argparse.ArgumentParser(description="Block-first PAM orthology engine.")
    parser.add_argument("--ref-cfd", required=True, type=Path)
    parser.add_argument("--query-cfd", required=True, type=Path)
    parser.add_argument("--delta", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path, help="Directory to store all CSV outputs")
    parser.add_argument("--tol", type=int, default=1000)
    
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    return args

REQUIRED_CFD_COLUMNS = {"protospacerID", "ref_chr", "ref_start", "ref_end"}

def load_cfd_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"protospacerID": str, "ref_chr": str, "ref_start": int, "ref_end": int})
    if REQUIRED_CFD_COLUMNS - set(df.columns):
        raise ValueError(f"Missing columns in {path}")
    df["midpoint"] = (df["ref_start"] + df["ref_end"]) / 2.0
    return df

def assign_synteny(df, trees, mode="reference"):
    df_out = df.copy() 
    block_ids, yhats = [], []

    for row in df.itertuples(index=False):
        tree = trees.get(row.ref_chr)
        x = row.midpoint 

        if tree is None:
            block_ids.append(None)
            if mode == "reference": yhats.append(None)
            continue

        hits = tree.at(x)
        if not hits:
            block_ids.append(None)
            if mode == "reference": yhats.append(None)
        elif len(hits) == 1:
            hit_data = next(iter(hits)).data
            block_ids.append(hit_data.index) 
            if mode == "reference": yhats.append((hit_data.slope * x) + hit_data.y_intercept)
        else:
            # Drop duplications/paralogs to enforce 1:1 strict orthology
            block_ids.append(None)
            if mode == "reference": yhats.append(None)

    df_out["syntenic_block"] = block_ids
    if mode == "reference": df_out["yhat"] = yhats
    return df_out

def merge_and_filter(ref_df, qry_df, out_dir, tol):
    print("[INFO] Executing block-wise merge and tolerance filtering...")
    
    ref_mapped = ref_df.dropna(subset=['syntenic_block'])
    qry_mapped = qry_df.dropna(subset=['syntenic_block'])

    # 1. Perform the merge (this temporarily creates the combinatorial explosion)
    merged = ref_mapped.merge(qry_mapped, on="protospacerID", suffixes=('_ref', '_qry'))
    
    # 2. Calculate projection distances
    merged['proj_dist'] = (merged['midpoint_qry'] - merged['yhat']).abs()
    
    # 3. Collapse down to the single best hit per PAM
    # Sort so the smallest distance is first, then keep only that first row for each ID
    merged_best = merged.sort_values('proj_dist').drop_duplicates(subset='protospacerID', keep='first')
    
    # 4. Split based on tolerance using the deduplicated dataframe
    conserved = merged_best[merged_best['proj_dist'] <= tol]
    failed_tol = merged_best[merged_best['proj_dist'] > tol]
    
    # Optional: You might want to save the "rescued" multi-mappers vs strict 1:1s 
    # but for your current SVMU benchmarking, this keeps it clean.
    conserved.to_csv(out_dir / "conserved_pams.csv", index=False)
    failed_tol.to_csv(out_dir / "failed_tolerance_pams.csv", index=False)
    
    return conserved, failed_tol
    
def main():
    args = parse_args()
    ref_cfd = load_cfd_table(args.ref_cfd)
    qry_cfd = load_cfd_table(args.query_cfd)
    
    _, primary_alns = parse(args.delta)
    resolve_synteny(primary_alignments=primary_alns, breakpoint_map=None)

    ref_trees = {aln.reference: aln.reference_synteny_tree for aln in primary_alns.values()}
    qry_trees = {aln.query: aln.query_synteny_tree for aln in primary_alns.values()}

    ref_df = assign_synteny(ref_cfd, ref_trees, mode="reference")
    qry_df = assign_synteny(qry_cfd, qry_trees, mode="query")

    # Isolate unmapped PAMs for diagnostics
    unmapped_ref = ref_df[ref_df["syntenic_block"].isna()]
    unmapped_qry = qry_df[qry_df["syntenic_block"].isna()]
    
    unmapped_ref.to_csv(args.out_dir / "unmapped_ref.csv", index=False)
    unmapped_qry.to_csv(args.out_dir / "unmapped_qry.csv", index=False)

    print(f"Diagnostics: {len(unmapped_ref)} Ref PAMs and {len(unmapped_qry)} Query PAMs lacked synteny.")

    merge_and_filter(ref_df, qry_df, args.out_dir, args.tol)

if __name__ == "__main__":
    main()