from svmu2.orchestration.parse import parse
from svmu2.orchestration.synteny import resolve_synteny
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import argparse
from pathlib import Path
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="Minimal multi-tolerance Venn plotter.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--ref-cfd", required=True, type=Path, help="Reference/ISO1 CFD CSV")
    parser.add_argument("--query-cfd", required=True, type=Path, help="Query/BL54591 CFD CSV")
    parser.add_argument("--delta", required=True, type=Path, help="nucmer delta file")
    parser.add_argument("--figures", required=True, type=Path, help="Output dir to store the figures in")

    args = parser.parse_args()
    args.figures.mkdir(parents=True, exist_ok=True)
    return args

REQUIRED_CFD_COLUMNS = {"protospacerID", "ref_chr", "ref_start", "ref_end"}

def load_cfd_table(path: Path) -> pd.DataFrame:
    dtypes = {"protospacerID": str, "ref_chr": str, "ref_start": int, "ref_end": int}
    df = pd.read_csv(path, dtype=dtypes)

    if missing := REQUIRED_CFD_COLUMNS - set(df.columns):
        raise ValueError(f"{path} CFD table is missing columns: {sorted(missing)}")
    
    df["midpoint"] = (df["ref_start"] + df["ref_end"]) / 2.0
    return df

def assign_synteny(df, trees, mode="reference"):
    df_out = df.copy() 
    block_ids = []
    
    if mode == "reference":
        yhats = []

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
            block_ids.append([hit.data.index for hit in hits])
            if mode == "reference": yhats.append([(hit.data.slope * x) + hit.data.y_intercept for hit in hits])

    df_out["syntenic_block"] = block_ids
    if mode == "reference":
        df_out["yhat"] = yhats
    return df_out

def pam_sets_by_block(df, block_col="syntenic_block", pam_col="protospacerID"):
    tmp = df.dropna(subset=[block_col]).explode(block_col)
    return tmp.groupby(block_col)[pam_col].apply(set).to_dict()

def compare_pams_by_block(ref_df, qry_df, tol, block_col="syntenic_block"):
    ref_sets = pam_sets_by_block(ref_df, block_col=block_col)
    qry_sets = pam_sets_by_block(qry_df, block_col=block_col)

    ref_mids = dict(zip(ref_df['protospacerID'], ref_df['midpoint']))
    qry_mids = dict(zip(qry_df['protospacerID'], qry_df['midpoint']))
    ref_yhats = dict(zip(ref_df['protospacerID'], ref_df['yhat']))

    rows = []
    all_blocks = set(ref_sets.keys()) | set(qry_sets.keys()) 

    for block in all_blocks:
        if pd.isna(block): continue
            
        ref_pams = ref_sets.get(block, set())
        qry_pams = qry_sets.get(block, set())
        
        shared = set()
        for pam in ref_pams & qry_pams:
            y_actual = qry_mids[pam]
            y_proj = ref_yhats[pam]
            
            if y_proj is not None:
                if abs(y_proj - y_actual) <= tol:
                    shared.add(pam)
            else:
                shared.add(pam) 

        ref_only = ref_pams - shared
        qry_only = qry_pams - shared

        rows.append({
            "syntenic_pams": sorted(shared),
            "ref_only_pams": sorted(ref_only),
            "qry_only_pams": sorted(qry_only),
        })
    return pd.DataFrame(rows)

def plot_tolerance_venns(ref_df, qry_df, output_dir):
    tolerances = [100, 300, 1000, 10000, 100000, 1000000]
    
    fig, axes = plt.subplots(2, 3, figsize=(20, 20))
    axes = axes.flatten()
    
    for i, tol in enumerate(tolerances):
        print(f"[INFO] Calculating overlaps for tolerance: ±{tol:,} bp")
        
        summary = compare_pams_by_block(ref_df, qry_df, tol)
        
        shared_set = set(pam for sublist in summary["syntenic_pams"] for pam in sublist)
        ref_only_set = set(pam for sublist in summary["ref_only_pams"] for pam in sublist)
        qry_only_set = set(pam for sublist in summary["qry_only_pams"] for pam in sublist)

        true_ref_only = len(ref_only_set - shared_set)
        true_qry_only = len(qry_only_set - shared_set)
        true_shared = len(shared_set)
        
        v = venn2(
            subsets=(true_ref_only, true_qry_only, true_shared),
            set_labels=("Reference PAMs", "Query PAMs"),
            ax=axes[i]
        )
        
        total_pams = true_ref_only + true_qry_only + true_shared
        subset_data = {'10': true_ref_only, '01': true_qry_only, '11': true_shared}
        
        for subset_id, value in subset_data.items():
            label = v.get_label_by_id(subset_id)
            if label:
                pct = (value / total_pams) * 100
                label.set_text(f"{value:,}\n({pct:.1f}%)")
                label.set_fontsize(11)
                
        axes[i].set_title(f"Tolerance: ±{tol:,} bp", fontsize=14, pad=10)

    plt.suptitle("Block-Aware PAM Conservation Across Spatial Tolerances", fontsize=18, y=0.95)
    plt.tight_layout(rect=[0, 0, 1, 0.93]) 
    
    output = output_dir / "multi_tolerance_venn.png"
    plt.savefig(output)
    plt.close()

def main():
    args = parse_args()
    print("[INFO] Loading inputs...")
    ref_cfd = load_cfd_table(args.ref_cfd)
    qry_cfd = load_cfd_table(args.query_cfd)
    _, primary_alns = parse(args.delta)

    print("[INFO] Resolving synteny...")
    resolve_synteny(primary_alignments=primary_alns, breakpoint_map=None)

    reference_trees = {aln.reference: aln.reference_synteny_tree for aln in primary_alns.values()}
    query_trees = {aln.query: aln.query_synteny_tree for aln in primary_alns.values()}

    print("[INFO] Assigning synteny and calculating projections...")
    ref_df = assign_synteny(ref_cfd, reference_trees, mode="reference")
    qry_df = assign_synteny(qry_cfd, query_trees, mode="query")

    ref_df = ref_df.explode(['syntenic_block', 'yhat'], ignore_index=True)
    ref_df = ref_df.dropna(subset=['yhat'])

    print("[INFO] Plotting multi-tolerance Venn diagrams...")
    plot_tolerance_venns(ref_df, qry_df, args.figures)
    print(f"[INFO] Success. Check the {args.figures} directory.")

if __name__ == "__main__":
    main()