''' 
AI is not performing well at adapting this script. The results are far too verbose and complex for what should be minimal
Free hand implementation like the good old days before the cyborgs. 
'''

from svmu2.orchestration.parse import parse
from svmu2.orchestration.synteny import resolve_synteny
from svmu2.models.line import build_alignment_primitives ## temporary
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import numpy as np
import seaborn as sns
import argparse
from pathlib import Path
import pandas as pd
import os
import gc


def parse_args():
    parser = argparse.ArgumentParser(
        description="Block-first PAM/protospacer orthology classifier using svmu2 synteny.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--ref-cfd", required=True, type=Path, help="Reference/ISO1 CFD CSV")
    parser.add_argument("--query-cfd", required=True, type=Path, help="Query/BL54591 CFD CSV")
    parser.add_argument("--delta", required=True, type=Path, help="nucmer delta file")
    parser.add_argument("--out", required=True, type=Path, help="Output orthology summary CSV")
    parser.add_argument("--figures", required = True, type=Path, help="Output dir to store the figures in")
    parser.add_argument("--diagnostics-out", type=Path, default=None, help="Optional CSV for no-block/ambiguous-block hit diagnostics")
    parser.add_argument("--tol",type=int,default=1000,help="Maximum bp distance between projected ref midpoint and query midpoint")
    parser.add_argument("--include-query-only",action="store_true",help="Also emit diagnostic rows for query hits in a block lacking same-ID ref hits")

    args = parser.parse_args()

    path_obj = Path(args.figures)

    if not path_obj.exists():
        path_obj.mkdir(parents=True, exist_ok=True)
        print(f"Created directory: {path_obj}")
    

    return args

REQUIRED_CFD_COLUMNS = {
    "protospacerID",
    "ref_chr",
    "ref_start",
    "ref_end",
}

def build_projection_models(primary_alns):
    print("[INFO] Extracting linear projection models from syntenic blocks...")
    block_models = {}

    for aln in primary_alns.values():
        tree = aln.reference_synteny_tree
        if not tree:
            continue
            
        for interval in tree:
            block = interval.data
            b_id = block_id(block)
            
            if b_id not in block_models:
                # Directly grab the pre-calculated attributes
                block_models[b_id] = {
                    'm': block.slope, 
                    'b': block.y_intercept
                }

    return block_models

def load_cfd_table(path: str) -> pd.DataFrame:
    # 1. Define types upfront to bypass inference and memory copies
    dtypes = {
        "protospacerID": str,
        "ref_chr": str,
        "ref_start": int, # Or "Int32" to handle NaN
        "ref_end": int,
    }
    
    df = pd.read_csv(path, dtype=dtypes)

    missing = REQUIRED_CFD_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"{path} CFD table is missing columns: {sorted(missing)}")
    
    # Vectorized midpoint calculation across the entire dataframe
    df["midpoint"] = (df["ref_start"] + df["ref_end"]) / 2.0

    # 2. Vectorized string ops (already good, just kept clean)
    if "name" not in df.columns:
        df["name"] = df["protospacerID"] + "_row" + df.index.astype(str)

    if "CFD" not in df.columns:
        df["CFD"] = pd.NA

    return df

def block_id(block):
    return f"{block.reference}_{block.query}_{block.index}"

def assign_synteny(df, trees):
    # Depending on memory constraints, we might want to drop df.copy() 
    # here and just assign the new lists directly to the incoming df.
    df_out = df.copy() 

    block_ids = []
    hit_counts = []

    for row in df.itertuples(index=False):
        # Direct attribute access instead of getattr()
        tree = trees.get(row.ref_chr)

        if tree is None:
            block_ids.append(None)
            hit_counts.append(0)
            continue

        hits = tree[row.ref_start:row.ref_end]
        hit_counts.append(len(hits))

        if not hits:
            block_ids.append(None)
        elif len(hits) == 1:
            block_ids.append(block_id(next(iter(hits)).data))
        else:
            block_ids.append([block_id(hit.data) for hit in hits])

    df_out["syntenic_block"] = block_ids
    df_out["n_hits"] = hit_counts

    return df_out

def plot(primary):
    targets = primary.values()

    for aln in targets:
        primitives = build_alignment_primitives(aln)
        xlabel = aln.reference
        ylabel = aln.query
        from svmu2.visualization.renderers import render_plotly
        title = f"{aln.reference}_{aln.query}"
        fig = render_plotly(primitives, title=title)
        fig.show()

def pam_sets_by_block(df, block_col="syntenic_block", pam_col="protospacerID"):
    tmp = (
        df.dropna(subset=[block_col])
          .explode(block_col)
    )

    return (
        tmp.groupby(block_col)[pam_col]
           .apply(set)
           .to_dict()
    )

def compare_pams_by_block(ref_df, qry_df, block_models, tol, block_col="syntenic_block"):
    ref_sets = pam_sets_by_block(ref_df, block_col=block_col)
    qry_sets = pam_sets_by_block(qry_df, block_col=block_col)

    # Fast O(1) lookups using the pre-calculated vectorized column
    ref_mids = dict(zip(ref_df['protospacerID'], ref_df['midpoint']))
    qry_mids = dict(zip(qry_df['protospacerID'], qry_df['midpoint']))

    rows = []
    all_blocks = set(ref_sets) | set(qry_sets) 

    for block in all_blocks:
        ref_pams = ref_sets.get(block, set())
        qry_pams = qry_sets.get(block, set())

        potential_shared = ref_pams & qry_pams
        shared = set()
        
        model = block_models.get(block)
        for pam in potential_shared:
            if model and model['m'] != 0:
                m, b = model['m'], model['b']
                # Project query coordinate onto reference space
                # THIS HERE NEEDS SOME ATTENTION
                ref_proj = (qry_mids[pam] - b) / m
                
                if abs(ref_mids[pam] - ref_proj) <= tol:
                    shared.add(pam)
            else:
                shared.add(pam) 

        ref_only = ref_pams - shared
        qry_only = qry_pams - shared

        rows.append({
            "syntenic_block": block,
            "n_ref_pams": len(ref_pams),
            "n_qry_pams": len(qry_pams),
            "n_syntenic_pams": len(shared),
            "n_ref_only_pams": len(ref_only),
            "n_qry_only_pams": len(qry_only),
            "syntenic_pams": sorted(shared),
            "ref_only_pams": sorted(ref_only),
            "qry_only_pams": sorted(qry_only),
        })

    return pd.DataFrame(rows)

def merge_by_block_to_disk(ref_df, qry_df, temp_csv, block_models, tol, merge_col="protospacerID"):
    print("[INFO] Setting up block-wise merge with spatial tolerance...")
    
    ref_df = ref_df.dropna(subset=['syntenic_block']).explode('syntenic_block')
    qry_df = qry_df.dropna(subset=['syntenic_block']).explode('syntenic_block')

    ref_grouped = ref_df.groupby('syntenic_block')
    qry_grouped = qry_df.groupby('syntenic_block')
    
    first_chunk = True
    
    for block_name, ref_sub in ref_grouped:
        if block_name in qry_grouped.groups:
            qry_sub = qry_grouped.get_group(block_name)
            merged_sub = ref_sub.merge(qry_sub, on=merge_col, suffixes=('_ref', '_qry'))
            
            # Apply the projection tolerance filter using pre-calculated midpoints
            if not merged_sub.empty and block_name in block_models:
                m, b = block_models[block_name]['m'], block_models[block_name]['b']
                if m != 0:
                    ref_proj = (merged_sub['midpoint_qry'] - b) / m
                    dist = (merged_sub['midpoint_ref'] - ref_proj).abs()
                    
                    merged_sub = merged_sub[dist <= tol]
            
            if merged_sub.empty:
                continue

            merged_sub.to_csv(
                temp_csv, 
                mode='w' if first_chunk else 'a', 
                header=first_chunk, 
                index=False
            )
            first_chunk = False

    print(f"[INFO] Block-wise merge complete. Temporary data written to {temp_csv}")

def plot_cfd_shifts(temp_csv, output_dir):

    print("[INFO] Loading merged data for vectorized math...")
    # Load the merged CSV. This is safe because it only contains the shared rows now.
    conserved = pd.read_csv(temp_csv)
    
    # Your vectorized math
    conserved['delta_CFD'] = conserved['CFD_qry'] - conserved['CFD_ref']
    
    print("[INFO] Rendering CFD plot...")
    plt.figure(figsize=(8, 6))
    
    sns.histplot(
        data=conserved, 
        x='delta_CFD', 
        bins=50, 
        color='#7b85ba', 
        edgecolor='black',
        linewidth=0.5,     # Thins the border so it doesn't overpower the color
        shrink=0.9,        # Creates a 10% gap between each bar
        log_scale=(False, True) 
    )
    
    plt.title("Distribution of CFD Shifts for Conserved PAMs", pad=15)
    plt.xlabel(r'$\Delta$CFD (BL54591 - ISO1)', fontsize=12)
    plt.ylabel('Count (Log Scale)', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xlim(-1.05, 1.05)
    
    plt.tight_layout()
    output = output_dir / "cfd_shift.png"
    plt.savefig(output)
    plt.close()
    
    # Clean up the temporary file so we don't leave trash on the drive
    if temp_csv.exists():
        pass
        #os.remove(temp_csv)
    
    del conserved
    gc.collect()

def plot_mismatch_distribution(ref_df, qry_df, output_dir):
    # Count occurrences of each NM (mismatch) value
    ref_counts = ref_df['mismatches'].value_counts().sort_index()
    qry_counts = qry_df['mismatches'].value_counts().sort_index()
    
    # Align the indices in case one genome is missing a specific NM category
    all_nms = sorted(list(set(ref_counts.index) | set(qry_counts.index)))
    
    # Extract aligned values, filling missing categories with 0
    ref_vals = [ref_counts.get(nm, 0) for nm in all_nms]
    qry_vals = [qry_counts.get(nm, 0) for nm in all_nms]
    
    # Setup for grouped bars
    x = np.arange(len(all_nms))
    width = 0.35  
    
    plt.figure(figsize=(8, 6))
    
    # Plotting
    plt.bar(x - width/2, ref_vals, width, label='ISO1 (Reference)', color='#5ab48c', edgecolor='black')
    plt.bar(x + width/2, qry_vals, width, label='BL54591 (Query)', color='#f28e62', edgecolor='black')
    
    # Formatting
    plt.title("Mismatch Counts from SAM (per Perfect PAM Search)", pad=15)
    plt.xlabel('Number of Mismatches (NM)', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.xticks(x, all_nms)
    
    # Dress up the y-axis with comma formatting for thousands/millions
    ax = plt.gca()
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
    
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    output = output_dir / "mismatch_distribution.png"
    plt.savefig(output)

def plot_venn(summary, output_dir):
    # Plot the aggregated, block-aware Venn diagram
    # Bumped figure size slightly to accommodate the extra text lines

    # Aggregate absolute counts across all syntenic blocks
    total_ref_only = summary["n_ref_only_pams"].sum()
    total_qry_only = summary["n_qry_only_pams"].sum()
    total_shared = summary["n_syntenic_pams"].sum()
    plt.figure(figsize=(8, 8)) 
    
    v = venn2(
        subsets=(total_ref_only, total_qry_only, total_shared),
        set_labels=("Reference PAMs", "Query PAMs")
    )

    # Calculate total for percentages
    total_pams = total_ref_only + total_qry_only + total_shared

    # Map the venn subset IDs to their calculated values
    # '10' = Left only, '01' = Right only, '11' = Intersection
    subset_data = {
        '10': total_ref_only,
        '01': total_qry_only,
        '11': total_shared
    }

    # Iterate through the subsets and format the text
    for subset_id, value in subset_data.items():
        label = v.get_label_by_id(subset_id)
        if label:
            pct = (value / total_pams) * 100
            # f-string formatting: 
            # {value:,} adds the thousands separator
            # {pct:.1f} formats the float to 1 decimal place
            label.set_text(f"{value:,}\n({pct:.1f}%)")
            label.set_fontsize(11)

    plt.title("Block-Aware PAM Conservation\n(Aggregated across all syntenic blocks)")
    plt.tight_layout()
    output = output_dir / "venn_diagram.png"
    plt.savefig(output)

def main():
    args = parse_args()
    ref_cfd = load_cfd_table(args.ref_cfd)
    qry_cfd = load_cfd_table(args.query_cfd)
    all_alns, primary_alns = parse(args.delta)

    #plot(primary=primary_alns)
    resolve_synteny(primary_alignments=primary_alns, breakpoint_map=None)
    # aln objs now have .reference_synteny_tree and .query_synteny_tree attributes

    blocks_dictionary = build_projection_models(primary_alns)

    reference_trees = {aln.reference: aln.reference_synteny_tree for aln in primary_alns.values()}
    query_trees = {aln.query: aln.query_synteny_tree for aln in primary_alns.values()}

    ref_df = assign_synteny(ref_cfd, reference_trees) ## assign synteny information to each 23mer
    qry_df = assign_synteny(qry_cfd, query_trees)

    summary = compare_pams_by_block(ref_df, qry_df, blocks_dictionary, args.tol)
    summary.to_csv(args.out)

    # Diagnostic tally for PAMs falling completely outside blocks
    ref_no_synteny = ref_df["syntenic_block"].isna().sum()
    qry_no_synteny = qry_df["syntenic_block"].isna().sum()
    print(f"Diagnostics — No synteny detected: {ref_no_synteny} Ref PAMs, {qry_no_synteny} Query PAMs.")

    plot_venn(summary, output_dir = args.figures)
    plot_mismatch_distribution(ref_df, qry_df, output_dir = args.figures)

    temp_csv = args.out.parent / "temp_merged_pams.csv"
    
    # Execute your chunked merge strategy
    merge_by_block_to_disk(ref_df, qry_df, temp_csv, blocks_dictionary, args.tol)
    
    # Nuke the original heavy dataframes from RAM -- 64 Gb system ram does not handle drosophila sized comparisons if we aren't careful
    del ref_df
    del qry_df
    gc.collect()

    plot_cfd_shifts(temp_csv, output_dir=args.figures)

if __name__ == "__main__":
    main()