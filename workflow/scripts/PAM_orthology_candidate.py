from svmu2.orchestration.parse import parse
from svmu2.orchestration.synteny import resolve_synteny
from svmu2.models.line import build_alignment_primitives
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

def load_cfd_table(path: str) -> pd.DataFrame:
    dtypes = {
        "protospacerID": str,
        "ref_chr": str,
        "ref_start": int, 
        "ref_end": int,
    }
    
    df = pd.read_csv(path, dtype=dtypes)

    missing = REQUIRED_CFD_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"{path} CFD table is missing columns: {sorted(missing)}")
    
    df["midpoint"] = (df["ref_start"] + df["ref_end"]) / 2.0

    if "name" not in df.columns:
        df["name"] = df["protospacerID"] + "_row" + df.index.astype(str)

    if "CFD" not in df.columns:
        df["CFD"] = pd.NA

    return df

def block_id(block):
    return f"{block.reference}_{block.query}_{block.index}"

def assign_synteny(df, trees, mode="reference"):
    """
    Queries the IntervalTree to assign syntenic blocks to PAMs.
    If mode='reference', calculates the projected coordinate (yhat) immediately
    and stores the slope (m) and intercept (b) for debugging.
    """
    df_out = df.copy() 

    block_ids = []
    hit_counts = []
    
    if mode == "reference":
        yhats = []
        slopes = []
        intercepts = []
    else:
        yhats = slopes = intercepts = None

    for row in df.itertuples(index=False):
        tree = trees.get(row.ref_chr)
        
        # We use midpoint for both the tree lookup and the projection math
        x = row.midpoint 

        if tree is None:
            block_ids.append(None)
            hit_counts.append(0)
            if mode == "reference": 
                yhats.append(None)
                slopes.append(None)
                intercepts.append(None)
            continue

        # Use .at() to query the exact midpoint, preventing boundary overlaps
        hits = tree.at(x)
        hit_counts.append(len(hits))

        if not hits:
            block_ids.append(None)
            if mode == "reference": 
                yhats.append(None)
                slopes.append(None)
                intercepts.append(None)
            
        elif len(hits) == 1:
            hit_data = next(iter(hits)).data
            block_ids.append(hit_data.index) 
            
            if mode == "reference":
                yhats.append((hit_data.slope * x) + hit_data.y_intercept)
                slopes.append(hit_data.slope)
                intercepts.append(hit_data.y_intercept)
            
        else:
            block_ids.append([hit.data.index for hit in hits])
            
            if mode == "reference":
                yhats.append([(hit.data.slope * x) + hit.data.y_intercept for hit in hits])
                slopes.append([hit.data.slope for hit in hits])
                intercepts.append([hit.data.y_intercept for hit in hits])

    # Assign base columns
    df_out["syntenic_block"] = block_ids
    df_out["n_hits"] = hit_counts
    
    # Only assign math columns if we are in reference mode
    if mode == "reference":
        df_out["yhat"] = yhats
        df_out["m"] = slopes
        df_out["b"] = intercepts

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

def compare_pams_by_block(ref_df, qry_df, tol, block_col="syntenic_block"):
    ref_sets = pam_sets_by_block(ref_df, block_col=block_col)
    qry_sets = pam_sets_by_block(qry_df, block_col=block_col)

    ref_mids = dict(zip(ref_df['protospacerID'], ref_df['midpoint']))
    qry_mids = dict(zip(qry_df['protospacerID'], qry_df['midpoint']))
    ref_yhats = dict(zip(ref_df['protospacerID'], ref_df['yhat']))

    rows = []
    plot_data = {} 
    all_blocks = set(ref_sets.keys()) | set(qry_sets.keys()) 

    for block in all_blocks:
        if pd.isna(block):
            continue
            
        ref_pams = ref_sets.get(block, set())
        qry_pams = qry_sets.get(block, set())

        potential_shared = ref_pams & qry_pams
        shared = set()
        
        plot_data[block] = {'x': [], 'y_actual': [], 'y_proj': []}

        for pam in potential_shared:
            x = ref_mids[pam]
            y_actual = qry_mids[pam]
            y_proj = ref_yhats[pam]
            
            if y_proj is not None:
                plot_data[block]['x'].append(x)
                plot_data[block]['y_actual'].append(y_actual)
                plot_data[block]['y_proj'].append(y_proj)

                if abs(y_proj - y_actual) <= tol:
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

    return pd.DataFrame(rows), plot_data

def plot_synteny_projections(plot_data, tol, output_dir):
    fig, ax = plt.subplots(figsize=(10, 8))
    
    colors = plt.get_cmap('tab10', max(10, len(plot_data)))
    
    for idx, (block, data) in enumerate(plot_data.items()):
        if not data['x']:
            continue
            
        xs = data['x']
        ys_actual = data['y_actual']
        ys_proj = data['y_proj']
        
        ax.scatter(xs, ys_actual, color=colors(idx), alpha=0.7, label=f'Block {block} Actual', zorder=3)
        
        sorted_indices = sorted(range(len(xs)), key=lambda k: xs[k])
        xs_sorted = [xs[i] for i in sorted_indices]
        ys_proj_sorted = [ys_proj[i] for i in sorted_indices]
        
        ax.plot(xs_sorted, ys_proj_sorted, color=colors(idx), linestyle='-', linewidth=2, zorder=2)
        ax.plot(xs_sorted, [y + tol for y in ys_proj_sorted], color=colors(idx), linestyle='--', alpha=0.5)
        ax.plot(xs_sorted, [y - tol for y in ys_proj_sorted], color=colors(idx), linestyle='--', alpha=0.5)

    ax.set_xlabel('Reference Midpoint (Absolute bp)')
    ax.set_ylabel('Query Midpoint (Absolute bp)')
    ax.set_title(f'PAM Coordinate Projection (Tolerance = $\pm${tol} bp)')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, linestyle=':', alpha=0.6)
    
    plt.tight_layout()
    output = output_dir / "synteny_projections.png"
    plt.savefig(output)
    plt.close()

def merge_by_block_to_disk(ref_df, qry_df, temp_csv, tol, merge_col="protospacerID"):
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
            
            # Using the pre-calculated inline yhat projection
            if not merged_sub.empty:
                dist = (merged_sub['midpoint_qry'] - merged_sub['yhat']).abs()
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
    conserved = pd.read_csv(temp_csv)
    
    conserved['delta_CFD'] = conserved['CFD_qry'] - conserved['CFD_ref']
    
    print("[INFO] Rendering CFD plot...")
    plt.figure(figsize=(8, 6))
    
    sns.histplot(
        data=conserved, 
        x='delta_CFD', 
        bins=50, 
        color='#7b85ba', 
        edgecolor='black',
        linewidth=0.5,
        shrink=0.9,
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
    
    if temp_csv.exists():
        pass 
    
    del conserved
    gc.collect()

def plot_tolerance_venns(ref_df, qry_df, output_dir):
    tolerances = [1000, 10000, 100000, 1000000]
    
    # Create a 2x2 grid for the plots
    fig, axes = plt.subplots(2, 2, figsize=(20, 20))
    axes = axes.flatten()
    
    for i, tol in enumerate(tolerances):
        print(f"[INFO] Calculating overlaps for tolerance: ±{tol:,} bp")
        
        # We can reuse your existing comparison function here
        summary, _ = compare_pams_by_block(ref_df, qry_df, tol)
        
        # Apply the strict deduplication logic
        shared_set = set(pam for sublist in summary["syntenic_pams"] for pam in sublist)
        ref_only_set = set(pam for sublist in summary["ref_only_pams"] for pam in sublist)
        qry_only_set = set(pam for sublist in summary["qry_only_pams"] for pam in sublist)

        true_ref_only = len(ref_only_set - shared_set)
        true_qry_only = len(qry_only_set - shared_set)
        true_shared = len(shared_set)
        
        # Plot to the specific subplot axis
        v = venn2(
            subsets=(true_ref_only, true_qry_only, true_shared),
            set_labels=("Reference PAMs", "Query PAMs"),
            ax=axes[i]
        )
        
        total_pams = true_ref_only + true_qry_only + true_shared
        subset_data = {'10': true_ref_only, '01': true_qry_only, '11': true_shared}
        
        # Dress up the labels with comma formatting and percentages
        for subset_id, value in subset_data.items():
            label = v.get_label_by_id(subset_id)
            if label:
                pct = (value / total_pams) * 100
                label.set_text(f"{value:,}\n({pct:.1f}%)")
                label.set_fontsize(11)
                
        # Title each subplot with its specific tolerance
        axes[i].set_title(f"Tolerance: ±{tol:,} bp", fontsize=14, pad=10)

    # Add a master title to the entire figure
    plt.suptitle("Block-Aware PAM Conservation Across Spatial Tolerances", fontsize=18, y=0.95)
    
    plt.tight_layout(rect=[0, 0, 1, 0.93]) # Adjust layout to fit the suptitle
    output = output_dir / "multi_tolerance_venn.png"
    plt.savefig(output)
    plt.close()

def plot_mismatch_distribution(ref_df, qry_df, output_dir):
    ref_counts = ref_df['mismatches'].value_counts().sort_index()
    qry_counts = qry_df['mismatches'].value_counts().sort_index()
    
    all_nms = sorted(list(set(ref_counts.index) | set(qry_counts.index)))
    
    ref_vals = [ref_counts.get(nm, 0) for nm in all_nms]
    qry_vals = [qry_counts.get(nm, 0) for nm in all_nms]
    
    x = np.arange(len(all_nms))
    width = 0.35  
    
    plt.figure(figsize=(8, 6))
    
    plt.bar(x - width/2, ref_vals, width, label='ISO1 (Reference)', color='#5ab48c', edgecolor='black')
    plt.bar(x + width/2, qry_vals, width, label='BL54591 (Query)', color='#f28e62', edgecolor='black')
    
    plt.title("Mismatch Counts from SAM (per Perfect PAM Search)", pad=15)
    plt.xlabel('Number of Mismatches (NM)', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.xticks(x, all_nms)
    
    ax = plt.gca()
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
    
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    output = output_dir / "mismatch_distribution.png"
    plt.savefig(output)
    plt.close()

def plot_venn(summary, output_dir):
    # 1. Flatten the lists of IDs into global sets to eliminate double-counting
    shared_set = set(pam for sublist in summary["syntenic_pams"] for pam in sublist)
    ref_only_set = set(pam for sublist in summary["ref_only_pams"] for pam in sublist)
    qry_only_set = set(pam for sublist in summary["qry_only_pams"] for pam in sublist)

    # 2. Enforce global categorization
    # If a PAM was shared in one block, it cannot be considered "only" in another
    true_ref_only = len(ref_only_set - shared_set)
    true_qry_only = len(qry_only_set - shared_set)
    true_shared = len(shared_set)

    plt.figure(figsize=(8, 8)) 
    
    v = venn2(
        subsets=(true_ref_only, true_qry_only, true_shared),
        set_labels=("Reference PAMs", "Query PAMs")
    )

    total_pams = true_ref_only + true_qry_only + true_shared

    subset_data = {
        '10': true_ref_only,
        '01': true_qry_only,
        '11': true_shared
    }

    for subset_id, value in subset_data.items():
        label = v.get_label_by_id(subset_id)
        if label:
            pct = (value / total_pams) * 100
            label.set_text(f"{value:,}\n({pct:.1f}%)")
            label.set_fontsize(11)

    plt.title("Unique PAM Conservation\n(Deduplicated across all syntenic blocks)")
    plt.tight_layout()
    output = output_dir / "venn_diagram.png"
    plt.savefig(output)
    plt.close()

def main():
    args = parse_args()
    ref_cfd = load_cfd_table(args.ref_cfd)
    qry_cfd = load_cfd_table(args.query_cfd)
    all_alns, primary_alns = parse(args.delta)

    resolve_synteny(primary_alignments=primary_alns, breakpoint_map=None)

    reference_trees = {aln.reference: aln.reference_synteny_tree for aln in primary_alns.values()}
    query_trees = {aln.query: aln.query_synteny_tree for aln in primary_alns.values()}

    ref_df = assign_synteny(ref_cfd, reference_trees, mode="reference")
    qry_df = assign_synteny(qry_cfd, query_trees, mode="query")

    # 1 & 2: Reassign the dataframe and specify the columns to explode
    ref_df = ref_df = ref_df.explode(['syntenic_block', 'yhat', 'm', 'b'], ignore_index=True)

    # 3: Drop PAMs that didn't align to a syntenic block
    plot_df = ref_df.dropna(subset=['yhat'])
    
    # 4: Isolate chr2L
    plot_df = plot_df[plot_df['ref_chr'] == '2L']
    sub = plot_df[(plot_df['midpoint'] > 9780000) & (plot_df['midpoint'] < 9880000)]
    sub.to_csv('subset_debugging.csv')

    # --- THE FIX ---
    # Filter for strict 1:1 orthology (no overlapping duplicate blocks)
    strict_orthologs = plot_df[plot_df['n_hits'] == 1]
    
    # Optional: Look at the paralogs/duplications to prove the theory!
    duplicates = plot_df[plot_df['n_hits'] > 1]
    duplicates.to_csv('debug_duplicates.csv')

    # Plot the strict 1:1 alignments in blue
    plt.scatter(strict_orthologs["midpoint"], strict_orthologs["yhat"], label="1:1 Orthologs", color="blue")
    
    # Plot the overlapping/duplicated regions in red
    plt.scatter(duplicates["midpoint"], duplicates["yhat"], label="Duplications (n_hits > 1)", color="red", alpha=0.5)
    
    plt.legend()
    plt.show()

    summary, plot_data = compare_pams_by_block(ref_df, qry_df, args.tol)
    summary.to_csv(args.out)

    # Save the absolute coordinate projection scatterplot to the figures directory
    plot_synteny_projections(plot_data, tol=args.tol, output_dir=args.figures)

    ref_no_synteny = ref_df["syntenic_block"].isna().sum()
    qry_no_synteny = qry_df["syntenic_block"].isna().sum()
    print(f"Diagnostics — No synteny detected: {ref_no_synteny} Ref PAMs, {qry_no_synteny} Query PAMs.")

    plot_venn(summary, output_dir = args.figures)
    plot_tolerance_venns(ref_df, qry_df, output_dir=args.figures)
    plot_mismatch_distribution(ref_df, qry_df, output_dir = args.figures)

    temp_csv = args.out.parent / "temp_merged_pams.csv"
    
    merge_by_block_to_disk(ref_df, qry_df, temp_csv, args.tol)
    
    del ref_df
    del qry_df
    gc.collect()

    plot_cfd_shifts(temp_csv, output_dir=args.figures)

if __name__ == "__main__":
    main()