import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import argparse
from pathlib import Path
import numpy as np
import gc

# Note: Since the discovery logic was separated, ensure you import your 
# synteny functions here if you are still dynamically calling compare_pams_by_block.
# from synteny_logic import compare_pams_by_block 

def parse_args():
    parser = argparse.ArgumentParser(description="Plotter for PAM synteny diagnostics.")
    parser.add_argument("--data-dir", required=True, type=Path, help="Directory containing the mapper CSV outputs")
    parser.add_argument("--figures", required=True, type=Path, help="Directory to save plots")
    parser.add_argument("--ref-cfd", required=True, type=Path, help="Reference/ISO1 CFD CSV")
    parser.add_argument("--query-cfd", required=True, type=Path, help="Query/BL54591 CFD CSV")
    
    args = parser.parse_args()
    args.figures.mkdir(parents=True, exist_ok=True)
    return args

def plot_chromosome_traces(data_dir, fig_dir):
    print("[INFO] Loading mapping data...")
    conserved = pd.read_csv(data_dir / "conserved_pams.csv")
    failed = pd.read_csv(data_dir / "failed_tolerance_pams.csv")

    # Get all unique chromosomes present in the failed data
    chromosomes = pd.concat([failed['ref_chr_ref']]).dropna().unique()

    for chrom in chromosomes:
        print(f"[INFO] Rendering standard static trace for {chrom}...")
        
        # Create a side-by-side layout for the trace and the histogram
        fig, (ax_trace, ax_hist) = plt.subplots(1, 2, figsize=(18, 8))
        
        sub = conserved[conserved['ref_chr_ref'] == chrom].copy()
        sub['difference'] = sub['yhat'] - sub['midpoint_qry']
        
        # --- Synteny Trace Scatter Plot ---
        # Kept alpha low and marker size (s) small to help visualize dense regions
        ax_trace.scatter(sub['midpoint_ref'], sub['midpoint_qry'], color='red', label='Mapped Midpoint (Truth)', alpha=0.6, s=10)
        ax_trace.scatter(sub['midpoint_ref'], sub['yhat'], color='blue', label='yhat (Interpolated)', alpha=0.6, s=10)
        
        # Formatting - Keeping absolute terms
        ax_trace.set_title(f"Synteny Trace: Chromosome {chrom}", pad=15)
        ax_trace.set_xlabel("Reference Midpoint (Absolute bp)", fontsize=12)
        ax_trace.set_ylabel("Query Midpoint (Absolute bp)", fontsize=12)
        
        ax_trace.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax_trace.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda y, p: format(int(y), ',')))
        
        ax_trace.legend(loc='upper left')
        ax_trace.grid(True, linestyle='--', alpha=0.4)
        
        # --- Difference Histogram ---
        ax_hist.hist(sub['difference'].dropna(), bins=100, log=True, color='#7b85ba', edgecolor='black', alpha=0.8)
        
        # Formatting - Keeping absolute terms
        ax_hist.set_title(f"Tolerance Diagnostics (yhat - truth): Chromosome {chrom}", pad=15)
        ax_hist.set_xlabel("Difference (Absolute bp)", fontsize=12)
        ax_hist.set_ylabel("Count (Log Scale)", fontsize=12)
        ax_hist.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax_hist.grid(True, linestyle='--', alpha=0.4)
        
        plt.tight_layout()
        
        # Save as a high-resolution PNG
        output_file = fig_dir / f"trace_{chrom}.png"
        plt.savefig(output_file, dpi=300)
        plt.close(fig)
        
        print(f"[INFO] Saved static trace to {output_file}")
        
def plot_cfd_shifts(data_dir, output_dir):
    print("[INFO] Loading merged data for vectorized math...")
    conserved = pd.read_csv(data_dir / "conserved_pams.csv")
    
    conserved['delta_CFD'] = conserved['CFD_qry'] - conserved['CFD_ref']
    
    print("[INFO] Rendering CFD plot...")
    plt.figure(figsize=(8, 6))
    
    # Implemented with strict Matplotlib, maintaining absolute counts and log scaling
    plt.hist(
        conserved['delta_CFD'].dropna(), 
        bins=50, 
        color='#7b85ba', 
        edgecolor='black',
        linewidth=0.5,
        rwidth=0.9,
        log=True
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
    
    del conserved
    gc.collect()

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

def main():
    args = parse_args()
    
    # Base synteny traces
    plot_chromosome_traces(args.data_dir, args.figures)
    print(f"[INFO] Diagnostic traces saved to {args.figures}")

    plot_cfd_shifts(args.data_dir, args.figures)
    print(f"[INFO] CFD shifts plot saved to {args.figures}")

if __name__ == "__main__":
    main()