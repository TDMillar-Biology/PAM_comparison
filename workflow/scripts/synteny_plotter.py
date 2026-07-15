import matplotlib.pyplot as plt
import pandas as pd
import argparse
from pathlib import Path
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="Plotter for PAM synteny diagnostics.")
    parser.add_argument("--data-dir", required=True, type=Path, help="Directory containing the mapper CSV outputs")
    parser.add_argument("--figures", required=True, type=Path, help="Directory to save plots")
    
    args = parser.parse_args()
    args.figures.mkdir(parents=True, exist_ok=True)
    return args

def plot_chromosome_traces(data_dir, fig_dir):
    print("[INFO] Loading mapping data...")
    conserved = pd.read_csv(data_dir / "conserved_pams.csv")
    failed = pd.read_csv(data_dir / "failed_tolerance_pams.csv")
    unmapped_ref = pd.read_csv(data_dir / "unmapped_ref.csv")
    unmapped_qry = pd.read_csv(data_dir / "unmapped_qry.csv")

    # Get all unique chromosomes present in the reference data
    chromosomes = pd.concat([conserved['ref_chr_ref'], unmapped_ref['ref_chr']]).dropna().unique()

    for chrom in chromosomes:
        print(f"[INFO] Rendering trace for {chrom}...")
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # 1. Plot the unmapped points along the axes
        chr_unmapped_ref = unmapped_ref[unmapped_ref['ref_chr'] == chrom]
        chr_unmapped_qry = unmapped_qry[unmapped_qry['ref_chr'] == chrom] # Assuming CFD has ref_chr populated for query
        
        ax.scatter(chr_unmapped_ref['midpoint'], np.zeros(len(chr_unmapped_ref)), 
                   color='grey', marker='|', alpha=0.5, label='Unmapped Ref PAMs (X-axis)')
        
        ax.scatter(np.zeros(len(chr_unmapped_qry)), chr_unmapped_qry['midpoint'], 
                   color='orange', marker='_', alpha=0.5, label='Unmapped Query PAMs (Y-axis)')

        # 2. Plot the valid, within-tolerance traces
        chr_conserved = conserved[conserved['ref_chr_ref'] == chrom]
        ax.scatter(chr_conserved['midpoint_ref'], chr_conserved['yhat'], 
                   color='#7b85ba', s=2, alpha=0.8, label='Conserved Trace (Projected)')

        # 3. Plot the out-of-tolerance misses
        chr_failed = failed[failed['ref_chr_ref'] == chrom]
        
        # Plot projected coordinate (where we expected it)
        ax.scatter(chr_failed['midpoint_ref'], chr_failed['yhat'], 
                   color='blue', marker='x', s=15, alpha=0.7, label='Failed: Expected (yhat)')
        
        # Plot actual coordinate (where it actually mapped)
        ax.scatter(chr_failed['midpoint_ref'], chr_failed['midpoint_qry'], 
                   color='red', marker='+', s=15, alpha=0.7, label='Failed: Actual')
        
        # Draw a faint line connecting the expectation to the reality
        for _, row in chr_failed.iterrows():
            ax.plot([row['midpoint_ref'], row['midpoint_ref']], 
                    [row['yhat'], row['midpoint_qry']], 
                    color='k', linestyle=':', alpha=0.3, linewidth=0.5)

        # Formatting - Keeping absolute terms
        ax.set_title(f"Synteny Trace & Tolerance Diagnostics: Chromosome {chrom}", pad=15)
        ax.set_xlabel("Reference Midpoint (Absolute bp)", fontsize=12)
        ax.set_ylabel("Query Midpoint (Absolute bp)", fontsize=12)
        
        # Format axes to handle millions gracefully without scientific notation
        ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda y, p: format(int(y), ',')))
        
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, linestyle='--', alpha=0.4)
        
        plt.tight_layout()
        plt.savefig(fig_dir / f"trace_{chrom}.png", dpi=300)
        plt.close()

def main():
    args = parse_args()
    plot_chromosome_traces(args.data_dir, args.figures)
    print(f"[INFO] Diagnostic plots saved to {args.figures}")

if __name__ == "__main__":
    main()