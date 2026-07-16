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
    #chromosomes = pd.concat([conserved['ref_chr_ref'], unmapped_ref['ref_chr']]).dropna().unique()
    import pdb 
    chromosomes = pd.concat([failed['ref_chr_ref']]).dropna().unique()

    value_counts = conserved['protospacerID'].value_counts()
    pdb.set_trace()
    fig, ax = plt.subplots(figsize= (12,10))
    plt.hist(value_counts)
    plt.show()

    for chrom in chromosomes:
        print(f"[INFO] Rendering trace for {chrom}...")
        fig, ax = plt.subplots(figsize=(12, 10))
        
        sub = conserved[conserved['ref_chr_ref'] == chrom]

        plt.scatter(sub['midpoint_ref'], sub['midpoint_qry'], color = 'Red', label = 'mapped midpoint - truth')
        plt.scatter(sub['midpoint_ref'], sub['yhat'], color = 'Blue', label = 'yhat - interpolated')
        plt.show()

        sub['difference'] = sub['yhat'] - sub['midpoint_qry']
        plt.hist(sub['yhat'] - sub['midpoint_qry'], bins = 100, log=True)
        plt.show()
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
        #plt.legend()
        #plt.show()
        plt.savefig(fig_dir / f"trace_{chrom}.png", dpi=300)
        plt.close()

def main():
    args = parse_args()
    plot_chromosome_traces(args.data_dir, args.figures)
    print(f"[INFO] Diagnostic plots saved to {args.figures}")

if __name__ == "__main__":
    main()