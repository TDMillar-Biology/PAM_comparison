#!/usr/bin/env python3
import csv
import argparse
import matplotlib.pyplot as plt
from collections import Counter
from pathlib import Path

def plot_mapping_distribution(csv_path, out_path, tail_cutoff=15):
    print(f"[INFO] Streaming {csv_path} to count hits per protospacer...")
    
    # Step 1: Stream the CSV and count hits per protospacerID
    # We use Counter to keep memory footprint tiny, even for 16M rows.
    hits_per_guide = Counter()
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            hits_per_guide[row['protospacerID']] += 1

    print(f"[INFO] Processed {len(hits_per_guide)} unique protospacers.")
    print("[INFO] Aggregating distribution...")

    # Step 2: Count the frequency of those mapping counts
    # (e.g., how many guides mapped exactly 1 time, 2 times, etc.)
    distribution = Counter(hits_per_guide.values())

    # Step 3: Bin the long tail for a cleaner plot
    binned_dist = Counter()
    for hits, count in distribution.items():
        if hits >= tail_cutoff:
            binned_dist[f"{tail_cutoff}+"] += count
        else:
            binned_dist[str(hits)] += count

    # Prepare data for plotting in sorted order (1 to tail_cutoff+)
    labels = [str(i) for i in range(1, tail_cutoff)] + [f"{tail_cutoff}+"]
    counts = [binned_dist[label] for label in labels]

    # Step 4: Generate the visualization
    print("[INFO] Generating plot...")
    plt.figure(figsize=(10, 6))
    
    bars = plt.bar(labels, counts, color='#2c7fb8', edgecolor='black', zorder=3)
    
    #plt.yscale('log')
    plt.xlabel('Number of Mapped Genomic Loci (Up to 3 Mismatches)', fontsize=12)
    plt.ylabel('Number of Unique Protospacers (Log Scale)', fontsize=12)
    plt.title('Distribution of Mapped Locations per Discovered PAM', fontsize=14)
    
    # Add a subtle grid behind the bars for easier log-scale reading
    plt.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)

    # Add text labels on top of each bar so the exact numbers are readable
    for bar, count in zip(bars, counts):
        if count > 0:
            plt.text(
                bar.get_x() + bar.get_width() / 2, 
                count * 1.2,  # Offset slightly above the bar
                f"{count:,}", 
                ha='center', va='bottom', fontsize=9, rotation=45
            )

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    
    print(f"[SUCCESS] Plot saved to {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Plot the distribution of mapped loci per protospacer.")
    parser.add_argument("-i", "--input", required=True, type=Path, help="Input raw CFD CSV file")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output PNG file path")
    parser.add_argument("-c", "--cutoff", type=int, default=30, help="Bin counts >= this number into a single tail bin")
    
    args = parser.parse_args()

    if not args.input.exists():
        print(f"[FATAL] Input file not found: {args.input}")
        return

    # Ensure output directory exists (using the Pathlib division operator!)
    if not args.output.parent.exists():
        args.output.parent.mkdir(parents=True, exist_ok=True)

    plot_mapping_distribution(args.input, args.output, args.cutoff)

if __name__ == "__main__":
    main()