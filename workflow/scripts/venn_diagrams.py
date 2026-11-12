#!/usr/bin/env python3
"""Plot mapper-computed, mutually exclusive PAM conservation categories.

This module deliberately contains NO synteny-resolution, coordinate projection,
or alignment-pairing code. Conservation counts are produced exclusively by
synteny_mapper.py in conservation_by_tolerance.tsv.
"""

import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REQUIRED_COLUMNS = {
    "tolerance_bp", "total_candidates", "shared", "ref_only", "qry_only",
    "within_tolerance", "exceeds_tolerance", "synteny_unresolved",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Plot PAM conservation counts from the synteny mapper")
    parser.add_argument("--summary", required=True, type=Path,
                        help="conservation_by_tolerance.tsv produced by synteny_mapper.py")
    parser.add_argument("--figures", required=True, type=Path,
                        help="Directory for multi_tolerance_sunburst.png")
    args = parser.parse_args()
    args.figures.mkdir(parents=True, exist_ok=True)
    return args


def load_summary(path):
    df = pd.read_csv(path, sep="\t")
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns {sorted(missing)}")
    if df["tolerance_bp"].duplicated().any():
        raise ValueError("Tolerance summary contains duplicate tolerance thresholds")
    for row in df.itertuples(index=False):
        if row.shared != row.within_tolerance + row.exceeds_tolerance + row.synteny_unresolved:
            raise ValueError(f"Shared-candidate categories do not sum at {row.tolerance_bp} bp")
        if row.total_candidates != row.ref_only + row.qry_only + row.shared:
            raise ValueError(f"Candidate categories do not sum at {row.tolerance_bp} bp")
    return df.sort_values("tolerance_bp")


def plot_tolerance_venns(summary, output_dir):
    n = len(summary)
    cols = 3
    nrows = (n + cols - 1) // cols
    fig, axes = plt.subplots(nrows, cols, figsize=(20, 10 * nrows))
    axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

    # Keep the double-ring presentation, but add an explicit unresolved category.
    palette = plt.get_cmap("Pastel1")
    ref_color, qry_color, shared_color = palette(0), palette(1), palette(2)
    blues = plt.get_cmap("Blues")
    for ax, row in zip(axes, summary.itertuples(index=False)):
        total = int(row.total_candidates)
        def pct(value):
            return f"{value / total * 100:.1f}%" if total else "0.0%"

        ref_only = int(row.ref_only)
        qry_only = int(row.qry_only)
        shared = int(row.shared)
        within = int(row.within_tolerance)
        exceeds = int(row.exceeds_tolerance)
        unresolved = int(row.synteny_unresolved)

        inner_sizes = [ref_only, qry_only, shared]
        outer_sizes = [ref_only, qry_only, within, exceeds, unresolved]
        inner_labels = ["", "", f"Shared\n{shared:,}\n({pct(shared)})"]
        outer_labels = [
            f"Ref only\n{ref_only:,}\n({pct(ref_only)})",
            f"Query only\n{qry_only:,}\n({pct(qry_only)})",
            f"Within tol.\n{within:,}\n({pct(within)})",
            f"Exceeds tol.\n{exceeds:,}\n({pct(exceeds)})",
            f"Unresolved\n{unresolved:,}\n({pct(unresolved)})",
        ]
        if total:
            ax.pie(outer_sizes, radius=1, labels=outer_labels,
                   labeldistance=1.15,
                   colors=[ref_color, qry_color, blues(0.4), blues(0.2), palette(7)],
                   wedgeprops=dict(width=0.35, edgecolor="white", linewidth=2),
                   textprops={"fontsize": 9})
            ax.pie(inner_sizes, radius=0.65, labels=inner_labels,
                   labeldistance=0.75,
                   colors=[ref_color, qry_color, shared_color],
                   wedgeprops=dict(width=0.35, edgecolor="white", linewidth=2),
                   textprops={"fontsize": 10, "weight": "bold"})
        else:
            ax.text(0.5, 0.5, "No detected candidates", ha="center", va="center",
                    transform=ax.transAxes)
        ax.set_title(f"Tolerance: ±{int(row.tolerance_bp):,} bp", fontsize=13, pad=30)

    for ax in axes[n:]:
        ax.set_visible(False)
    fig.suptitle("Synteny-Guided PAM Conservation Across Spatial Tolerances", fontsize=18, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    output = output_dir / "multi_tolerance_sunburst.png"
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"[INFO] Figure saved to {output}")


def main():
    args = parse_args()
    summary = load_summary(args.summary)
    plot_tolerance_venns(summary, args.figures)


if __name__ == "__main__":
    main()
