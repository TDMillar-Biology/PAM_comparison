#!/usr/bin/env python3
"""Map ISO1-derived CRISPR candidates and summarize conservation by tolerance.

Synteny projections and best alignment pairs are computed once here. Downstream
figures must use conservation_by_tolerance.tsv rather than rematching targets.

Definitions (count unique protospacerID, not Bowtie alignment rows):
  shared:           detected in both raw CFD alignment tables
  within_tolerance: shared with a selected projected match at distance <= tol
  exceeds_tolerance: shared with a selected projected match at distance > tol
  synteny_unresolved: shared, but no usable mapped reference/query pair
  ref_only/qry_only: detected in only one raw CFD alignment table

Block agreement is diagnostic, NOT a hard matching criterion: block boundaries
may split otherwise concordant genomic positions.
"""

from svmu2.orchestration.parse import parse
from svmu2.orchestration.synteny import resolve_synteny
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

DEFAULT_TOLERANCES = (100, 300, 1000, 10000, 100000, 1000000)
REQUIRED_CFD_COLUMNS = {"protospacerID", "ref_chr", "ref_start", "ref_end"}


def parse_args():
    parser = argparse.ArgumentParser(description="Synteny-guided PAM orthology and conservation statistics")
    parser.add_argument("--ref-cfd", required=True, type=Path)
    parser.add_argument("--query-cfd", required=True, type=Path)
    parser.add_argument("--delta", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--tol", type=int, default=1000,
                        help="Primary distance threshold for conserved_pams.csv (bp)")
    parser.add_argument("--tolerances", nargs="+", type=int, default=DEFAULT_TOLERANCES,
                        help="Distance thresholds to include in conservation_by_tolerance.tsv (bp)")
    args = parser.parse_args()
    if args.tol < 0 or any(t < 0 for t in args.tolerances):
        parser.error("Distance tolerances must be non-negative")
    args.tolerances = sorted(set(args.tolerances) | {args.tol})
    args.out_dir.mkdir(parents=True, exist_ok=True)
    return args


def load_cfd_table(path: Path) -> pd.DataFrame:
    # Read only once; keep the original columns in matched CSVs for downstream CFD plots.
    df = pd.read_csv(path, dtype={"protospacerID": str, "ref_chr": str,
                                  "ref_start": int, "ref_end": int})
    missing = REQUIRED_CFD_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if df["protospacerID"].isna().any():
        raise ValueError(f"{path}: missing protospacerID values")
    df["midpoint"] = (df["ref_start"] + df["ref_end"]) / 2.0
    return df


def assign_synteny(df, trees, mode="reference"):
    # Preserve the existing rule: only unambiguous single-interval assignments
    # are projected. We do not reject pairs solely for differing block IDs.
    df_out = df.copy()
    block_ids, yhats = [], []
    for row in df.itertuples(index=False):
        tree = trees.get(row.ref_chr)
        x = row.midpoint
        if tree is None:
            block_ids.append(None)
            if mode == "reference":
                yhats.append(None)
            continue
        hits = tree.at(x)
        if len(hits) == 1:
            hit_data = next(iter(hits)).data
            block_ids.append(hit_data.index)
            if mode == "reference":
                yhats.append((hit_data.slope * x) + hit_data.y_intercept)
        else:
            # No hit or multiple overlapping syntenic intervals: unresolved.
            block_ids.append(None)
            if mode == "reference":
                yhats.append(None)
    df_out["syntenic_block"] = block_ids
    if mode == "reference":
        df_out["yhat"] = yhats
    return df_out


def candidate_counts(ref_df: pd.DataFrame, qry_df: pd.DataFrame) -> dict:
    """Count unique candidate IDs before synteny filtering."""
    ref_ids = pd.Index(ref_df["protospacerID"].unique())
    qry_ids = pd.Index(qry_df["protospacerID"].unique())
    n_shared = len(ref_ids.intersection(qry_ids, sort=False))
    return {
        "ref_detected": len(ref_ids),
        "qry_detected": len(qry_ids),
        "shared": n_shared,
        "ref_only": len(ref_ids) - n_shared,
        "qry_only": len(qry_ids) - n_shared,
        "total_candidates": len(ref_ids) + len(qry_ids) - n_shared,
    }


def summarize_tolerances(merged_best: pd.DataFrame, counts: dict, tolerances) -> pd.DataFrame:
    """Exhaustive, mutually exclusive conservation counts for every tolerance."""
    matched = len(merged_best)  # One row per protospacerID after best-pair selection.
    if matched > counts["shared"]:
        raise ValueError("More selected pairs than candidate IDs detected in both tables")
    unresolved = counts["shared"] - matched
    distances = merged_best["proj_dist"].to_numpy(dtype=float)
    if not np.isfinite(distances).all():
        raise ValueError("Selected orthology pairs have nonfinite projection distances")

    ref_blocks = merged_best["syntenic_block_ref"].to_numpy()
    qry_blocks = merged_best["syntenic_block_qry"].to_numpy()
    same_block = ref_blocks == qry_blocks
    # Cross-chromosome matches are diagnostic; a block-ID match alone is not proof of orthology.
    ref_chrom = merged_best["ref_chr_ref"].to_numpy()
    qry_chrom = merged_best["ref_chr_qry"].to_numpy()
    different_chrom = ref_chrom != qry_chrom
    rows = []
    for tol in tolerances:
        within = distances <= tol
        n_within = int(within.sum())
        n_exceeds = matched - n_within
        row = {
            "tolerance_bp": int(tol),
            "ref_detected": counts["ref_detected"],
            "qry_detected": counts["qry_detected"],
            "total_candidates": counts["total_candidates"],
            "shared": counts["shared"],
            "ref_only": counts["ref_only"],
            "qry_only": counts["qry_only"],
            "within_tolerance": n_within,
            "exceeds_tolerance": n_exceeds,
            "synteny_unresolved": unresolved,
            "same_block_within": int((same_block & within).sum()),
            "different_block_within": int((~same_block & within).sum()),
            "same_block_exceeds": int((same_block & ~within).sum()),
            "different_block_exceeds": int((~same_block & ~within).sum()),
            "different_chromosome_within": int((different_chrom & within).sum()),
            "different_chromosome_exceeds": int((different_chrom & ~within).sum()),
        }
        if n_within + n_exceeds + unresolved != counts["shared"]:
            raise AssertionError("Shared-candidate accounting failed")
        if (row["ref_only"] + row["qry_only"] + n_within + n_exceeds + unresolved
                != counts["total_candidates"]):
            raise AssertionError("Mutually exclusive category accounting failed")
        rows.append(row)
    return pd.DataFrame(rows)


def merge_and_filter(ref_df, qry_df, out_dir, tol, tolerances, counts):
    print("[INFO] Selecting nearest projected reference/query alignment per candidate")
    ref_mapped = ref_df.dropna(subset=["syntenic_block", "yhat"])
    qry_mapped = qry_df.dropna(subset=["syntenic_block"])

    # Preserve the original candidate-ID join and minimum-distance selection.
    # Different block IDs are allowed (and counted), because neighboring blocks
    # may differ at genuine syntenic boundaries.


    # Distribution of Cartesian relationships per protospacer
    ref_counts = ref_mapped.groupby("protospacerID").size().rename("n_ref")
    qry_counts = qry_mapped.groupby("protospacerID").size().rename("n_qry")

    pair_counts = pd.concat([ref_counts, qry_counts], axis=1, join="inner")
    pair_counts["n_pairs"] = pair_counts["n_ref"] * pair_counts["n_qry"]

    distribution = (
        pair_counts["n_pairs"]
        .value_counts()
        .sort_index()
        .rename_axis("relationships")
        .reset_index(name="protospacers")
    )

    distribution["percent"] = (
        100 * distribution["protospacers"] / len(pair_counts)
    )

    print("\nCartesian relationship distribution:")
    print(distribution.to_string(index=False))
    print(f"\nMaximum relationships per protospacer: {pair_counts['n_pairs'].max()}")
    print(f"Total candidate relationships: {pair_counts['n_pairs'].sum():,}")
    print("\nReference × query alignment multiplicity:")
    print(pd.crosstab(pair_counts["n_ref"], pair_counts["n_qry"]))

    # Optional: retain the distribution as a TSV
    distribution.to_csv(
        "cartesian_relationship_distribution.tsv",
        sep="\t",
        index=False
    )


    merged = ref_mapped.merge(qry_mapped, on="protospacerID", suffixes=("_ref", "_qry"))
    merged["proj_dist"] = (merged["midpoint_qry"] - merged["yhat"]).abs() # this isnt chromosome aware !!!!


    import numpy as np

    merged = ref_mapped.merge(
        qry_mapped,
        on="protospacerID",
        suffixes=("_ref", "_qry")
    )

    # Chromosome and syntenic-block agreement
    merged["same_chr"] = (
        merged["ref_chr_ref"] == merged["ref_chr_qry"]
    )

    blocks_known = (
        merged["syntenic_block_ref"].notna()
        & merged["syntenic_block_qry"].notna()
    )

    merged["same_block"] = (
        merged["same_chr"]
        & blocks_known
        & (
            merged["syntenic_block_ref"]
            == merged["syntenic_block_qry"]
        )
    )

        # Coordinate projection only makes sense on compatible chromosomes
    merged["proj_dist"] = (
        merged["midpoint_qry"] - merged["yhat"]
    ).abs()

    merged.loc[~merged["same_chr"], "proj_dist"] = np.nan

    # Classify each candidate relationship
    merged["relationship"] = "unresolved"

    merged.loc[
        ~merged["same_chr"] | (blocks_known & ~merged["same_block"]),
        "relationship"
    ] = "nonsyntenic"

    merged.loc[
        merged["same_block"] & (merged["proj_dist"] > args.tol),
        "relationship"
    ] = "same_block_divergent"

    merged.loc[
        merged["same_block"] & (merged["proj_dist"] <= args.tol),
        "relationship"
    ] = "syntenic"

    ### This merge isnt chromosome aware, we could merge first on chrom, then protospacer ID? 



    # Number of alignments in each assembly
    ref_counts = ref_mapped.groupby("protospacerID").size().rename("n_ref")
    qry_counts = qry_mapped.groupby("protospacerID").size().rename("n_qry")

    # Number of relationships of each type per protospacer
    relationship_counts = (
        merged.groupby(["protospacerID", "relationship"])
        .size()
        .unstack(fill_value=0)
    )

    pair_summary = pd.concat(
        [ref_counts, qry_counts, relationship_counts],
        axis=1,
        join="inner"
    ).fillna(0)

    pair_summary["n_pairs"] = (
        pair_summary["n_ref"] * pair_summary["n_qry"]
    )

    # Aggregate by reference/query alignment multiplicity
    summary = (
        pair_summary.groupby(["n_ref", "n_qry"])
        .agg(
            protospacers=("n_pairs", "size"),
            total_pairs=("n_pairs", "sum"),
            syntenic_pairs=("syntenic", "sum"),
            nonsyntenic_pairs=("nonsyntenic", "sum"),
            divergent_pairs=("same_block_divergent", "sum"),
            unresolved_pairs=("unresolved", "sum")
        )
        .reset_index()
    )

    print(summary.to_string(index=False))

    pair_summary.to_csv(
        "protospacer_relationship_summary.tsv",
        sep="\t"
    )

    summary.to_csv(
        "relationship_multiplicity_summary.tsv",
        sep="\t",
        index=False
    )


    # Preserve the publication pipeline's best-distance selection exactly.
    # (Tied distances inherit the row ordering of the input CFD tables.)
    merged_best = (merged.sort_values("proj_dist")
                   .drop_duplicates(subset="protospacerID", keep="first"))

    conserved = merged_best[merged_best["proj_dist"] <= tol]
    failed_tol = merged_best[merged_best["proj_dist"] > tol]
    conserved.to_csv(out_dir / "conserved_pams.csv", index=False)
    failed_tol.to_csv(out_dir / "failed_tolerance_pams.csv", index=False)

    summary = summarize_tolerances(merged_best, counts, tolerances)
    summary.to_csv(out_dir / "conservation_by_tolerance.tsv", sep="\t", index=False)
    print("[INFO] Unique-candidate conservation counts (not Bowtie alignment rows):")
    print(summary[["tolerance_bp", "ref_only", "qry_only", "within_tolerance",
                   "exceeds_tolerance", "synteny_unresolved"]].to_string(index=False))
    import pdb
    pdb.set_trace()
    return summary


def main():
    args = parse_args()
    ref_cfd = load_cfd_table(args.ref_cfd)
    qry_cfd = load_cfd_table(args.query_cfd)
    counts = candidate_counts(ref_cfd, qry_cfd)
    print(f"[INFO] Unique candidate IDs: {counts['ref_detected']:,} ref; "
          f"{counts['qry_detected']:,} query; {counts['shared']:,} shared")

    _, primary_alns = parse(args.delta)
    resolve_synteny(primary_alignments=primary_alns, breakpoint_map=None)
    ref_trees = {aln.reference: aln.reference_synteny_tree for aln in primary_alns.values()}
    qry_trees = {aln.query: aln.query_synteny_tree for aln in primary_alns.values()}

    ref_df = assign_synteny(ref_cfd, ref_trees, mode="reference") # tracks a unique name and protospacerID
    qry_df = assign_synteny(qry_cfd, qry_trees, mode="query")

    unmapped_ref = ref_df[ref_df["syntenic_block"].isna()]
    unmapped_qry = qry_df[qry_df["syntenic_block"].isna()]
    unmapped_ref.to_csv(args.out_dir / "unmapped_ref.csv", index=False)
    unmapped_qry.to_csv(args.out_dir / "unmapped_qry.csv", index=False)
    print(f"[INFO] Alignment rows without unique synteny: "
          f"{len(unmapped_ref):,} reference; {len(unmapped_qry):,} query")

    merge_and_filter(ref_df, qry_df, args.out_dir, args.tol, args.tolerances, counts)

if __name__ == "__main__":
    main()
