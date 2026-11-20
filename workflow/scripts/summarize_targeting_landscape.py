#!/usr/bin/env python3
"""Summarize guide-level target multiplicity and site-level correspondence.

Input: guide_targeting_summary.tsv[.gz] from target_site_reporting.py.
Each row represents one reference-derived guide (protospacerID). Its n_ref
and n_qry fields count *genomic alignments*, not Cartesian ref/query pairs.

The script streams the input in chunks so a full-cohort summary does not need
all guide records in memory. It produces a readable report and tidy TSV for
reproducible downstream plotting. Neither assembly name is hard-coded.
"""

import argparse
from pathlib import Path

import pandas as pd


# Keep the upstream column names intact: these are part of the output contract
# of workflow/scripts/target_site_reporting.py.
REQUIRED_COLUMNS = (
    "n_ref", "n_qry", "delta_hits",
    "n_ref_matched", "n_qry_matched",
    "n_ref_unmatched", "n_qry_unmatched",
    "n_ref_unresolved", "n_qry_unresolved",
)

GUIDE_METRICS = (
    "total_guides", "equal_multiplicity", "query_enriched", "query_depleted",
    "guides_with_query_unmatched", "guides_with_ref_unmatched",
    "guides_with_unresolved",
)
SITE_METRICS = (
    "ref_sites", "qry_sites", "ref_matched", "qry_matched",
    "ref_unmatched", "qry_unmatched", "ref_unresolved", "qry_unresolved",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--guide-summary", type=Path, required=True,
                        help="Per-guide TSV from map_pam_orthology (.tsv or .tsv.gz)")
    parser.add_argument("--report", type=Path, required=True,
                        help="Destination for the human-readable text report")
    parser.add_argument("--table", type=Path, required=True,
                        help="Destination for the tidy, machine-readable TSV")
    parser.add_argument("--ref-label", default="Reference",
                        help="Reference assembly label (display only)")
    parser.add_argument("--query-label", default="Query",
                        help="Query assembly label (display only)")
    parser.add_argument("--chunksize", type=int, default=200_000,
                        help="Number of guide rows processed per chunk (default 200000)")
    args = parser.parse_args()
    if args.chunksize < 1:
        parser.error("--chunksize must be a positive integer")
    return args


def count_summary(path: Path, chunksize: int) -> dict:
    """Accumulate integer counts, keeping a constant-memory working set."""
    stats = {name: 0 for name in (*GUIDE_METRICS, *SITE_METRICS)}

    # Pandas automatically detects gzip by filename; an uncompressed .tsv is
    # read with exactly the same code. Only the needed columns are loaded.
    for df in pd.read_csv(path, sep="\t", usecols=list(REQUIRED_COLUMNS),
                          chunksize=chunksize):
        # A missing count would make it impossible to reconcile site totals.
        if df[list(REQUIRED_COLUMNS)].isna().any().any():
            raise ValueError("Input contains missing values in required count columns")

        # Site counts must partition into matched, resolved-unmatched, and
        # unresolved. Check this before trusting the manuscript summary.
        for side in ("ref", "qry"):
            accounted = (df[f"n_{side}_matched"] + df[f"n_{side}_unmatched"]
                         + df[f"n_{side}_unresolved"])
            if not accounted.eq(df[f"n_{side}"]).all():
                raise ValueError(f"{side} alignment counts do not reconcile in input")
        if not df["delta_hits"].eq(df["n_qry"] - df["n_ref"]).all():
            raise ValueError("delta_hits does not equal n_qry - n_ref")

        stats["total_guides"] += len(df)
        stats["equal_multiplicity"] += int((df["n_ref"] == df["n_qry"]).sum())
        stats["query_enriched"] += int((df["delta_hits"] > 0).sum())
        stats["query_depleted"] += int((df["delta_hits"] < 0).sum())
        stats["guides_with_query_unmatched"] += int((df["n_qry_unmatched"] > 0).sum())
        stats["guides_with_ref_unmatched"] += int((df["n_ref_unmatched"] > 0).sum())
        stats["guides_with_unresolved"] += int(
            ((df["n_ref_unresolved"] + df["n_qry_unresolved"]) > 0).sum()
        )

        # Each alignment is counted once, regardless of how many candidate
        # ref/query pairs it appears in during the Cartesian synteny merge.
        for metric, column in (
            ("ref_sites", "n_ref"), ("qry_sites", "n_qry"),
            ("ref_matched", "n_ref_matched"), ("qry_matched", "n_qry_matched"),
            ("ref_unmatched", "n_ref_unmatched"),
            ("qry_unmatched", "n_qry_unmatched"),
            ("ref_unresolved", "n_ref_unresolved"),
            ("qry_unresolved", "n_qry_unresolved"),
        ):
            stats[metric] += int(df[column].sum())

    if not stats["total_guides"]:
        raise ValueError(f"No guide rows found in {path}")
    if (stats["equal_multiplicity"] + stats["query_enriched"]
            + stats["query_depleted"] != stats["total_guides"]):
        raise AssertionError("Multiplicity categories do not sum to total guides")
    for side in ("ref", "qry"):
        if (stats[f"{side}_matched"] + stats[f"{side}_unmatched"]
                + stats[f"{side}_unresolved"] != stats[f"{side}_sites"]):
            raise AssertionError(f"{side} site categories do not sum to site total")
    return stats


def tidy_table(stats: dict) -> pd.DataFrame:
    """Use guide or site totals as denominators, never Cartesian pair counts."""
    records = []
    for metric in GUIDE_METRICS:
        records.append(("guides", metric, stats[metric], stats["total_guides"]))
    for side in ("ref", "qry"):
        total = stats[f"{side}_sites"]
        for metric in (f"{side}_sites", f"{side}_matched",
                       f"{side}_unmatched", f"{side}_unresolved"):
            records.append((f"{side}_sites", metric, stats[metric], total))
    table = pd.DataFrame(records, columns=["group", "metric", "count", "denominator"])
    table["percent"] = (
        100 * table["count"] / table["denominator"].where(table["denominator"] != 0)
    )
    return table


def readable_report(stats: dict, ref_label: str, query_label: str) -> str:
    lines = [
        "GUIDE-LEVEL AND SITE-LEVEL TARGETING SUMMARY",
        f"Reference: {ref_label}",
        f"Query:     {query_label}",
        "",
        "COUNTS",
    ]
    for key in (*GUIDE_METRICS, *SITE_METRICS):
        lines.append(f"{key:35s} {stats[key]:>12,}")
    lines += ["", "PERCENTAGE OF GUIDES (denominator: total_guides)"]
    for key in GUIDE_METRICS[1:]:
        pct = 100 * stats[key] / stats["total_guides"]
        lines.append(f"{key:35s} {pct:>11.3f}%")
    lines += ["", "PERCENTAGE OF ALIGNMENTS (within each assembly)"]
    for side in ("ref", "qry"):
        total = stats[f"{side}_sites"]
        for key in (f"{side}_matched", f"{side}_unmatched", f"{side}_unresolved"):
            pct = 100 * stats[key] / total if total else float("nan")
            lines.append(f"{key:35s} {pct:>11.3f}%")
    lines += [
        "",
        "NOTES",
        "- Multiplicity compares retained target alignments per guide, not pairings.",
        "- 'Unmatched' means resolved but without a qualifying syntenic partner",
        "  under the upstream mapper's tolerance and block criteria.",
        "- 'Unresolved' sites are not automatically query-specific or reference-specific.",
        "- Query-enriched/depleted describes detection relative to the reference,",
        "  not an evolutionary gain or loss.",
    ]
    return "\n".join(lines) + "\n"


def main():
    args = parse_args()
    stats = count_summary(args.guide_summary, args.chunksize)
    table = tidy_table(stats)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.table.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(readable_report(stats, args.ref_label,
                                          args.query_label))
    table.to_csv(args.table, sep="\t", index=False, float_format="%.6f", na_rep="NA")
    print(f"[INFO] Summarized {stats['total_guides']:,} guides")
    print(f"[INFO] Report: {args.report}")
    print(f"[INFO] Table:  {args.table}")


if __name__ == "__main__":
    main()
