"""Build site-level CRISPR targeting tables from the existing mapper results.

This module deliberately DOES NOT change Bowtie alignments, CFD scores, or the
legacy one-best-pair summary. It projects the relationship evidence back onto
individual alignment rows (the unit needed to find unexpected query targets).
"""

import numpy as np
import pandas as pd


RELATIONSHIP_TYPES = (
    "syntenic_within_tolerance",
    "same_block_outside_tolerance",
    "same_chr_different_block",
    "different_chr",
    "unresolved",
)


def mark_design_source(ref_df, manifest):
    """Tag the ISO1 alignment at the exact original guide discovery locus.

    Source BED coordinates are *local to extracted euchromatin*, exactly like
    the reference CFD coordinates. The absolute source coordinates are kept
    separately in the manifest for manuscript reporting.

    We also compare oriented 23-mer sequences, not merely coordinates, to
    avoid assigning a hit on the opposite strand to the original source.
    """
    source = manifest.set_index("protospacerID", verify_integrity=True)
    ref_df["is_design_source"] = (
        ref_df["ref_chr"].eq(ref_df["protospacerID"].map(source["source_chr"]))
        & ref_df["ref_start"].eq(ref_df["protospacerID"].map(source["source_start_local"]))
        & ref_df["ref_end"].eq(ref_df["protospacerID"].map(source["source_end_local"]))
        & ref_df["genome_target"].eq(ref_df["protospacer"])
    )
    return source


def write_analysis_tables(ref_df, qry_df, merged, source, out_dir, ref_name, qry_name):
    """Preserve every alignment and derive per-guide summaries.

    A qualifying relationship is a same-chromosome, same-block pair with
    projection distance within the selected tolerance. A hit is 'matched'
    if it participates in >=1 such pair; there is NO forced one-to-one match.
    Thus 10x10 repetitive cases may retain all 100 qualifying pairs, while
    each of the 20 genomic alignment rows is counted once in site statistics.
    """
    required = {"protospacerID", "name", "ref_chr", "ref_start", "ref_end",
                "CFD", "pam", "genome_target", "syntenic_block"}
    for side, df in (("reference", ref_df), ("query", qry_df)):
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{side} CFD table missing {sorted(missing)}")
        if not df["name"].is_unique:
            raise ValueError(f"{side} alignment names are not unique")

    qualifying = merged["relationship"].eq("syntenic_within_tolerance")
    matching = merged.loc[qualifying, ["name_ref", "name_qry", "is_design_source_ref"]]
    # Names are unique ONLY within their assembly, so never join ref and query
    # by the '_1', '_2', ... suffix. These are independent identifier spaces.
    matched_ref = pd.Index(matching["name_ref"].unique())
    matched_qry = pd.Index(matching["name_qry"].unique())
    source_qry = pd.Index(matching.loc[matching["is_design_source_ref"], "name_qry"].unique())

    ref_df["has_syntenic_partner"] = ref_df["name"].isin(matched_ref)
    qry_df["has_syntenic_partner"] = qry_df["name"].isin(matched_qry)
    qry_df["matches_design_source"] = qry_df["name"].isin(source_qry)

    # Site classes are independent of the number of Cartesian relationships.
    # A 'query_unmatched' site may be in a changed block, far from a projection,
    # or simply absent from the eligible pair table; it is NOT proof of an
    # evolutionary insertion or proof of no sequence anywhere in ISO1.
    ref_df["site_class"] = np.select(
        [ref_df["is_design_source"], ref_df["syntenic_block"].isna(),
         ref_df["has_syntenic_partner"]],
        ["design_source", "unresolved", "conserved_alternative"],
        default="reference_unmatched",
    )
    qry_df["site_class"] = np.select(
        [qry_df["syntenic_block"].isna(), qry_df["matches_design_source"],
         qry_df["has_syntenic_partner"]],
        ["unresolved", "corresponding_source", "conserved_alternative"],
        default="query_unmatched",
    )

    # Independent per-alignment CFD values (rather than CFD of each pair).
    # The source-site reference score is normally 1.0, but derive it from data.
    design_ref_cfd = (ref_df.loc[ref_df["is_design_source"]]
                      .groupby("protospacerID")["CFD"].max())
    design_qry_cfd = (qry_df.loc[qry_df["matches_design_source"]]
                      .groupby("protospacerID")["CFD"].max())

    # Include *unassigned* alignments too. Assigning them 'unresolved' is safer
    # than incorrectly interpreting failed synteny resolution as a new target.
    provenance_cols = ["source_chr", "source_start_local", "source_end_local",
                       "source_start", "source_end", "source_strand", "source_pam"]
    output_cols = ["protospacerID", "name", "assembly", "ref_chr", "ref_start",
                   "ref_end", "protospacer", "genome_target", "wt", "sg", "pam", "CFD",
                   "mismatches", "midpoint", "syntenic_block",
                   "site_class", "has_syntenic_partner", "matches_design_source",
                   "design_source_CFD_ref", "corresponding_source_CFD_query"] + provenance_cols
    target_path = out_dir / "target_sites.tsv.gz"
    import gzip
    with gzip.open(target_path, "wt", newline="") as fh:
        first = True
        for label, df in ((ref_name, ref_df), (qry_name, qry_df)):
            for start in range(0, len(df), 200_000):
                chunk = df.iloc[start:start + 200_000].copy()
                chunk["assembly"] = label
                if "matches_design_source" not in chunk:
                    chunk["matches_design_source"] = False
                chunk["design_source_CFD_ref"] = chunk["protospacerID"].map(design_ref_cfd)
                chunk["corresponding_source_CFD_query"] = chunk["protospacerID"].map(design_qry_cfd)
                for col in provenance_cols:
                    chunk[col] = chunk["protospacerID"].map(source[col])
                # Source and target PAMs are distinct; 'pam' is the dinucleotide
                # fetched from the *aligned genome* and used in CFD scoring.
                chunk.reindex(columns=output_cols).to_csv(
                    fh, sep="\t", index=False, header=first, na_rep="NA"
                )
                first = False

    # Keep a compact *edge table* in addition to the site table. All Cartesian
    # pairs are retained; its rows are possible relationships, NOT independent
    # genomic sites. Readers can use this to audit ambiguous/repetitive hits.
    edge_cols = ["protospacerID", "name_ref", "name_qry", "ref_chr_ref",
                 "ref_start_ref", "ref_end_ref", "ref_chr_qry", "ref_start_qry",
                 "ref_end_qry", "syntenic_block_ref", "syntenic_block_qry",
                 "same_chr", "same_block", "comparable_proj_dist", "relationship",
                 "is_design_source_ref", "CFD_ref", "CFD_qry"]
    with gzip.open(out_dir / "site_relationships.tsv.gz", "wt", newline="") as fh:
        first = True
        for start in range(0, len(merged), 200_000):
            # Avoid copying the complete 23M-row join a second time.
            chunk = merged.iloc[start:start + 200_000][edge_cols].copy()
            chunk["delta_CFD"] = chunk["CFD_qry"] - chunk["CFD_ref"]
            chunk.to_csv(fh, sep="\t", index=False, header=first, na_rep="NA")
            first = False

    # Aggregate *distinct alignments*; counting nonsyntenic Cartesian edges
    # would inflate the number of novel potential query sites (e.g. PAM_1).
    def counts(df, label):
        # Groupby-sum is vectorized. Avoid a Python lambda over ~8M guide IDs.
        counts_df = df.groupby("protospacerID", sort=False).agg(
            **{f"n_{label}": ("name", "size"),
               f"n_{label}_matched": ("has_syntenic_partner", "sum")}
        )
        counts_df[f"n_{label}_unresolved"] = (
            df.loc[df["syntenic_block"].isna()].groupby("protospacerID").size()
        ).reindex(counts_df.index, fill_value=0)
        resolved_unmatched = (df["syntenic_block"].notna()
                              & ~df["has_syntenic_partner"])
        counts_df[f"n_{label}_unmatched"] = (
            df.loc[resolved_unmatched].groupby("protospacerID").size()
        ).reindex(counts_df.index, fill_value=0)
        return counts_df

    per_guide = counts(ref_df, "ref").join(counts(qry_df, "qry"), how="outer").fillna(0)
    count_cols = [c for c in per_guide if c.startswith("n_")]
    per_guide[count_cols] = per_guide[count_cols].astype("int32")
    per_guide["delta_hits"] = per_guide["n_qry"] - per_guide["n_ref"]
    per_guide["design_source_CFD_ref"] = design_ref_cfd
    per_guide["corresponding_source_CFD_query"] = design_qry_cfd
    query_unmatched = (qry_df["syntenic_block"].notna()
                       & ~qry_df["has_syntenic_partner"])
    per_guide["max_unmatched_qry_CFD"] = (
        qry_df.loc[query_unmatched].groupby("protospacerID")["CFD"].max()
    )
    per_guide.index.name = "protospacerID"
    per_guide.to_csv(out_dir / "guide_targeting_summary.tsv.gz", sep="\t",
                     compression="gzip", na_rep="NA")

    if not ((per_guide["n_ref_matched"] + per_guide["n_ref_unmatched"]
             + per_guide["n_ref_unresolved"] == per_guide["n_ref"]).all()
            and (per_guide["n_qry_matched"] + per_guide["n_qry_unmatched"]
                 + per_guide["n_qry_unresolved"] == per_guide["n_qry"]).all()):
        raise AssertionError("Alignment-level counts do not reconcile")
    print(f"[INFO] Site-level reports: {target_path} and site_relationships.tsv.gz")
    print(f"[INFO] Per-guide report: guide_targeting_summary.tsv.gz ({len(per_guide):,} guides)")
