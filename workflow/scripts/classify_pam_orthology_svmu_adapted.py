#!/usr/bin/env python3
"""
classify_pam_orthology_block_first.py

Standalone synteny-aware PAM orthology classifier.

This script uses svmu2 only for delta parsing and synteny resolution. It keeps
CRISPR/PAM-specific logic outside of svmu2.

Core algorithm
--------------
1. Parse a nucmer delta file with svmu2 and resolve primary synteny.
2. Convert svmu2 primary synteny blocks into lightweight ref/query blocks.
3. Build interval indexes in reference space and query space.
4. Assign every reference CFD hit to a reference-side synteny block.
5. Assign every query CFD hit to a query-side synteny block.
6. Compare hits block-by-block and protospacerID-by-protospacerID.
7. For each reference hit, project its reference midpoint into query space and
   choose the query hit closest to that projected coordinate, requiring that it
   be within --tol bp.

This replaces the older annotate-then-merge pattern with one biologically
explicit PAM orthology classifier.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd
from intervaltree import IntervalTree

from svmu2.orchestration.parse import parse
from svmu2.orchestration.synteny import resolve_synteny


# -----------------------------------------------------------------------------
# Data models
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class SyntenyBlock:
    """Lightweight reference/query synteny block wrapper."""

    ref_chr: str
    ref_start: int
    ref_end: int
    query_chr: str
    query_start: int
    query_end: int
    block_id: str
    original: object | None = None

    @property
    def ref_min(self) -> int:
        return min(self.ref_start, self.ref_end)

    @property
    def ref_max(self) -> int:
        return max(self.ref_start, self.ref_end)

    @property
    def query_min(self) -> int:
        return min(self.query_start, self.query_end)

    @property
    def query_max(self) -> int:
        return max(self.query_start, self.query_end)

    @property
    def orientation(self) -> int:
        """Return +1 for co-oriented blocks and -1 for inverted blocks."""
        return 1 if self.query_end >= self.query_start else -1

    def project_ref_to_query(self, ref_pos: int) -> float:
        """Project a reference coordinate into query coordinate space.

        The stored query_start/query_end retain orientation. The interval-tree
        indexes may normalize spans, but projection should preserve inversion
        direction.
        """
        ref_span = self.ref_max - self.ref_min
        if ref_span == 0:
            return float(self.query_start)

        frac = (ref_pos - self.ref_min) / ref_span
        frac = max(0.0, min(1.0, frac))

        if self.orientation == 1:
            return self.query_min + frac * (self.query_max - self.query_min)
        return self.query_max - frac * (self.query_max - self.query_min)


@dataclass(frozen=True)
class Hit:
    """A compact CFD hit record assigned to zero/one/multiple synteny blocks."""

    row_index: int
    name: str
    protospacer_id: str
    chrom: str
    start: int
    end: int
    midpoint: int
    cfd: Optional[float]
    source: str  # "ref" or "query"


@dataclass(frozen=True)
class AssignedHit:
    hit: Hit
    block: SyntenyBlock


@dataclass(frozen=True)
class QueryCandidate:
    hit: Hit
    distance_to_expected: float


# -----------------------------------------------------------------------------
# Robust adapters around svmu2 objects
# -----------------------------------------------------------------------------

def _block_from_svmu2_object(block: object, fallback_index: int | str) -> SyntenyBlock:
    """Convert an svmu2 primary synteny block into a local SyntenyBlock.

    This is the only place where this script depends on the concrete svmu2
    block schema. Keep it strict: if svmu2 changes these attributes, fail
    here rather than fishing for alternate names throughout the PAM code.
    """
    required = [
        "reference",
        "query",
        "reference_start",
        "reference_end",
        "query_start",
        "query_end",
    ]
    missing = [attr for attr in required if not hasattr(block, attr)]
    if missing:
        raise AttributeError(
            f"SVMU2 block is missing required attribute(s): {missing}. "
            "Update _block_from_svmu2_object() for this svmu2 version."
        )

    ref_chr = str(block.reference)
    query_chr = str(block.query)
    ref_start = int(block.reference_start)
    ref_end = int(block.reference_end)
    query_start = int(block.query_start)
    query_end = int(block.query_end)

    # index is useful if svmu2 provides it; otherwise use deterministic
    # enumeration supplied by the caller. This keeps block IDs stable within
    # an alignment without requiring svmu2 to expose an ID field.
    block_index = getattr(block, "index", fallback_index)
    block_id = f"{ref_chr}_{query_chr}_{block_index}"

    return SyntenyBlock(
        ref_chr=ref_chr,
        ref_start=ref_start,
        ref_end=ref_end,
        query_chr=query_chr,
        query_start=query_start,
        query_end=query_end,
        block_id=str(block_id),
        original=block,
    )


def _iter_alignment_blocks(aln: object) -> Iterable[object]:
    """Return the svmu2-resolved primary synteny blocks for one alignment."""
    blocks = getattr(aln, "primary_synteny_blocks", None)
    if blocks is None:
        raise AttributeError(
            f"Alignment {getattr(aln, 'reference', '<unknown>')} -> "
            f"{getattr(aln, 'query', '<unknown>')} has no primary_synteny_blocks. "
            "The svmu2 synteny resolver must populate this attribute."
        )
    return blocks


def _resolve_primary_synteny(delta_file: str | Path) -> dict:
    """Parse delta and ask svmu2 to resolve primary synteny.

    Preferred path: svmu2.orchestration.synteny.resolve_synteny().
    Fallback path: if immutable svmu2 fails while building its interval trees
    due to reversed query intervals, re-parse and call the core synteny resolver
    directly, then build safe local interval indexes in this script.
    """
    _, primary = parse(str(delta_file))

    try:
        return resolve_synteny(primary)
    except ValueError as exc:
        # Some svmu2 versions build query interval trees without normalizing
        # inverted blocks. The synteny blocks themselves are still useful, so
        # avoid modifying svmu2 and build safe local indexes below.
        if "Null Interval" not in str(exc) and "IntervalTree" not in str(exc):
            raise

        print(
            "[WARN] svmu2 resolve_synteny failed while building an interval tree. "
            "Falling back to svmu2.core.synteny.resolve_alignment_synteny() and "
            "building normalized local indexes."
        )

        _, primary = parse(str(delta_file))
        try:
            from svmu2.core.synteny import resolve_alignment_synteny
        except ImportError as import_exc:
            raise RuntimeError(
                "Could not import svmu2.core.synteny.resolve_alignment_synteny "
                "for fallback synteny resolution."
            ) from import_exc

        for aln in primary.values():
            resolve_alignment_synteny(aln, expected_breakpoints=0)

        return primary


def build_synteny_blocks_from_delta(delta_file: str | Path) -> list[SyntenyBlock]:
    """Parse delta with svmu2, resolve synteny, and return lightweight blocks."""
    primary = _resolve_primary_synteny(delta_file)

    blocks: list[SyntenyBlock] = []
    used_ids: set[str] = set()

    for aln_key, aln in primary.items():
        aln_ref = str(getattr(aln, "reference", aln_key))
        aln_query = str(getattr(aln, "query", "query"))
        for i, block in enumerate(_iter_alignment_blocks(aln)):
            fallback_index = f"{aln_ref}_{aln_query}_{i}"
            b = _block_from_svmu2_object(block, fallback_index=fallback_index)

            if b.block_id in used_ids:
                raise ValueError(
                    f"Duplicate synteny block ID generated: {b.block_id}. "
                    "Block IDs must be unique for block-local PAM matching."
                )

            used_ids.add(b.block_id)
            blocks.append(b)

    if not blocks:
        raise ValueError(f"No synteny blocks found in delta file: {delta_file}")

    return blocks


# -----------------------------------------------------------------------------
# Interval indexes
# -----------------------------------------------------------------------------

def _safe_interval(start: int, end: int) -> tuple[int, int]:
    """Return a non-empty half-open interval for intervaltree.

    SVMU/genomic coordinates are treated as inclusive endpoints here, while
    intervaltree stores half-open intervals [begin, end). Therefore the upper
    coordinate is incremented by one after normalization.
    """
    lo = int(min(start, end))
    hi = int(max(start, end)) + 1
    if hi <= lo:
        hi = lo + 1
    return lo, hi


def build_block_indexes(
    blocks: list[SyntenyBlock],
) -> tuple[dict[str, IntervalTree], dict[str, IntervalTree], dict[str, SyntenyBlock]]:
    """Build reference-space and query-space interval trees keyed by chromosome."""
    ref_trees: dict[str, IntervalTree] = defaultdict(IntervalTree)
    query_trees: dict[str, IntervalTree] = defaultdict(IntervalTree)
    block_by_id: dict[str, SyntenyBlock] = {}

    for block in blocks:
        r0, r1 = _safe_interval(block.ref_min, block.ref_max)
        q0, q1 = _safe_interval(block.query_min, block.query_max)
        ref_trees[block.ref_chr][r0:r1] = block
        query_trees[block.query_chr][q0:q1] = block
        block_by_id[block.block_id] = block

    return dict(ref_trees), dict(query_trees), block_by_id


# -----------------------------------------------------------------------------
# CFD table handling and block assignment
# -----------------------------------------------------------------------------

REQUIRED_CFD_COLUMNS = {
    "protospacerID",
    "ref_chr",
    "ref_start",
    "ref_end",
    "ref_midpoint",
}


def load_cfd_table(path: str | Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = REQUIRED_CFD_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"{label} CFD table is missing columns: {sorted(missing)}")

    df = df.copy()
    df["ref_chr"] = df["ref_chr"].astype(str)
    df["protospacerID"] = df["protospacerID"].astype(str)
    for col in ["ref_start", "ref_end", "ref_midpoint"]:
        df[col] = df[col].astype(int)

    if "name" not in df.columns:
        df["name"] = df["protospacerID"].astype(str) + "_row" + df.index.astype(str)

    if "CFD" not in df.columns:
        df["CFD"] = pd.NA

    return df


def row_to_hit(row_index: int, row: pd.Series, source: str) -> Hit:
    return Hit(
        row_index=int(row_index),
        name=str(row.get("name", row_index)),
        protospacer_id=str(row["protospacerID"]),
        chrom=str(row["ref_chr"]),
        start=int(row["ref_start"]),
        end=int(row["ref_end"]),
        midpoint=int(row["ref_midpoint"]),
        cfd=None if pd.isna(row.get("CFD", pd.NA)) else float(row["CFD"]),
        source=source,
    )


def assign_hits_to_blocks(
    df: pd.DataFrame,
    trees: dict[str, IntervalTree],
    source: str,
) -> tuple[dict[str, list[AssignedHit]], list[dict[str, object]]]:
    """
    Assign CFD hits to synteny blocks using interval overlap.

    Returns:
        assigned_by_block: block_id -> list[AssignedHit]
        unassigned_rows: diagnostic rows for no-block or multi-block hits
    """
    assigned_by_block: dict[str, list[AssignedHit]] = defaultdict(list)
    diagnostics: list[dict[str, object]] = []

    for row_index, row in df.iterrows():
        hit = row_to_hit(int(row_index), row, source=source)
        tree = trees.get(hit.chrom)

        if tree is None:
            diagnostics.append(hit_diagnostic(hit, status=f"NO_{source.upper()}_BLOCK"))
            continue

        start, end = _safe_interval(hit.start, hit.end)
        overlaps = sorted(tree.overlap(start, end), key=lambda iv: (iv.begin, iv.end))

        if not overlaps:
            # Fallback to midpoint in case coordinate conventions differ subtly.
            overlaps = sorted(tree.at(hit.midpoint), key=lambda iv: (iv.begin, iv.end))

        if not overlaps:
            diagnostics.append(hit_diagnostic(hit, status=f"NO_{source.upper()}_BLOCK"))
            continue

        if len(overlaps) > 1:
            block_ids = ";".join(str(iv.data.block_id) for iv in overlaps)
            diagnostics.append(
                hit_diagnostic(
                    hit,
                    status=f"AMBIGUOUS_{source.upper()}_BLOCK",
                    block_id=block_ids,
                    n_blocks=len(overlaps),
                )
            )
            continue

        block = overlaps[0].data
        assigned_by_block[block.block_id].append(AssignedHit(hit=hit, block=block))

    return dict(assigned_by_block), diagnostics


def hit_diagnostic(
    hit: Hit,
    status: str,
    block_id: object = pd.NA,
    n_blocks: int = 0,
) -> dict[str, object]:
    prefix = hit.source
    return {
        "status": status,
        "protospacerID": hit.protospacer_id,
        f"{prefix}_name": hit.name,
        f"{prefix}_chr": hit.chrom,
        f"{prefix}_start": hit.start,
        f"{prefix}_end": hit.end,
        f"{prefix}_midpoint": hit.midpoint,
        f"{prefix}_CFD": hit.cfd,
        "block_id": block_id,
        "n_blocks": n_blocks,
    }


# -----------------------------------------------------------------------------
# Block-local comparison
# -----------------------------------------------------------------------------

def group_query_hits_by_pid(
    assigned_query_hits: list[AssignedHit],
) -> dict[str, list[Hit]]:
    grouped: dict[str, list[Hit]] = defaultdict(list)
    for assigned in assigned_query_hits:
        grouped[assigned.hit.protospacer_id].append(assigned.hit)
    return dict(grouped)


def compare_block_hits(
    block: SyntenyBlock,
    ref_hits: list[AssignedHit],
    query_hits: list[AssignedHit],
    tol: int,
) -> list[dict[str, object]]:
    """Compare reference/query hits assigned to the same synteny block."""
    query_by_pid = group_query_hits_by_pid(query_hits)
    rows: list[dict[str, object]] = []

    for assigned_ref in ref_hits:
        ref_hit = assigned_ref.hit
        expected = block.project_ref_to_query(ref_hit.midpoint)
        same_pid_hits = query_by_pid.get(ref_hit.protospacer_id, [])

        candidates: list[QueryCandidate] = []
        for query_hit in same_pid_hits:
            distance = abs(query_hit.midpoint - expected)
            if distance <= tol:
                candidates.append(QueryCandidate(hit=query_hit, distance_to_expected=distance))

        candidates.sort(
            key=lambda c: (
                c.distance_to_expected,
                float("inf") if c.hit.cfd is None else -c.hit.cfd,
            )
        )

        base = reference_output_base(ref_hit, block, expected, len(same_pid_hits))

        if not same_pid_hits:
            rows.append(empty_query_output(base, status="LOST_IN_QUERY_BLOCK", n_within_tolerance=0))
            continue

        if not candidates:
            rows.append(
                empty_query_output(
                    base,
                    status="NON_POSITIONAL_MATCH",
                    n_within_tolerance=0,
                )
            )
            continue

        best = candidates[0]
        status = "SYNTENIC_MULTI" if len(candidates) > 1 else "SYNTENIC_SINGLE"
        rows.append(
            {
                **base,
                "status": status,
                "query_name": best.hit.name,
                "query_chr": best.hit.chrom,
                "query_start": best.hit.start,
                "query_end": best.hit.end,
                "query_midpoint": best.hit.midpoint,
                "query_CFD": best.hit.cfd,
                "distance_to_expected": best.distance_to_expected,
                "n_within_tolerance": len(candidates),
                "all_candidate_names": ";".join(c.hit.name for c in candidates),
                "all_candidate_midpoints": ";".join(str(c.hit.midpoint) for c in candidates),
            }
        )

    return rows


def reference_output_base(
    ref_hit: Hit,
    block: SyntenyBlock,
    expected_query_pos: float,
    n_same_block_pid_candidates: int,
) -> dict[str, object]:
    return {
        "protospacerID": ref_hit.protospacer_id,
        "ref_name": ref_hit.name,
        "ref_chr": ref_hit.chrom,
        "ref_start": ref_hit.start,
        "ref_end": ref_hit.end,
        "ref_midpoint": ref_hit.midpoint,
        "ref_CFD": ref_hit.cfd,
        "block_id": block.block_id,
        "block_ref_chr": block.ref_chr,
        "block_ref_start": block.ref_min,
        "block_ref_end": block.ref_max,
        "block_query_chr": block.query_chr,
        "block_query_start": block.query_min,
        "block_query_end": block.query_max,
        "block_orientation": block.orientation,
        "expected_query_chr": block.query_chr,
        "expected_query_pos": expected_query_pos,
        "n_same_block_pid_candidates": n_same_block_pid_candidates,
    }


def empty_query_output(
    base: dict[str, object],
    status: str,
    n_within_tolerance: int,
) -> dict[str, object]:
    return {
        **base,
        "status": status,
        "query_name": pd.NA,
        "query_chr": pd.NA,
        "query_start": pd.NA,
        "query_end": pd.NA,
        "query_midpoint": pd.NA,
        "query_CFD": pd.NA,
        "distance_to_expected": pd.NA,
        "n_within_tolerance": n_within_tolerance,
        "all_candidate_names": pd.NA,
        "all_candidate_midpoints": pd.NA,
    }


def query_only_rows(
    ref_by_block: dict[str, list[AssignedHit]],
    query_by_block: dict[str, list[AssignedHit]],
) -> list[dict[str, object]]:
    """Optional diagnostic rows for query hits in blocks lacking same-ID ref hits."""
    rows: list[dict[str, object]] = []

    for block_id, query_hits in query_by_block.items():
        ref_pids = {assigned.hit.protospacer_id for assigned in ref_by_block.get(block_id, [])}
        for assigned_query in query_hits:
            hit = assigned_query.hit
            if hit.protospacer_id in ref_pids:
                continue
            rows.append(
                {
                    "status": "QUERY_ONLY_IN_BLOCK",
                    "protospacerID": hit.protospacer_id,
                    "ref_name": pd.NA,
                    "ref_chr": pd.NA,
                    "ref_start": pd.NA,
                    "ref_end": pd.NA,
                    "ref_midpoint": pd.NA,
                    "ref_CFD": pd.NA,
                    "block_id": block_id,
                    "query_name": hit.name,
                    "query_chr": hit.chrom,
                    "query_start": hit.start,
                    "query_end": hit.end,
                    "query_midpoint": hit.midpoint,
                    "query_CFD": hit.cfd,
                }
            )

    return rows


# -----------------------------------------------------------------------------
# Public workflow
# -----------------------------------------------------------------------------

def classify_pam_orthology(
    ref_cfd: str | Path,
    query_cfd: str | Path,
    delta: str | Path,
    out_csv: str | Path,
    tol: int = 1000,
    include_query_only: bool = False,
    diagnostics_csv: str | Path | None = None,
) -> pd.DataFrame:
    
    print(f"[INFO] Loading reference CFD: {ref_cfd}")
    ref_df = load_cfd_table(ref_cfd, label="reference")

    print(f"[INFO] Loading query CFD: {query_cfd}")
    query_df = load_cfd_table(query_cfd, label="query")

    print(f"[INFO] Building synteny blocks from delta: {delta}")
    blocks = build_synteny_blocks_from_delta(delta)
    ref_trees, query_trees, block_by_id = build_block_indexes(blocks)
    print(f"[INFO] Loaded {len(blocks):,} synteny blocks")

    print("[INFO] Assigning reference hits to reference-space synteny blocks")
    ref_by_block, ref_diagnostics = assign_hits_to_blocks(ref_df, ref_trees, source="ref")

    print("[INFO] Assigning query hits to query-space synteny blocks")
    query_by_block, query_diagnostics = assign_hits_to_blocks(query_df, query_trees, source="query")

    import pdb
    pdb.set_trace()

    print("[INFO] Comparing hits block-by-block")
    rows: list[dict[str, object]] = []
    all_block_ids = sorted(set(ref_by_block) | set(query_by_block))

    for block_id in all_block_ids:
        ref_hits = ref_by_block.get(block_id, [])
        if not ref_hits:
            continue
        block = block_by_id[block_id]
        query_hits = query_by_block.get(block_id, [])
        rows.extend(compare_block_hits(block, ref_hits, query_hits, tol=tol))

    if include_query_only:
        rows.extend(query_only_rows(ref_by_block, query_by_block))

    out = pd.DataFrame(rows)
    out.to_csv(out_csv, index=False)

    diagnostics = pd.DataFrame(ref_diagnostics + query_diagnostics)
    if diagnostics_csv is not None:
        diagnostics.to_csv(diagnostics_csv, index=False)

    print(f"[INFO] Reference hits: {len(ref_df):,}")
    print(f"[INFO] Query hits:     {len(query_df):,}")
    print(f"[INFO] Assigned ref blocks:   {len(ref_by_block):,}")
    print(f"[INFO] Assigned query blocks: {len(query_by_block):,}")
    print(f"[INFO] Diagnostics rows: {len(diagnostics):,}")
    print(f"[INFO] Wrote {len(out):,} orthology rows -> {out_csv}")

    return out


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Block-first PAM/protospacer orthology classifier using svmu2 synteny.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--ref-cfd", required=True, type=Path, help="Reference/ISO1 CFD CSV")
    parser.add_argument("--query-cfd", required=True, type=Path, help="Query/BL54591 CFD CSV")
    parser.add_argument("--delta", required=True, type=Path, help="nucmer delta file")
    parser.add_argument("--out", required=True, type=Path, help="Output orthology summary CSV")
    parser.add_argument(
        "--diagnostics-out",
        type=Path,
        default=None,
        help="Optional CSV for no-block/ambiguous-block hit diagnostics",
    )
    parser.add_argument(
        "--tol",
        type=int,
        default=1000,
        help="Maximum bp distance between projected ref midpoint and query midpoint",
    )
    parser.add_argument(
        "--include-query-only",
        action="store_true",
        help="Also emit diagnostic rows for query hits in a block lacking same-ID ref hits",
    )

    args = parser.parse_args()

    for path in [args.ref_cfd, args.query_cfd, args.delta]:
        if not path.exists():
            raise FileNotFoundError(path)

    out = classify_pam_orthology(
        ref_cfd=args.ref_cfd,
        query_cfd=args.query_cfd,
        delta=args.delta,
        out_csv=args.out,
        tol=args.tol,
        include_query_only=args.include_query_only,
        diagnostics_csv=args.diagnostics_out,
    )

    if "status" in out.columns:
        print("[INFO] Status counts:")
        print(out["status"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
