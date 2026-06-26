''' 
AI is not performing well at adapting this script. The results are far too verbose and complex for what should be minimal
Free hand implementation like the good old days before the cyborgs. 
'''

from svmu2.orchestration.parse import parse
from svmu2.orchestration.synteny import resolve_synteny
from svmu2.models.line import build_alignment_primitives ## temporary

import argparse
from pathlib import Path
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="Block-first PAM/protospacer orthology classifier using svmu2 synteny.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--ref-cfd", required=True, type=Path, help="Reference/ISO1 CFD CSV")
    parser.add_argument("--query-cfd", required=True, type=Path, help="Query/BL54591 CFD CSV")
    parser.add_argument("--delta", required=True, type=Path, help="nucmer delta file")
    parser.add_argument("--out", required=True, type=Path, help="Output orthology summary CSV")
    parser.add_argument("--diagnostics-out", type=Path, default=None, help="Optional CSV for no-block/ambiguous-block hit diagnostics")
    parser.add_argument("--tol",type=int,default=1000,help="Maximum bp distance between projected ref midpoint and query midpoint")
    parser.add_argument("--include-query-only",action="store_true",help="Also emit diagnostic rows for query hits in a block lacking same-ID ref hits")

    args = parser.parse_args()

    return args

REQUIRED_CFD_COLUMNS = {
    "protospacerID",
    "ref_chr",
    "ref_start",
    "ref_end",
    "ref_midpoint",
}

def load_cfd_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = REQUIRED_CFD_COLUMNS - set(df.columns)

    if missing:
        raise ValueError(f"{path} CFD table is missing columns: {sorted(missing)}")

    # coerce columns
    for col in ["ref_chr", "protospacerID"]:
        df[col] = df[col].astype(str)
    for col in ["ref_start", "ref_end", "ref_midpoint"]:
        df[col] = df[col].astype(int)

    if "name" not in df.columns:
        df["name"] = df["protospacerID"].astype(str) + "_row" + df.index.astype(str)

    if "CFD" not in df.columns:
        df["CFD"] = pd.NA

    return df

def block_id(block):
    return f"{block.reference}_{block.query}_{block.index}"


def assign_synteny(df, trees):
    df = df.copy()

    block_ids = []
    hit_counts = []

    for row in df.itertuples(index=False):
        chrom = getattr(row, "ref_chr")
        start = getattr(row, "ref_start")
        end = getattr(row, "ref_end")

        tree = trees.get(chrom)

        if tree is None:
            block_ids.append(None)
            hit_counts.append(0)
            continue

        hits = tree[start:end]
        hit_counts.append(len(hits))

        if not hits:
            block_ids.append(None)

        elif len(hits) == 1:
            block = next(iter(hits)).data
            block_ids.append(block_id(block))

        else:
            blocks = [hit.data for hit in hits]
            block_ids.append([block_id(block) for block in blocks])

    df["syntenic_block"] = block_ids
    df["n_hits"] = hit_counts

    return df

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



def main():
    args = parse_args()
    ref_cfd = load_cfd_table(args.ref_cfd)
    qry_cfd = load_cfd_table(args.query_cfd)
    all_alns, primary_alns = parse(args.delta)
    plot(primary=primary_alns)
    resolve_synteny(primary_alignments=primary_alns, breakpoint_map=None)
    # aln objs now have .reference_synteny_tree and .query_synteny_tree attributes

    reference_trees = {aln.reference: aln.reference_synteny_tree for aln in primary_alns.values()}
    query_trees = {aln.query: aln.query_synteny_tree for aln in primary_alns.values()}

    ref_df = assign_synteny(ref_cfd, reference_trees)
    qry_df = assign_synteny(qry_cfd, query_trees)

    import pdb
    pdb.set_trace()

if __name__ == "__main__":
    main()