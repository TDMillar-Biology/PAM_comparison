#!/usr/bin/env python3
"""
January 20th 2026 refactor
Parallel annotation of CFD rows with synteny blocks
IMPORTANT -- fixed serialization of tree objects that was extremely IO intensive and slowed the parallelization drastically
"""

import argparse
import pandas as pd
from multiprocessing import Pool, cpu_count

from svmu2.orchestration.synteny import resolve_synteny
from svmu2.orchestration.parse import parse

############################################################
# Worker globals (set once per process)
############################################################
_SYNTENY_TREES = None

def init_worker(trees):
    """
    Initializer for worker processes.
    Each worker gets its own deserialized copy of the synteny trees.
    """
    global _SYNTENY_TREES
    _SYNTENY_TREES = trees

def lookup_row(args):
    """
    Look up a single CFD row against synteny interval trees.

    args: (ref_chr, ref_pos)
    """
    ref_chr, ref_pos = args

    tree = _SYNTENY_TREES.get(ref_chr)
    if tree is None:
        return None

    hits = tree.at(ref_pos)
    if not hits:
        return None

    # If multiple blocks overlap, return the block_id of the first
    # (or customize selection logic here)
    best = max(hits, key=lambda iv: (iv.end - iv.begin, -iv.data.index)) # longest with overlap
    return best.data.index

############################################################
# Main pipeline
############################################################
def annotate_cfd_parallel(cfd_file, delta_file, out_file, nproc):
    print(f"[INFO] Load CFD table: {cfd_file}")
    df = pd.read_csv(cfd_file)

    _, alns = parse(delta_file)
    primary = resolve_synteny(alns)

    trees = {}
    for ref_chrom, aln in primary.items():
        print(f"[INFO] Processing alignment: {ref_chrom}")

        if getattr(aln, "primary_synteny_tree", None) is not None:
            trees[ref_chrom] = aln.primary_synteny_tree

    if not trees:
        print("[WARN] No usable synteny trees built — exiting")
        return

    print(f"[INFO] Parallel synteny lookup for {len(df)} rows")
    print(f"[INFO] Using {nproc} worker processes")

    args = [
        (df.loc[i, "ref_chr"], df.loc[i, "ref_midpoint"])
        for i in range(len(df))
    ]

    with Pool(
        processes=nproc,
        initializer=init_worker,
        initargs=(trees,),
    ) as pool:
        result = pool.map(lookup_row, args)

    df["block_id"] = result

    df.to_csv(out_file, index=False)
    print(f"[✔] Saved output → {out_file}")

    return df


############################################################
# CLI facing main
############################################################
def main():
    parser = argparse.ArgumentParser(
        description="Annotate CFD table with synteny block IDs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-c", "--cfd",
        required=True,
        help="Input CFD CSV file"
    )
    parser.add_argument(
        "-d", "--delta",
        required=True,
        help="Delta alignment file"
    )
    parser.add_argument(
        "-o", "--out",
        required=True,
        help="Output CSV file"
    )
    parser.add_argument(
        "--nproc",
        type=int,
        default=None,
        help="Number of worker processes to use (default: safe auto)"
    )

    args = parser.parse_args()

    ### SAFE NPROC HANDLING THAT MAY BE OVERKILL BUT WORKS WELL AT STOPPING MISTAKES ON SHARED COMPUTE NODES
    avail = cpu_count()

    if args.nproc is None:
        # Conservative default
        nproc = max(1, min(4, avail - 1))
        print(f"[INFO] --nproc not specified, defaulting to {nproc}")
    else:
        if args.nproc < 1:
            raise ValueError("--nproc must be >= 1")
        if args.nproc > avail:
            print(
                f"[WARN] Requested nproc={args.nproc} exceeds available CPUs ({avail}); "
                f"capping to {avail - 1}"
            )
            nproc = max(1, avail - 1)
        else:
            nproc = args.nproc

    ### END OF SAFE NPROC HANDLING

    annotate_cfd_parallel(
        cfd_file=args.cfd,
        delta_file=args.delta,
        out_file=args.out,
        nproc=nproc
    )

if __name__ == "__main__":
    main()
