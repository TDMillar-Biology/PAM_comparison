#!/usr/bin/env python3
'''
Compute CFD scores (Doench et al.) for all protospacer–genome alignments.

This script:
- reads SAM alignments (all mappings preserved)
- extracts strand-aware genomic target sequences via BEDTools
- computes CFD scores per alignment hit
- outputs one row per genomic hit (no collapsing or orthology inference)

Assumes a Linux environment with bedtools available.

last reviewed June 2026 - fully functional, could benefit from modularizing compute_cfd_from_sam
'''
import argparse
import pickle
import pandas as pd
import pysam
import sys
import subprocess
import tempfile
from Bio.Seq import Seq
from collections import defaultdict
from pathlib import Path

############################################################
# Utility functions
############################################################

def revcom(s):
    '''
    accepts string with ATCG dictionary, returns RNA reverse complement of the string
    with Watson Crick base pairing rules
    '''
    basecomp = {'A':'T','C':'G','G':'C','T':'A','U':'A'}
    return ''.join(basecomp[b] for b in s[::-1])

def get_mm_pam_scores(mismatch_scores, pam_scores):
    '''
    Load mismatch and PAM penalty scores from Doench et al. CFD model.
    '''
    mm_scores = pickle.load(open(mismatch_scores, "rb"))
    pam_scores = pickle.load(open(pam_scores, "rb"))
    return mm_scores, pam_scores

def calc_cfd(wt, sg, pam, mm_scores, pam_scores):
    '''
    Compute CFD score (Doench et al.) as a multiplicative penalty model:
    product of per-position mismatch penalties multiplied by PAM penalty.
    '''

    score = 1.0

    wt_u = wt.replace("T", "U")
    sg_u = sg.replace("T", "U")

    for i, (w, s) in enumerate(zip(wt_u, sg_u), start=1):
        if w != s:
            key = f"r{w}:d{revcom(s)},{i}"
            score *= mm_scores[key]

    score *= pam_scores[pam]
    return score

def load_fasta_strip_coords(path):
    """
    Load FASTA created by bedtools -name, stripping ::coords.
    """
    names, seqs = [], []
    current_name = None
    seq_parts = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    names.append(current_name)
                    seqs.append("".join(seq_parts))
                clean = line[1:].split("::")[0]
                current_name = clean
                seq_parts = []
            else:
                seq_parts.append(line)

    if current_name:
        names.append(current_name)
        seqs.append("".join(seq_parts))

    return pd.DataFrame({"name": names, "genome_target": seqs})

############################################################
# Main CFD workflow
############################################################

def compute_cfd_from_sam(samfile, fastafile, outfile, mismatch_scores, pam_scores):
    sam = pysam.AlignmentFile(samfile)

    records = []
    bed_rows = []
    hit_counter = defaultdict(int)

    # -------------------------
    # Extract protospacers + BED
    # -------------------------
    lost_pams = []  ## store unmapped PAMs

    for aln in sam:
        # ---- Handle unmapped reads (lost PAMs) ----
        if aln.is_unmapped: ## this is SAM flag 4
            lost_pams.append({
                "name": aln.query_name,
                "protospacer": aln.query_sequence,
                "note": "unmapped/lost"
            })
            continue

        # Extra safety – malformed coordinate (probably not going to happen if a good aligner is used)
        if aln.reference_start is None or aln.reference_end is None:
            lost_pams.append({
                "name": aln.query_name,
                "protospacer": aln.query_sequence,
                "note": "malformed"
            })
            continue
        
        # maintain an index as to give unique names to alignments
        # this is done to accomodate one to many mapping
        hit_counter[aln.query_name] += 1
        suffix = hit_counter[aln.query_name]
        name = f"{aln.query_name}_{suffix}"

        qseq = aln.query_sequence # this is the sequence we mapped
        if aln.is_reverse: # this is SAM flag 16
            qseq = str(Seq(qseq).reverse_complement())

        # Coordinates of the hit in the reference genome
        ref_chr   = aln.reference_name
        ref_start = aln.reference_start
        ref_end   = aln.reference_end
        ref_mid   = ref_start + (ref_end - ref_start) // 2 # May not be required under new write up, consider removal

        records.append({
            "name"         : aln.query_name + "_" + str(hit_counter[aln.query_name]),
            "protospacer"  : qseq,
            "protospacerID" : aln.query_name,
            # reference genome coordinates for synteny
            "ref_chr"      : ref_chr,
            "ref_start"    : ref_start,
            "ref_end"      : ref_end,
            "ref_midpoint" : ref_mid
        })

        # create bed entry of mapped coord to facilitate reference sequence extraction
        bed_rows.append([
            aln.reference_name,
            aln.reference_start,
            aln.reference_end,
            name,
            0,
            "-" if aln.is_reverse else "+"
        ])

    # -------------------------
    # prep bed out
    # -------------------------

    bed_df = pd.DataFrame(bed_rows)

    ## This is implemented as below to prevent a race condition when run in parallel
    with tempfile.TemporaryDirectory() as tmpdir: 

        tmpdir = Path(tmpdir)

        bed_path = tmpdir / "targets.bed"
        fasta_out = tmpdir / "targets.fasta"
        bed_df.to_csv(
            bed_path,
            sep="\t",
            index=False,
            header=False
        )

        # -------------------------
        # Run bedtools
        # -------------------------
        ## the limitation of this approach is that we lose ambiguity awareness
        ## this means that N can be introduced into sequence info if it's recorded as such in the reference genome
        ## this could be from true ambiguity in ref or Poly N scars from scaffolding the reference

        subprocess.run([
            "bedtools",
            "getfasta",
            "-fi", fastafile,
            "-bed", str(bed_path),
            "-s",
            "-name",
            "-fo", str(fasta_out) 
        ], check=True)


        # -------------------------
        # Load FASTA output
        # -------------------------
        fasta_df = load_fasta_strip_coords(fasta_out)

    # -------------------------
    # Merge protospacers + WT
    # -------------------------
    summary_df = pd.DataFrame(records)
    master = summary_df.merge(fasta_df, on="name")

    # -------------------------
    # Extract WT, PAM, SG for CFD
    # -------------------------
    # ASSUMPTION:
    # genome_target is 23 bp: 20 bp protospacer + NGG PAM
    # indices correspond to SpCas9 targeting (Doench et al.)

    master["wt"]  = master["genome_target"].str[:20]
    master["pam"] = master["genome_target"].str[21:23]
    master["sg"]  = master["protospacer"].str[:20]

    ## skip if N is in target / guide
    wt_N_mask = master["wt"].str.contains("N")
    sg_N_mask = master["sg"].str.contains("N")
    pam_N_mask = master["pam"].str.contains("N")

    print(f"Skipped {wt_N_mask.sum()} ambiguous WT rows.")
    print(f"Skipped {sg_N_mask.sum()} ambiguous SG rows.")
    print(f"Skipped {pam_N_mask.sum()} ambiguous PAM rows.")

    master = master[~wt_N_mask & ~sg_N_mask & ~pam_N_mask].copy()

    mm_scores, pam_scores = get_mm_pam_scores(mismatch_scores, pam_scores)

    master["CFD"] = master.apply(
        lambda r: calc_cfd(r.wt, r.sg, r.pam, mm_scores, pam_scores),
        axis=1
    )

    # -------------------------
    # Save output
    # -------------------------
    master.to_csv(outfile, index=False)
    # ------------------------------------------------
    # Write lost PAMs table
    # ------------------------------------------------
    if lost_pams:
        lost_df = pd.DataFrame(lost_pams)
        lost_out = outfile.replace(".csv", "_lostPAMs.csv")
        lost_df.to_csv(lost_out, index=False)
        print(f"[INFO] Wrote {len(lost_df)} lost PAMs → {lost_out}")

    return master


############################################################
# Run from command line
############################################################

def main():
    parser = argparse.ArgumentParser(
        description="Compute CFD scores (Doench et al.) for CRISPR protospacer-genome alignments.",
        epilog="Example: python compute_CFD.py -s alignments.sam -f reference.fa -o cfd_results.csv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-s", "--sam",
        required=True,
        type=Path,
        help="Path to the input SAM file containing protospacer alignments."
    )
    parser.add_argument(
        "-f", "--fasta",
        required=True,
        type=Path,
        help="Path to the reference genome FASTA file. (Ensure it is indexed with .fai)."
    )
    parser.add_argument(
        "-o", "--out",
        required=True,
        type=Path,
        help="Path to the output CSV file for CFD scores."
    )

    parser.add_argument(
    "--mismatch-scores",
    required=True,
    type=Path,
    help="Path to mismatch_score.pkl"
    )

    parser.add_argument(
    "--pam-scores",
    required=True,
    type=Path,
    help="Path to pam_scores.pkl"
    )

    args = parser.parse_args()

    # Fail early: Check if input files actually exist before spinning up any processes
    if not args.sam.exists():
        sys.exit(f"[FATAL] SAM file not found at: {args.sam}")
    if not args.fasta.exists():
        sys.exit(f"[FATAL] FASTA file not found at: {args.fasta}")
    if not args.mismatch_scores.exists():
        sys.exit(f"[FATAL] mismatch score file not found at: {args.mismatch_scores}")
    if not args.pam_scores.exists():
        sys.exit(f"[FATAL] PAM score file not found at: {args.pam_scores}")

    print(f"[INFO] Initializing CFD calculation...")
    print(f"[INFO] SAM:   {args.sam}")
    print(f"[INFO] FASTA: {args.fasta}")
    print("-" * 40)

    # Execute the workflow (casting Paths back to strings for pysam/bedtools compatibility)
    compute_cfd_from_sam(
        str(args.sam),
        str(args.fasta),
        str(args.out),
        str(args.mismatch_scores),
        str(args.pam_scores)
    )
    
    print(f"\n[INFO] CFD results successfully saved to: {args.out}")

if __name__ == "__main__":
    main()
