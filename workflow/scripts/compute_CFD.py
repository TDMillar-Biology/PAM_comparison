#!/usr/bin/env python3
'''
Compute CFD scores (Doench et al.) for all protospacer–genome alignments.

This script:
- Streams SAM alignments to maintain O(1) memory complexity.
- Utilizes multiprocessing.Pool to distribute CFD calculations across CPU cores.
- Workers independently query the reference FASTA to avoid IPC serialization errors.
- Yields output asynchronously via imap_unordered for maximum throughput.
'''
import argparse
import pickle
import pysam
import sys
import csv
import multiprocessing as mp
from Bio.Seq import Seq
from collections import defaultdict
from pathlib import Path

############################################################
# Global Worker Resources
############################################################
# These variables will be initialized once per worker process 
# rather than being passed repeatedly through IPC.
worker_fasta = None
worker_mm_scores = None
worker_pam_scores = None

def worker_init(fasta_path, mismatch_path, pam_path):
    """Initialize read-only objects in each worker's memory space."""
    global worker_fasta, worker_mm_scores, worker_pam_scores
    
    worker_fasta = pysam.FastaFile(fasta_path)
    worker_mm_scores = pickle.load(open(mismatch_path, "rb"))
    worker_pam_scores = pickle.load(open(pam_path, "rb"))

############################################################
# Utility functions
############################################################

def revcom(s):
    basecomp = {'A':'T','C':'G','G':'C','T':'A','U':'A', 'N':'N'}
    return ''.join(basecomp.get(b, 'N') for b in s[::-1])

def calc_cfd(wt, sg, pam):
    """Computes CFD score using worker-level global scoring dicts."""
    score = 1.0
    wt_u = wt.replace("T", "U")
    sg_u = sg.replace("T", "U")

    for i, (w, s) in enumerate(zip(wt_u, sg_u), start=1):
        if w != s:
            key = f"r{w}:d{revcom(s)},{i}"
            score *= worker_mm_scores.get(key, 1.0) 

    score *= worker_pam_scores.get(pam, 1.0)
    return score

############################################################
# Worker Execution
############################################################

def process_chunk(chunk):
    """
    Worker function: takes a chunk of primitive read data, fetches genomic
    context, calculates CFD, and returns formatted rows.
    """
    results = []
    skipped = {'wt': 0, 'sg': 0, 'pam': 0}

    for read in chunk:
        name, qseq, qname, is_reverse, ref_chr, ref_start, ref_end = read

        try:
            genomic_seq = worker_fasta.fetch(ref_chr, ref_start, ref_end).upper()
        except (KeyError, ValueError):
            continue 

        if len(genomic_seq) < 23:
            continue 

        if is_reverse:
            genomic_seq = revcom(genomic_seq)

        wt  = genomic_seq[:20]
        pam = genomic_seq[21:23] 
        sg  = qseq[:20]

        # Tally and skip ambiguities
        if "N" in wt:
            skipped['wt'] += 1
            continue
        if "N" in sg:
            skipped['sg'] += 1
            continue
        if "N" in pam:
            skipped['pam'] += 1
            continue

        cfd = calc_cfd(wt, sg, pam)

        # 1. Calculate number of mismatches between wt and sg
        mm_count = sum(1 for w, s in zip(wt, sg) if w != s)

        # 2. Determine differences between the found PAM and 'GG'
        if pam[0] != 'G':
            mm_count += 1
        if pam[1] != 'G':
            mm_count += 1

        results.append([
            name, qseq, qname, ref_chr,
            ref_start, ref_end,
            genomic_seq, wt, pam, sg, cfd,
            mm_count
        ])

    return results, skipped

############################################################
# Main Dispatcher
############################################################

def compute_cfd_parallel(samfile, fastafile, outfile, mismatch_scores, pam_scores, threads, chunk_size=50000):
    sam = pysam.AlignmentFile(samfile)
    hit_counter = defaultdict(int)
    
    lost_pams = []
    total_skipped = {'wt': 0, 'sg': 0, 'pam': 0}

    # Set up the worker pool
    pool = mp.Pool(
        processes=threads,
        initializer=worker_init,
        initargs=(fastafile, mismatch_scores, pam_scores)
    )

    def read_generator():
        """Generator to pack SAM records into chunks of primitive data."""
        chunk = []
        for aln in sam:
            if aln.is_unmapped or aln.reference_start is None or aln.reference_end is None:
                lost_pams.append({
                    "name": aln.query_name,
                    "protospacer": aln.query_sequence,
                    "note": "unmapped/malformed"
                })
                continue
            
            hit_counter[aln.query_name] += 1
            name = f"{aln.query_name}_{hit_counter[aln.query_name]}"
            
            qseq = aln.query_sequence
            if aln.is_reverse: 
                qseq = str(Seq(qseq).reverse_complement())
                
            # Append pure primitive types (no pysam objects)
            chunk.append((
                name, qseq, aln.query_name, aln.is_reverse, 
                aln.reference_name, aln.reference_start, aln.reference_end
            ))

            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        
        # Yield the final partial chunk
        if chunk:
            yield chunk

    # Stream out results
    with open(outfile, 'w', newline='') as out_f:
        writer = csv.writer(out_f)
        writer.writerow([
            "name", "protospacer", "protospacerID", "ref_chr",
            "ref_start", "ref_end",
            "genome_target", "wt", "pam", "sg", "CFD", 
            "mismatches"
        ])

        # imap_unordered yields results as soon as any worker finishes a chunk
        for chunk_results, skipped_counts in pool.imap_unordered(process_chunk, read_generator()):
            writer.writerows(chunk_results)
            
            # Aggregate skipped metrics
            total_skipped['wt'] += skipped_counts['wt']
            total_skipped['sg'] += skipped_counts['sg']
            total_skipped['pam'] += skipped_counts['pam']

    pool.close()
    pool.join()

    print(f"Skipped {total_skipped['wt']} ambiguous WT rows.")
    print(f"Skipped {total_skipped['sg']} ambiguous SG rows.")
    print(f"Skipped {total_skipped['pam']} ambiguous PAM rows.")

    if lost_pams:
        lost_out = str(outfile).replace(".csv", "_lostPAMs.csv")
        with open(lost_out, 'w', newline='') as f:
            dict_writer = csv.DictWriter(f, fieldnames=["name", "protospacer", "note"])
            dict_writer.writeheader()
            dict_writer.writerows(lost_pams)
        print(f"[INFO] Wrote {len(lost_pams)} lost PAMs → {lost_out}")

############################################################
# Execution
############################################################

def main():
    parser = argparse.ArgumentParser(
        description="Compute CFD scores (Doench et al.) for CRISPR alignments using multiprocessing.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-s", "--sam", required=True, type=Path)
    parser.add_argument("-f", "--fasta", required=True, type=Path)
    parser.add_argument("-o", "--out", required=True, type=Path)
    parser.add_argument("--mismatch-scores", required=True, type=Path)
    parser.add_argument("--pam-scores", required=True, type=Path)
    
    # New parallelization arguments
    parser.add_argument("-t", "--threads", type=int, default=mp.cpu_count(),
                        help="Number of CPU cores to use. Defaults to all available.")
    parser.add_argument("-c", "--chunk-size", type=int, default=50000,
                        help="Number of alignments to pass to a worker at once.")

    args = parser.parse_args()

    if not args.sam.exists(): sys.exit(f"[FATAL] SAM missing: {args.sam}")
    if not args.fasta.exists(): sys.exit(f"[FATAL] FASTA missing: {args.fasta}")

    print(f"[INFO] SAM: {args.sam} | FASTA: {args.fasta}")
    print(f"[INFO] Initiating parallel compute pool across {args.threads} threads...")
    
    compute_cfd_parallel(
        str(args.sam), str(args.fasta), str(args.out),
        str(args.mismatch_scores), str(args.pam_scores),
        args.threads, args.chunk_size
    )

    print(f"\n[INFO] Complete. Output saved to: {args.out}")

if __name__ == "__main__":
    main()