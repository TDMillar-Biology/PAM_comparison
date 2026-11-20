#!/usr/bin/env python3
"""Discover ISO1-derived 23-mers with provenance in genomic coordinates.

The input FASTA is now full-length and N-masked *outside* euchromatin.
The .fa guide IDs, .bed source sites, and .tsv manifest all share original
chromosome coordinates. The legacy `*_local` columns remain for compatibility
but are now identical to genomic coordinates (zero-based, half-open).
"""

import re
import argparse
import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

# Overlapping matches: lookahead does not consume the 23-mer, so adjacent PAMs
# are not accidentally skipped. Lowercase masked bases are intentionally rejected.
FWD_PATTERN = re.compile(r'(?=([ACGT]{21}GG))')
REV_PATTERN = re.compile(r'(?=(CC[ACGT]{21}))')


def discover_pams(fasta_path):
    """Yield each 23-mer with genomic coordinates and its orientation."""
    with open(fasta_path) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            chrom = record.id
            seq = str(record.seq)  # Do not .upper(): that would undo soft masking.
            for m in FWD_PATTERN.finditer(seq):
                yield {"chrom": chrom, "start": m.start(1), "end": m.end(1),
                       "sequence": m.group(1), "strand": "+"}
            for m in REV_PATTERN.finditer(seq):
                # On '-' strand the source interval is CCN+20nt; its reverse
                # complement is the usual 20nt-spacer + NGG guide orientation.
                oriented = str(Seq(m.group(1)).reverse_complement())
                yield {"chrom": chrom, "start": m.start(1), "end": m.end(1),
                       "sequence": oriented, "strand": "-"}

def source_offsets(source_bed):
    """Read offsets of intervals used to build the euchromatin FASTA.

    The current extraction rule removes the start offset from FASTA headers.
    Its chromosome names therefore refer to local subsequences, not full
    chromosome coordinates. We require one interval per chromosome to avoid
    silently assigning incorrect genomic coordinates.
    """
    if source_bed is None:
        return {}
    offsets = {}
    with open(source_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            if chrom in offsets:
                raise ValueError(
                    f"More than one interval for {chrom} in {source_bed}; "
                    "the current FASTA-header stripping rule cannot retain "
                    "multiple independent coordinate offsets per chromosome"
                )
            start, end = int(start), int(end)
            if start < 0 or end <= start:
                raise ValueError(f"Invalid source interval: {line.rstrip()}")
            offsets[chrom] = (start, end)
    return offsets


def write_outputs(pam_generator, prefix, offsets=None):
    """Write FASTA, genomic-coordinate BED and a joined provenance TSV."""
    offsets = offsets or {}
    fasta_path, bed_path, manifest_path = (
        f"{prefix}.fa", f"{prefix}.bed", f"{prefix}.tsv"
    )
    fields = ["protospacerID", "source_chr", "source_start_local",
              "source_end_local", "source_start", "source_end",
              "source_strand", "guide_23mer", "guide_20mer", "source_pam"]
    count = 0
    with (open(fasta_path, "w") as fa,
          open(bed_path, "w") as bed,
          open(manifest_path, "w", newline="") as manifest):
        writer = csv.DictWriter(manifest, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for p in pam_generator:
            count += 1
            name = f"PAM_{count}"
            seq = p["sequence"]
            if len(seq) != 23:
                raise AssertionError(f"Unexpected {len(seq)}-nt guide for {name}")
            if offsets and p["chrom"] not in offsets:
                raise ValueError(f"{p['chrom']} missing from source interval BED")
            offset, limit = offsets.get(p["chrom"], (0, None))
            if limit is not None and p["end"] > limit - offset:
                raise ValueError(f"{name} extends beyond {p['chrom']} source interval")
            fa.write(f">{name}\n{seq}\n")
            # Full-length masked FASTA: BED coordinates are now genomic.
            bed.write(f"{p['chrom']}\t{p['start']}\t{p['end']}\t"
                      f"{name}\t0\t{p['strand']}\n")
            writer.writerow({
                "protospacerID": name, "source_chr": p["chrom"],
                "source_start_local": p["start"],
                "source_end_local": p["end"],
                "source_start": offset + p["start"],
                "source_end": offset + p["end"],
                "source_strand": p["strand"],
                "guide_23mer": seq, "guide_20mer": seq[:20],
                "source_pam": seq[21:23],
            })
    print(f"Wrote {fasta_path}")
    print(f"Wrote {bed_path}")
    print(f"Wrote {manifest_path}")
    return count


def main():
    parser = argparse.ArgumentParser(
        description="Discover masked-aware Cas9 PAM-adjacent 23mers")
    parser.add_argument("-i", "--input", required=True,
                        help="Full-length ISO1 FASTA with non-euchromatin hard-masked")
    parser.add_argument("-o", "--outprefix", required=True)
    parser.add_argument("--source-bed", type=Path,
                        help="Legacy option for a BED-extracted FASTA; omit for full-length coordinate-preserving discovery")
    args = parser.parse_args()
    # Without --source-bed, offsets are zero and local/genomic positions match.
    total = write_outputs(discover_pams(args.input), args.outprefix,
                          source_offsets(args.source_bed))
    print(f"[INFO] Discovered {total:,} PAM sites")


if __name__ == "__main__":
    main()
