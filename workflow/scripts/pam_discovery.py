#!/usr/bin/env python3

import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

# ----------------------------
# PAM patterns (masked-aware)
# ----------------------------

# Finds overlapping 23-mers (20nt spacer + NGG) excluding all masked regions
FWD_PATTERN = re.compile(r'(?=([ACGT]{21}GG))')

# Finds overlapping 23-mers on the reverse strand (CCN + 20nt spacer)
REV_PATTERN = re.compile(r'(?=(CC[ACGT]{21}))')


def discover_pams(fasta_path):
    """
    PAM generator object read fasta stream output
    Yield PAM records with coordinates and strand.
    """
    for record in SeqIO.parse(fasta_path, "fasta"):
        chrom = record.id
        seq = str(record.seq) # don't upper to stay robust to soft masking

        # Forward strand NGG
        for m in FWD_PATTERN.finditer(seq):
            yield {
                "chrom": chrom,
                "start": m.start(1), # start index of sub string in greater string
                "end": m.end(1), # end index
                "sequence": m.group(1), # sequence matching the RE pattern (ie the 23mer)
                "strand": "+"
            }

        # Reverse strand CCN
        for m in REV_PATTERN.finditer(seq):
            # Reverse complement to orient as 20nt + NGG
            pam_seq = str(Seq(m.group(1)).reverse_complement())
            yield {
                "chrom": chrom,
                "start": m.start(1),
                "end": m.end(1),
                "sequence": pam_seq,
                "strand": "-"
            }


def write_outputs(pam_generator, prefix):
    """
    Write FASTA and BED outputs streaming directly from the generator.
    """
    fasta_path = f"{prefix}.fa"
    bed_path = f"{prefix}.bed"
    count = 0

    with open(fasta_path, "w") as fa, open(bed_path, "w") as bed:
        for p in pam_generator:
            count += 1
            name = f"PAM_{count}"

            # FASTA
            fa.write(f">{name}\n{p['sequence']}\n")

            # BED
            bed.write(
                f"{p['chrom']}\t{p['start']}\t{p['end']}\t"
                f"{name}\t0\t{p['strand']}\n"
            )

    print(f"Wrote {fasta_path}")
    print(f"Wrote {bed_path}")
    return count


def main():
    parser = argparse.ArgumentParser(
        description="Discover masked-aware Cas9 PAM-adjacent 23mers"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Masked FASTA (e.g. ISO1 euchromatin)"
    )
    parser.add_argument(
        "-o", "--outprefix",
        required=True,
        help="Output prefix. Outputs will be prefix.bed, prefix.fa"
    )

    args = parser.parse_args()

    # Pass the generator directly to the writer to save RAM
    pam_generator = discover_pams(args.input)
    total_found = write_outputs(pam_generator, args.outprefix)
    
    print(f"[INFO] Discovered {total_found:,} PAM sites")

if __name__ == "__main__":
    main()