#!/usr/bin/env bash
set -euo pipefail

############################################################
# CRISPR–CFD + Synteny Analysis Pipeline
# Author: Trevor Millar
############################################################

##############################
# Variables
##############################
threads=16
ref="../../ref_genome/ISO1-r6.58_main.fasta" # ISO1 release 6 only main chromosomes
qry="../../assemblies/final/BL54591.curated.broken.fasta"
qry_scaf="./ragtag_output/ragtag.scaffold.fasta"
qry_scaf_main="BL54591.canonical.scaffolded.fasta" # scaffolded main contigs (subset of fasta from qry_scaf)

repeat_masker_lib="../../repeat_masker_libraries/LarracuenteLab.Repeat.library.specieslib_mod2_Mel_042623.fasta"

ISO1_FA="ISO1-r6.58_euchromatin.fasta.masked"
BL_FA="BL54591.curated.broken.fasta.masked" # no longer necessary (map to full scaffolded de novo genome)
PAM_FASTA="masked_ISO1-r6.58_PAM.fa"
DELTA="ISO1_BL54591.delta"

# Output files
ISO1_SAM="pams_ISO1_10pct.sam"
BL_SAM="pams_BL54591_10pct.sam"
ISO1_CFD="iso1_raw_cfd.csv"
BL_CFD="bl_raw_cfd.csv"
ISO1_SYNT="iso1_synteny.csv"
BL_SYNT="bl_synteny.csv"
FINAL="pam_summary.csv"


############################################################
# Preprocessing steps (Run once)
############################################################
# These steps produce the masked FASTA files and the delta alignment.
# Comment these out after initial preprocessing.
echo "scaffolding BL54591 to ISO1 r6.58 main contigs prior to masking or euchromatin sub setting"
#ragtag.py scaffold "$ref" "$qry"

if [[ ! -s "$qry_scaf" ]]; then
    echo "Error: RagTag scaffolding failed or output not found: $qry_scaf"
    exit 1
fi

# restrict the scaffolded de novo genome to the canonical contigs (same as ISO1_main)
#awk '$1 ~ /_RagTag$/ {print $1}' ./ragtag_output/ragtag.scaffold.agp | sort -u > main_scaffolds.txt
#seqtk subseq "$qry_scaf" main_scaffolds.txt > "$qry_scaf_main"


echo "0a. Extracting euchromatin regions from ISO1"
#bedtools getfasta -fi "$ref" -bed euchromatin_boundaries.bed -fo ISO1-r6.58_euchromatin.fasta

echo "0b. Fixing chromosome names"
#sed -i s'/X:277911-18930000/X/'g  ISO1-r6.58_euchromatin.fasta
#sed -i s'/2L:82455-19570000/2L/'g ISO1-r6.58_euchromatin.fasta
#sed -i s'/2R:8860000-24684540/2R/'g ISO1-r6.58_euchromatin.fasta
#sed -i s'/3L:158639-18438500/3L/'g ISO1-r6.58_euchromatin.fasta
#sed -i s'/3R:9497000-31845060/3R/'g ISO1-r6.58_euchromatin.fasta

echo "0c. RepeatMasking ISO1 prior to PAM discovery"
#RepeatMasker -lib "$repeat_masker_lib" -pa "$threads" ISO1-r6.58_euchromatin.fasta
#python3 pam_discovery.py -i ISO1-r6.58_euchromatin.fasta.masked -o ISO1.masked.euchromatic.PAMs

echo "0d. Running nucmer for synteny"
#nucmer -p ISO1_BL54591 "$ref" "$qry_scaf_main"

############################################################
# Stage 1 — Bowtie indexing
############################################################
echo "1. Building Bowtie indexes"

#bowtie-build "$ref" ISO1
#bowtie-build "$qry_scaf_main" BL54591


############################################################
# Stage 2 — Align all protospacers
############################################################
echo "2. Aligning protospacers to both genomes (v=3 mismatches)"

#bowtie -f -v 3 -a --best --sam ISO1 ISO1.masked.euchromatic.PAMs.fa > "$ISO1_SAM"
#bowtie -f -v 3 -a --best --sam BL54591 ISO1.masked.euchromatic.PAMs.fa > "$BL_SAM"


############################################################
# Stage 3 — Compute CFD scores for every mapping
############################################################
echo "3. Computing CFD scores with compute_CFD.py"

#python3 compute_CFD.py "$ISO1_SAM" "$ref" "$ISO1_CFD"
#python3 compute_CFD.py "$BL_SAM"   "$qry_scaf_main"   "$BL_CFD"


############################################################
# Stage 4 — Annotate synteny blocks
############################################################
echo "4. Annotating synteny with annotate_synteny.py"

python3 annotate_synteny.py "$ISO1_CFD" "$DELTA" "$ISO1_SYNT"
python3 annotate_synteny.py "$BL_CFD"   "$DELTA" "$BL_SYNT"


############################################################
# Stage 5 — Merge ISO1 + BL54591 CFD tables by synteny
############################################################
echo "5. Merging synteny-corrected tables"

#python3 merge_synteny_corrected.py "$ISO1_SYNT" "$BL_SYNT" "$FINAL"


############################################################
# Done
############################################################
echo  "Pipeline complete. Final output: $FINAL"
