#!/bin/bash
#SBATCH --job-name=pam_bowtie
#SBATCH --array=0-3
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/pam_%A_%a.out
#SBATCH --error=logs/pam_%A_%a.err
#SBATCH --ntasks-per-node=1
set -euo pipefail

##############################
# Inputs
##############################
ISO1_IDX="ISO1"
BL_IDX="BL54591"

PAMS_ALL="ISO1.masked.euchromatic.PAMs.fa"
PAMS_SUB="ISO1.masked.euchromatic.PAMs.10pct.fa"

module load GCC/13.3.0 Bowtie/1.3.1

##############################
# Decode array index
##############################
case "$SLURM_ARRAY_TASK_ID" in
  0)
    GENOME="ISO1"
    IDX="$ISO1_IDX"
    READS="$PAMS_ALL"
    OUT="pams_ISO1_all.sam"
    ;;
  1)
    GENOME="ISO1"
    IDX="$ISO1_IDX"
    READS="$PAMS_SUB"
    OUT="pams_ISO1_10pct.sam"
    ;;
  2)
    GENOME="BL54591"
    IDX="$BL_IDX"
    READS="$PAMS_ALL"
    OUT="pams_BL54591_all.sam"
    ;;
  3)
    GENOME="BL54591"
    IDX="$BL_IDX"
    READS="$PAMS_SUB"
    OUT="pams_BL54591_10pct.sam"
    ;;
  *)
    echo "Invalid SLURM_ARRAY_TASK_ID"
    exit 1
    ;;
esac

##############################
# Run Bowtie
##############################
echo "[INFO] Genome: $GENOME"
echo "[INFO] Reads:  $READS"
echo "[INFO] Output: $OUT"

bowtie -f -v 3 -a --best --sam "$IDX" "$READS" > "$OUT"

echo "[INFO] Done: $OUT"
