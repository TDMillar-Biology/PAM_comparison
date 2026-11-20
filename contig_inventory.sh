# Adjust filenames if needed
REF="data/reference_genome/ISO1-r6.58_main.fasta"
QRY_FULL="results/01_scaffolding/BL54591.scaffolded.fasta"
QRY_MAIN="results/01_scaffolding/BL54591.canonical.scaffolded.fasta"

mkdir -p results/00_qc/contig_inventory

for entry in \
    "ISO1:$REF" \
    "BL54591_full:$QRY_FULL" \
    "BL54591_main:$QRY_MAIN"
do
    label="${entry%%:*}"
    fasta="${entry#*:}"

    if [[ ! -f "$fasta" ]]; then
        echo "ERROR: Missing $fasta" >&2
        continue
    fi

    # Generate FASTA index if missing or outdated
    if [[ ! -f "${fasta}.fai" || "$fasta" -nt "${fasta}.fai" ]]; then
        samtools faidx "$fasta"
    fi

    # Export contig names and lengths
    cut -f1,2 "${fasta}.fai" \
        > "results/00_qc/contig_inventory/${label}.tsv"

    echo
    echo "=== $label ==="
    awk '
        {n++; total += $2; if ($2 > longest) longest = $2}
        END {
            printf "Contigs: %d\n", n
            printf "Total:   %.3f Mb\n", total/1e6
            printf "Longest: %.3f Mb\n", longest/1e6
        }
    ' "${fasta}.fai"

    echo "Contig inventory:"
    cat "results/00_qc/contig_inventory/${label}.tsv"
done
