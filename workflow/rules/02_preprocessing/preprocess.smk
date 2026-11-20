# workflow/rules/02_preprocessing/preprocess.smk

## Keep the original ISO1 chromosome coordinates during guide discovery.
## The BED lists regions to KEEP; complement identifies bases to hard-mask.
## Unlike bedtools getfasta, maskfasta never changes sequence lengths/offsets.
rule mask_non_euchromatin:
    input:
        ref = rules.prepare_canonical_reference.output.fasta,
        bed = config["euchromatin_bed"]
    output:
        masked_fa = f"results/02_preprocessing/{config['ref_name']}.euchromatin_hardmasked.fasta",
        excluded_bed = f"results/02_preprocessing/{config['ref_name']}.non_euchromatin.bed",
        genome = f"results/02_preprocessing/{config['ref_name']}.genome"
    log:
        "logs/02_preprocessing/mask_non_euchromatin.log"
    container:
        "workflow/containers/images/02_preprocessing.sif"
    shell:
        """
        set -euo pipefail
        # The FASTA index supplies chromosome lengths and order for BEDTools.
        samtools faidx {input.ref:q} > {log:q} 2>&1
        cut -f1,2 {input.ref:q}.fai > {output.genome:q}
        bedtools sort -i {input.bed:q} -g {output.genome:q} | \
            bedtools complement -i stdin -g {output.genome:q} \
            > {output.excluded_bed:q} 2>> {log:q}
        bedtools maskfasta -fi {input.ref:q} \
            -bed {output.excluded_bed:q} \
            -fo {output.masked_fa:q} >> {log:q} 2>&1
        """

rule repeat_masker:
    input:
        fa = rules.mask_non_euchromatin.output.masked_fa,
        lib = config["repeat_lib"]
    output:
        masked = f"results/02_preprocessing/{config['ref_name']}.euchromatin_hardmasked.fasta.masked"
    log:
        "logs/02_preprocessing/repeat_masker.log"
    container:
        "workflow/containers/images/02_preprocessing.sif"
    threads: config.get("threads", 16)
    shell:
        """
        pa=$(( {threads} / 4 ))
        if [ "$pa" -lt 1 ]; then pa=1; fi

        FAMDB_DIR="" FAMDB_DATA_DIR="" \
        RepeatMasker \
            -lib {input.lib} \
            -pa "$pa" \
            -dir results/02_preprocessing/ \
            {input.fa} > {log} 2>&1
        """