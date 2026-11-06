# workflow/rules/02_preprocessing/preprocess.smk

rule extract_and_fix_euchromatin:
    input:
        ref = config["ref_genome"],
        bed = config["euchromatin_bed"]
    output:
        fixed_fa = f"results/02_preprocessing/{config['ref_name']}_euchromatin.fixed.fasta"
    log:
        "logs/02_preprocessing/extract_and_fix_euchromatin.log"
    container:
        "containers/02_preprocessing.sif"
    shell:
        """
        (bedtools getfasta -fi {input.ref} -bed {input.bed} | \
        sed -E 's/^>([^:]+):.*$/>\\1/' > {output.fixed_fa}) > {log} 2>&1
        """

rule repeat_masker:
    input:
        fa = rules.extract_and_fix_euchromatin.output.fixed_fa,
        lib = config["repeat_lib"]
    output:
        masked = f"results/02_preprocessing/{config['ref_name']}_euchromatin.fixed.fasta.masked"
    log:
        "logs/02_preprocessing/repeat_masker.log"
    container:
        "containers/02_preprocessing.sif"
    threads: config.get("threads", 16)
    shell:
        "RepeatMasker -lib {input.lib} -pa {threads} -dir results/02_preprocessing/ {input.fa} > {log} 2>&1"