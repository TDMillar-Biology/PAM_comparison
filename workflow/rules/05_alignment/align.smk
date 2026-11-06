# workflow/rules/05_alignment/align.smk

rule bowtie_index_ref:
    input:
        ref = rules.extract_and_fix_euchromatin.output.fixed_fa
    output:
        idx = "results/05_alignment/ISO1/ISO1.1.ebwt"
    log:
        "logs/05_alignment/bowtie_index_ref.log"
    params:
        prefix = "results/05_alignment/ISO1/ISO1"
    container:
        "containers/05_alignment.sif"
    shell:
        "bowtie-build {input.ref} {params.prefix} > {log} 2>&1"

rule bowtie_align_ref:
    input:
        idx = rules.bowtie_index_ref.output.idx,
        pams = get_pam_fasta
    output:
        sam = "results/05_alignment/ISO1_pams.sam"
    log:
        "logs/05_alignment/bowtie_align_ref.log"
    params:
        prefix = "results/05_alignment/ISO1/ISO1",
        mismatch = config["bowtie_mismatches"],
        max_hits = config["max_hits"]
    container:
        "containers/05_alignment.sif"
    threads: config.get("threads", 16)
    resources:
        mem_mb = 32000,
        runtime = 360
    shell:
        "bowtie -f -m {params.max_hits} -v {params.mismatch} -a --best --sam -p {threads} {params.prefix} {input.pams} > {output.sam} 2> {log}"

rule bowtie_index_qry:
    input:
        ref = rules.extract_main_scaffolds.output.main_fasta
    output:
        idx = "results/05_alignment/BL54591/BL54591.1.ebwt"
    log:
        "logs/05_alignment/bowtie_index_qry.log"
    params:
        prefix = "results/05_alignment/BL54591/BL54591"
    container:
        "containers/05_alignment.sif"
    shell:
        "bowtie-build {input.ref} {params.prefix} > {log} 2>&1"

rule bowtie_align_qry:
    input:
        idx = rules.bowtie_index_qry.output.idx,
        pams = get_pam_fasta
    output:
        sam = "results/05_alignment/BL54591_pams.sam"
    log:
        "logs/05_alignment/bowtie_align_qry.log"
    params:
        prefix = "results/05_alignment/BL54591/BL54591",
        mismatch = config["bowtie_mismatches"],
        max_hits = config["max_hits"]
    container:
        "containers/05_alignment.sif"
    threads: config.get("threads", 16)
    resources:
        mem_mb = 32000,
        runtime = 360
    shell:
        "bowtie -f -m {params.max_hits} -v {params.mismatch} -a --best --sam -p {threads} {params.prefix} {input.pams} > {output.sam} 2> {log}"
