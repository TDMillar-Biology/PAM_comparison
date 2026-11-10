# workflow/rules/05_alignment/align.smk

rule bowtie_index_ref:
    input:
        ref = rules.extract_and_fix_euchromatin.output.fixed_fa
    output:
        idx = f"results/05_alignment/{config['ref_name']}/{config['ref_name']}.1.ebwt"
    log:
        "logs/05_alignment/bowtie_index_ref.log"
    params:
        prefix = f"results/05_alignment/{config['ref_name']}/{config['ref_name']}"
    container:
        "workflow/containers/images/05_alignment.sif"
    shell:
        "bowtie-build {input.ref} {params.prefix} > {log} 2>&1"

rule bowtie_align_ref:
    input:
        idx = rules.bowtie_index_ref.output.idx,
        pams = get_pam_fasta
    output:
        sam = f"results/05_alignment/{config['ref_name']}_pams.sam"
    log:
        "logs/05_alignment/bowtie_align_ref.log"
    params:
        prefix = f"results/05_alignment/{config['ref_name']}/{config['ref_name']}",
        mismatch = config["bowtie_mismatches"],
        max_hits = config["max_hits"]
    container:
        "workflow/containers/images/05_alignment.sif"
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
        idx = f"results/05_alignment/{config['query_name']}/{config['query_name']}.1.ebwt"
    log:
        "logs/05_alignment/bowtie_index_qry.log"
    params:
        prefix = f"results/05_alignment/{config['query_name']}/{config['query_name']}"
    container:
        "workflow/containers/images/05_alignment.sif"
    shell:
        "bowtie-build {input.ref} {params.prefix} > {log} 2>&1"

rule bowtie_align_qry:
    input:
        idx = rules.bowtie_index_qry.output.idx,
        pams = get_pam_fasta
    output:
        sam = f"results/05_alignment/{config['query_name']}_pams.sam"
    log:
        "logs/05_alignment/bowtie_align_qry.log"
    params:
        prefix = f"results/05_alignment/{config['query_name']}/{config['query_name']}",
        mismatch = config["bowtie_mismatches"],
        max_hits = config["max_hits"]
    container:
        "workflow/containers/images/05_alignment.sif"
    threads: config.get("threads", 16)
    resources:
        mem_mb = 32000,
        runtime = 360
    shell:
        "bowtie -f -m {params.max_hits} -v {params.mismatch} -a --best --sam -p {threads} {params.prefix} {input.pams} > {output.sam} 2> {log}"