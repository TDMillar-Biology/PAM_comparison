rule bowtie_index_ref:
    input:
        ref = "results/02_preprocessing/ISO1-r6.58_euchromatin.fixed.fasta"
    output:
        idx = "results/05_alignment/ISO1/ISO1.1.ebwt" 
    params:
        prefix = "results/05_alignment/ISO1/ISO1"
    conda: 
        "../envs/alignment.yaml"
    shell:
        "bowtie-build {input.ref} {params.prefix}"

rule bowtie_align_ref:
    input:
        idx = "results/05_alignment/ISO1/ISO1.1.ebwt",
        pams = get_pam_fasta
    output:
        sam = "results/05_alignment/ISO1_pams.sam"
    log:
        "logs/05_alignment/bowtie_iso1.log"
    params:
        prefix = "results/05_alignment/ISO1/ISO1",
        mismatch = config["bowtie_mismatches"],
        max_hits = config["max_hits"]
    conda: 
        "../envs/alignment.yaml"
    threads: config["threads"]
    resources:
        mem_mb = 32000,   
        runtime = 360    
    shell:
        "bowtie -f -m {params.max_hits} -v {params.mismatch} -a --best --sam -p {threads} {params.prefix} {input.pams} > {output.sam} 2> {log}"

rule bowtie_index_qry:
    input:
        ref = "results/01_scaffolding/BL54591.canonical.scaffolded.fasta"
    output:
        idx = "results/05_alignment/BL54591/BL54591.1.ebwt" 
    params:
        prefix = "results/05_alignment/BL54591/BL54591"
    conda: 
        "../envs/alignment.yaml"
    shell:
        "bowtie-build {input.ref} {params.prefix}"

rule bowtie_align_qry:
    input:
        idx = "results/05_alignment/BL54591/BL54591.1.ebwt",
        pams = get_pam_fasta
    output:
        sam = "results/05_alignment/BL54591_pams.sam"
    log:
        "logs/05_alignment/bowtie_bl.log"
    params:
        prefix = "results/05_alignment/BL54591/BL54591",
        mismatch = config["bowtie_mismatches"],
        max_hits = config["max_hits"]
    conda: 
        "../envs/alignment.yaml"
    threads: config["threads"]
    resources:
        mem_mb = 32000,
        runtime = 360
    shell:
        "bowtie -f -m {params.max_hits} -v {params.mismatch} -a --best --sam -p {threads} {params.prefix} {input.pams} > {output.sam} 2> {log}"