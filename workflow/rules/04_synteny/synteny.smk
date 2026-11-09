# workflow/rules/04_synteny/synteny.smk

rule nucmer_synteny:
    input:
        ref = rules.extract_and_fix_euchromatin.output.fixed_fa,
        qry = rules.extract_main_scaffolds.output.main_fasta
    output:
        delta = f"results/04_synteny/{config['ref_name']}_{config['query_name']}.delta"
    params:
        prefix = f"results/04_synteny/{config['ref_name']}_{config['query_name']}"
    log:
        "logs/04_synteny/nucmer_synteny.log"
    container:
        "containers/images/04_synteny.sif"
    threads: config.get("threads", 8)
    resources:
        mem_mb = 16000,
        runtime = 120
    shell:
        "nucmer -p {params.prefix} {input.ref} {input.qry} -t {threads} > {log} 2>&1"
