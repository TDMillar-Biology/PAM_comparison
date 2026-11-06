# workflow/rules/04_synteny/synteny.smk

rule nucmer_synteny:
    input:
        ref = rules.extract_and_fix_euchromatin.output.fixed_fa,
        qry = rules.extract_main_scaffolds.output.main_fasta
    output:
        delta = "results/04_synteny/ISO1_BL54591.delta"
    log:
        "logs/04_synteny/nucmer_synteny.log"
    container:
        "containers/04_synteny.sif"
    threads: 1
    resources:
        mem_mb = 16000,
        runtime = 120
    shell:
        "nucmer -p results/04_synteny/ISO1_BL54591 {input.ref} {input.qry} > {log} 2>&1"
