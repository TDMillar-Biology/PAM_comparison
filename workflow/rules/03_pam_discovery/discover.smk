# workflow/rules/03_pam_discovery/discover.smk

rule discover_pams:
    input:
        masked = rules.repeat_masker.output.masked
    output:
        fa = "results/03_pam_discovery/ISO1.masked.euchromatic.PAMs.fa",
        bed = "results/03_pam_discovery/ISO1.masked.euchromatic.PAMs.bed"
    log:
        "logs/03_pam_discovery/discover_pams.log"
    container:
        "containers/python_analysis.sif"
    shell:
        "python3 workflow/scripts/pam_discovery.py -i {input.masked} -o results/03_pam_discovery/ISO1.masked.euchromatic.PAMs > {log} 2>&1"

rule sample_pams:
    input:
        all_pams = rules.discover_pams.output.fa
    output:
        sample = "results/03_pam_discovery/ISO1.masked.euchromatic.PAMs.10pct.fa"
    log:
        "logs/03_pam_discovery/sample_pams.log"
    params:
        percentage = "0.01",
        seed = "123"
    container:
        "containers/01_scaffolding.sif"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 15
    shell:
        "seqtk sample -s {params.seed} {input.all_pams} {params.percentage} > {output.sample} 2> {log}"
