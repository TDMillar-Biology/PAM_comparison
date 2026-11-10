# workflow/rules/03_pam_discovery/discover.smk

rule discover_pams:
    input:
        masked = rules.repeat_masker.output.masked
    output:
        fa  = f"results/03_pam_discovery/{config['ref_name']}.masked.euchromatic.PAMs.fa",
        bed = f"results/03_pam_discovery/{config['ref_name']}.masked.euchromatic.PAMs.bed"
    params:
        out_prefix = f"results/03_pam_discovery/{config['ref_name']}.masked.euchromatic.PAMs"
    log:
        "logs/03_pam_discovery/discover_pams.log"
    container:
        "workflow/containers/images/python_analysis.sif"
    shell:
        """
        python3 workflow/scripts/pam_discovery.py \
            -i {input.masked} \
            -o {params.out_prefix} > {log} 2>&1
        """

rule sample_pams:
    input:
        all_pams = rules.discover_pams.output.fa
    output:
        sample = f"results/03_pam_discovery/{config['ref_name']}.masked.euchromatic.PAMs.10pct.fa"
    params:
        fraction = config.get("test_fraction", 0.01),
        seed = config.get("seed", 123)
    log:
        "logs/03_pam_discovery/sample_pams.log"
    container:
        "workflow/containers/images/01_scaffolding.sif"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 15
    shell:
        "seqtk sample -s {params.seed} {input.all_pams} {params.fraction} > {output.sample} 2> {log}"