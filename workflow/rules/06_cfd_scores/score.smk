# workflow/rules/06_cfd_scores/score.smk

rule compute_cfd_ref:
    input:
        sam = rules.bowtie_align_ref.output.sam,
        ref = rules.extract_and_fix_euchromatin.output.fixed_fa,
        mismatch_scores = config["cfd"]["mismatch_scores"],
        pam_scores = config["cfd"]["pam_scores"]
    output:
        cfd = f"results/06_cfd_scores/{config['ref_name']}_raw_cfd.csv"
    log:
        f"logs/06_cfd_scores/compute_cfd_{config['ref_name']}.log"
    container:
        "workflow/containers/images/python_analysis.sif"
    threads: config.get("threads", 16)
    resources:
        mem_mb = 16000,
        runtime = 60
    shell:
        """
        python3 workflow/scripts/compute_CFD.py \
            --sam {input.sam} \
            --fasta {input.ref} \
            --out {output.cfd} \
            --mismatch-scores {input.mismatch_scores} \
            --pam-scores {input.pam_scores} \
            --threads {threads} > {log} 2>&1
        """

rule compute_cfd_qry:
    input:
        sam = rules.bowtie_align_qry.output.sam,
        ref = rules.extract_main_scaffolds.output.main_fasta,
        mismatch_scores = config["cfd"]["mismatch_scores"],
        pam_scores = config["cfd"]["pam_scores"]
    output:
        cfd = f"results/06_cfd_scores/{config['query_name']}_raw_cfd.csv"
    log:
        f"logs/06_cfd_scores/compute_cfd_{config['query_name']}.log"
    container:
        "workflow/containers/images/python_analysis.sif"
    threads: config.get("threads", 16)
    resources:
        mem_mb = 16000,
        runtime = 60
    shell:
        """
        python3 workflow/scripts/compute_CFD.py \
            --sam {input.sam} \
            --fasta {input.ref} \
            --out {output.cfd} \
            --mismatch-scores {input.mismatch_scores} \
            --pam-scores {input.pam_scores} \
            --threads {threads} > {log} 2>&1
        """