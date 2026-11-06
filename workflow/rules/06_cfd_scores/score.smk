# workflow/rules/06_cfd_scores/score.smk

rule compute_cfd_iso1:
    input:
        sam = rules.bowtie_align_ref.output.sam,
        ref = rules.extract_and_fix_euchromatin.output.fixed_fa,
        mismatch_scores = config["cfd"]["mismatch_scores"],
        pam_scores = config["cfd"]["pam_scores"]
    output:
        cfd = "results/06_cfd_scores/iso1_raw_cfd.csv"
    log:
        "logs/06_cfd_scores/compute_cfd_iso1.log"
    container:
        "containers/python_analysis.sif"
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
            --pam-scores {input.pam_scores} > {log} 2>&1
        """

rule compute_cfd_bl:
    input:
        sam = rules.bowtie_align_qry.output.sam,
        ref = rules.extract_main_scaffolds.output.main_fasta,
        mismatch_scores = config["cfd"]["mismatch_scores"],
        pam_scores = config["cfd"]["pam_scores"]
    output:
        cfd = "results/06_cfd_scores/bl_raw_cfd.csv"
    log:
        "logs/06_cfd_scores/compute_cfd_bl.log"
    container:
        "containers/python_analysis.sif"
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
            --pam-scores {input.pam_scores} > {log} 2>&1
        """
