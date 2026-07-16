rule classify_pam_orthology: ##candidate for removal
    input:
        ref_cfd = "results/06_cfd_scores/iso1_raw_cfd.csv",
        qry_cfd = "results/06_cfd_scores/bl_raw_cfd.csv",
        delta = "results/04_synteny/ISO1_BL54591.delta"
    output:
        summary = "results/07_summary/pam_orthology_summary.csv",
        figures = directory("results/08_figures")
    conda:
        "../envs/PAM_orthology.yaml"
    resources:
        mem_mb = 16000,
        runtime = 60
    shell:
        """
        python3 workflow/scripts/PAM_orthology_candidate.py \
            --ref-cfd {input.ref_cfd} \
            --query-cfd {input.qry_cfd} \
            --delta {input.delta} \
            --out {output.summary} \
            --figures {output.figures} \
            --tol 10_000
        """