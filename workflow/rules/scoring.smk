rule compute_cfd_iso1:
    input:
        sam = "results/05_alignment/ISO1_pams.sam",
        ref = "results/02_preprocessing/ISO1-r6.58_euchromatin.fixed.fasta",
        mismatch_scores = config["cfd"]["mismatch_scores"],
        pam_scores = config["cfd"]["pam_scores"]
    output:
        cfd = "results/06_cfd_scores/iso1_raw_cfd.csv"
    conda:
        "../envs/python_bio.yaml"
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
            --pam-scores {input.pam_scores}
        """

rule compute_cfd_bl:
    input:
        sam = "results/05_alignment/BL54591_pams.sam",
        ref = "results/01_scaffolding/BL54591.canonical.scaffolded.fasta",
        mismatch_scores = config["cfd"]["mismatch_scores"],
        pam_scores = config["cfd"]["pam_scores"]
    output:
        cfd = "results/06_cfd_scores/bl_raw_cfd.csv"
    conda:
        "../envs/python_bio.yaml"
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
            --pam-scores {input.pam_scores}
        """

rule map_pam_orthology:
    input:
        ref_cfd = "results/06_cfd_scores/iso1_raw_cfd.csv",
        qry_cfd = "results/06_cfd_scores/bl_raw_cfd.csv",
        delta = "results/04_synteny/ISO1_BL54591.delta"
    output:
        conserved = "results/07_summary/conserved_pams.csv",
        failed = "results/07_summary/failed_tolerance_pams.csv",
        unmapped_ref = "results/07_summary/unmapped_ref.csv",
        unmapped_qry = "results/07_summary/unmapped_qry.csv"
    params:
        out_dir = "results/07_summary",
        tol = 1000
    conda:
        "../envs/PAM_orthology.yaml"
    shell:
        """
        python3 workflow/scripts/synteny_mapper.py \
            --ref-cfd {input.ref_cfd} \
            --query-cfd {input.qry_cfd} \
            --delta {input.delta} \
            --out-dir {params.out_dir} \
            --tol {params.tol}
        """

rule plot_synteny_diagnostics:
    input:
        # Tying these as inputs guarantees this won't run until the mapper finishes successfully
        conserved = "results/07_summary/conserved_pams.csv",
        failed = "results/07_summary/failed_tolerance_pams.csv",
        unmapped_ref = "results/07_summary/unmapped_ref.csv",
        unmapped_qry = "results/07_summary/unmapped_qry.csv"
    output:
        figures = directory("results/08_figures")
    params:
        data_dir = "results/07_summary"
    conda:
        "../envs/PAM_orthology.yaml"
    shell:
        """
        python3 workflow/scripts/synteny_plotter.py \
            --data-dir {params.data_dir} \
            --figures {output.figures}
        """