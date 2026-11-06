# workflow/rules/07_summary/summary.smk

rule map_pam_orthology:
    input:
        ref_cfd = rules.compute_cfd_iso1.output.cfd,
        qry_cfd = rules.compute_cfd_bl.output.cfd,
        delta = rules.nucmer_synteny.output.delta
    output:
        conserved = "results/07_summary/conserved_pams.csv",
        failed = "results/07_summary/failed_tolerance_pams.csv",
        unmapped_ref = "results/07_summary/unmapped_ref.csv",
        unmapped_qry = "results/07_summary/unmapped_qry.csv"
    log:
        "logs/07_summary/map_pam_orthology.log"
    params:
        out_dir = "results/07_summary",
        tol = 1000
    container:
        "containers/python_analysis.sif"
    shell:
        """
        python3 workflow/scripts/synteny_mapper.py \
            --ref-cfd {input.ref_cfd} \
            --query-cfd {input.qry_cfd} \
            --delta {input.delta} \
            --out-dir {params.out_dir} \
            --tol {params.tol} > {log} 2>&1
        """
