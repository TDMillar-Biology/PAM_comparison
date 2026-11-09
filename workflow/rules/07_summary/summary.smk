# workflow/rules/07_summary/summary.smk

rule map_pam_orthology:
    input:
        ref_cfd = rules.compute_cfd_ref.output.cfd,
        qry_cfd = rules.compute_cfd_qry.output.cfd,
        delta   = rules.nucmer_synteny.output.delta
    output:
        conserved    = f"results/07_summary/{config['ref_name']}_{config['query_name']}/conserved_pams.csv",
        failed       = f"results/07_summary/{config['ref_name']}_{config['query_name']}/failed_tolerance_pams.csv",
        unmapped_ref = f"results/07_summary/{config['ref_name']}_{config['query_name']}/unmapped_ref.csv",
        unmapped_qry = f"results/07_summary/{config['ref_name']}_{config['query_name']}/unmapped_qry.csv"
    params:
        out_dir = f"results/07_summary/{config['ref_name']}_{config['query_name']}",
        tol     = config.get("cfd", {}).get("tolerance", 1000)
    log:
        f"logs/07_summary/map_pam_orthology_{config['ref_name']}_{config['query_name']}.log"
    container:
        "containers/images/python_analysis.sif"
    threads: config.get("threads", 8)
    resources:
        mem_mb = 16000,
        runtime = 60
    shell:
        """
        python3 workflow/scripts/synteny_mapper.py \
            --ref-cfd {input.ref_cfd} \
            --query-cfd {input.qry_cfd} \
            --delta {input.delta} \
            --out-dir {params.out_dir} \
            --tol {params.tol} > {log} 2>&1
        """