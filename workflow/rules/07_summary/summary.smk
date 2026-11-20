# workflow/rules/07_summary/summary.smk


rule map_pam_orthology:
    input:
        ref_cfd = rules.compute_cfd_ref.output.cfd,
        qry_cfd = rules.compute_cfd_qry.output.cfd,
        delta   = rules.nucmer_synteny.output.delta,
        # Each guide ID is tied to its original ISO1 locus and PAM here.
        source_manifest = rules.discover_pams.output.manifest
    output:
        conserved = f"results/07_summary/{config['ref_name']}_{config['query_name']}/conserved_pams.csv",
        failed = f"results/07_summary/{config['ref_name']}_{config['query_name']}/failed_tolerance_pams.csv",
        unmapped_ref = f"results/07_summary/{config['ref_name']}_{config['query_name']}/unmapped_ref.csv",
        unmapped_qry = f"results/07_summary/{config['ref_name']}_{config['query_name']}/unmapped_qry.csv",
        tolerance_summary = f"results/07_summary/{config['ref_name']}_{config['query_name']}/conservation_by_tolerance.tsv",
        target_sites = f"results/07_summary/{config['ref_name']}_{config['query_name']}/target_sites.tsv.gz",
        relationships = f"results/07_summary/{config['ref_name']}_{config['query_name']}/site_relationships.tsv.gz",
        guide_summary = f"results/07_summary/{config['ref_name']}_{config['query_name']}/guide_targeting_summary.tsv.gz"
    params:
        out_dir = f"results/07_summary/{config['ref_name']}_{config['query_name']}",
        tol     = config.get("cfd", {}).get("tolerance", 1000),
        ref_name = config["ref_name"],
        query_name = config["query_name"]
    log:
        f"logs/07_summary/map_pam_orthology_{config['ref_name']}_{config['query_name']}.log"
    container:
        "workflow/containers/images/python_analysis.sif"
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
            --source-manifest {input.source_manifest} \
            --ref-name {params.ref_name} \
            --query-name {params.query_name} \
            --out-dir {params.out_dir} \
            --tol {params.tol} > {log} 2>&1
        """


rule summarize_cfd:
    input:
        conserved = rules.map_pam_orthology.output.conserved,
        failed = rules.map_pam_orthology.output.failed
    output:
        report = f"results/07_summary/{config['ref_name']}_{config['query_name']}/cfd_summary.txt"
    log:
        f"logs/07_summary/summarize_cfd_{config['ref_name']}_{config['query_name']}.log"
    container:
        "workflow/containers/images/python_analysis.sif"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 60
    shell:
        """
        python3 workflow/scripts/crispr_cfd_summary.py \
            {input.conserved} \
            {input.failed} \
            > {output.report} 2> {log}
        """


rule summarize_targeting_landscape:
    input:
        # Use the mapper's declared output, not a hard-coded assembly pair.
        guide_summary = rules.map_pam_orthology.output.guide_summary
    output:
        report = f"results/07_summary/{config['ref_name']}_{config['query_name']}/targeting_landscape_summary.txt",
        table = f"results/07_summary/{config['ref_name']}_{config['query_name']}/targeting_landscape_summary.tsv"
    log:
        f"logs/07_summary/summarize_targeting_landscape_{config['ref_name']}_{config['query_name']}.log"
    params:
        ref_label = config["ref_name"],
        query_label = config["query_name"]
    container:
        "workflow/containers/images/python_analysis.sif"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 60
    shell:
        """
        python3 workflow/scripts/summarize_targeting_landscape.py \
            --guide-summary {input.guide_summary:q} \
            --report {output.report:q} \
            --table {output.table:q} \
            --ref-label {params.ref_label:q} \
            --query-label {params.query_label:q} \
            > {log:q} 2>&1
        """
