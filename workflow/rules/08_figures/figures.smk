# workflow/rules/08_figures/figures.smk

rule plot_synteny_diagnostics:
    input:
        conserved    = rules.map_pam_orthology.output.conserved,
        failed       = rules.map_pam_orthology.output.failed,
        unmapped_ref = rules.map_pam_orthology.output.unmapped_ref,
        unmapped_qry = rules.map_pam_orthology.output.unmapped_qry,
        ref_cfd      = rules.compute_cfd_ref.output.cfd,
        qry_cfd      = rules.compute_cfd_qry.output.cfd
    output:
        plot    = f"results/08_figures/{config['ref_name']}_{config['query_name']}/cfd_shift.png"
    log:
        f"logs/08_figures/plot_synteny_diagnostics_{config['ref_name']}_{config['query_name']}.log"
    params:
        data_dir = rules.map_pam_orthology.params.out_dir,
        figures = f"results/08_figures/{config['ref_name']}_{config['query_name']}"
    container:
        "containers/images/python_analysis.sif"
    shell:
        """
        python3 workflow/scripts/synteny_plotter.py \
            --data-dir {params.data_dir} \
            --figures {params.figures} \
            --ref-cfd {input.ref_cfd} \
            --query-cfd {input.qry_cfd} > {log} 2>&1
        """

rule plot_tolerance_venns:
    input:
        ref_cfd = rules.compute_cfd_ref.output.cfd,
        qry_cfd = rules.compute_cfd_qry.output.cfd,
        delta   = rules.nucmer_synteny.output.delta
    output:
        sunburst_plot = f"results/08_figures/{config['ref_name']}_{config['query_name']}/multi_tolerance_sunburst.png"
    log:
        f"logs/08_figures/plot_tolerance_venns_{config['ref_name']}_{config['query_name']}.log"
    params:
        figures_dir = f"results/08_figures/{config['ref_name']}_{config['query_name']}"
    container:
        "containers/python_analysis.sif"
    shell:
        """
        python3 workflow/scripts/venn_diagrams.py \
            --ref-cfd {input.ref_cfd} \
            --query-cfd {input.qry_cfd} \
            --delta {input.delta} \
            --figures {params.figures_dir} > {log} 2>&1
        """