# workflow/rules/08_figures/figures.smk

rule plot_synteny_diagnostics:
    input:
        conserved = rules.map_pam_orthology.output.conserved,
        failed = rules.map_pam_orthology.output.failed,
        unmapped_ref = rules.map_pam_orthology.output.unmapped_ref,
        unmapped_qry = rules.map_pam_orthology.output.unmapped_qry,
        ref_cfd = rules.compute_cfd_iso1.output.cfd,
        qry_cfd = rules.compute_cfd_bl.output.cfd
    output:
        figures = directory("results/08_figures"),
        plot = "results/08_figures/cfd_shift.png"
    log:
        "logs/08_figures/plot_synteny_diagnostics.log"
    params:
        data_dir = "results/07_summary"
    container:
        "containers/python_analysis.sif"
    shell:
        """
        python3 workflow/scripts/synteny_plotter.py \
            --data-dir {params.data_dir} \
            --figures {output.figures} \
            --ref-cfd {input.ref_cfd} \
            --query-cfd {input.qry_cfd} > {log} 2>&1
        """

rule plot_tolerance_venns:
    input:
        ref_cfd = rules.compute_cfd_iso1.output.cfd,
        qry_cfd = rules.compute_cfd_bl.output.cfd,
        delta = rules.nucmer_synteny.output.delta
    output:
        sunburst_plot = "results/08_figures/multi_tolerance_sunburst.png"
    log:
        "logs/08_figures/plot_tolerance_venns.log"
    params:
        figures_dir = "results/08_figures"
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
