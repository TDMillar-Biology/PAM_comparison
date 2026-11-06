# workflow/rules/01_scaffolding/scaffold.smk

rule mummer_for_scaffolding:
    input:
        reference = config["ref_genome"],
        query = config["query_assembly"]
    output:
        delta = f"results/01_scaffolding/{config['ref_name']}_{config['query_name']}.delta"
    params:
        prefix = f"results/01_scaffolding/{config['ref_name']}_{config['query_name']}"
    resources:
        mem_mb = 16000,
        runtime = 120,
        tasks = 1
    threads: config.get("threads", 8)
    log:
        "logs/01_scaffolding/mummer_for_scaffolding.log"
    container:
        "containers/01_scaffolding.sif"
    shell:
        "nucmer -p {params.prefix} {input.reference} {input.query} -t {threads} > {log} 2>&1"

rule daedalus_scaffold_assembly:
    input:
        delta = rules.mummer_for_scaffolding.output.delta,
        qry = config["query_assembly"]
    output:
        coords = f"results/01_scaffolding/{config['ref_name']}_{config['query_name']}.coords",
        scaffold = f"results/01_scaffolding/{config['query_name']}.scaffolded.fasta",
        agp = f"results/01_scaffolding/{config['query_name']}.scaffolded.agp"
    params:
        prefix = f"{config['query_name']}.scaffolded",
        chromosomes = " ".join(config["main_chromosomes"])
    log:
        "logs/01_scaffolding/daedalus_scaffold_assembly.log"
    container:
        "containers/01_scaffolding.sif"
    threads: 4
    resources:
        mem_mb = 16000,
        runtime = 120
    shell:
        """
        show-coords -rclT {input.delta} > {output.coords} 2> {log}

        daedalus scaffold \
            --coords {output.coords} \
            --fasta {input.qry} \
            --prefix {params.prefix} \
            --out_dir results/01_scaffolding \
            --chromosomes {params.chromosomes} \
            >> {log} 2>&1
        """

rule extract_main_scaffolds:
    input:
        scaffold = rules.daedalus_scaffold_assembly.output.scaffold,
        agp = rules.daedalus_scaffold_assembly.output.agp
    output:
        main_fasta = f"results/01_scaffolding/{config['query_name']}.canonical.scaffolded.fasta",
        scaffold_list = f"results/01_scaffolding/{config['query_name']}.main_scaffolds.txt"
    params:
        chrom_patterns = "|".join([str(chrom) for chrom in config.get("main_chromosomes", ["X", "2L", "2R", "3L", "3R", "4", "Y", "mitochondrion_genome"])])
    log:
        "logs/01_scaffolding/extract_main_scaffolds.log"
    container:
        "containers/01_scaffolding.sif"
    shell:
        """
        awk '{{print $1}}' {input.agp} | sort -u | grep -E "^({params.chrom_patterns})" > {output.scaffold_list} 2> {log} || true

        if [ ! -s {output.scaffold_list} ]; then
            awk '{{print $1}}' {input.agp} | sort -u > {output.scaffold_list} 2>> {log}
        fi

        seqtk subseq {input.scaffold} {output.scaffold_list} > {output.main_fasta} 2>> {log}
        """