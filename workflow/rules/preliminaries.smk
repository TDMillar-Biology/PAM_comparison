rule scaffold_assembly:
    input:
        ref = config["ref_genome"],
        qry = config["query_assembly"]
    output:
        scaffold = "results/01_scaffolding/ragtag.scaffold.fasta",
        agp = "results/01_scaffolding/ragtag.scaffold.agp"
    log:
        "logs/01_scaffolding/ragtag.log"
    conda: 
        "../envs/scaffolding.yaml"
    shell:
        "ragtag.py scaffold {input.ref} {input.qry} -o results/01_scaffolding/ 2> {log}"

rule extract_main_scaffolds:
    input:
        scaffold = "results/01_scaffolding/ragtag.scaffold.fasta",
        agp = "results/01_scaffolding/ragtag.scaffold.agp"
    output:
        main_fasta = "results/01_scaffolding/BL54591.canonical.scaffolded.fasta"
    log:
        "logs/01_scaffolding/extract_main.log"
    conda: 
        "../envs/utilities.yaml"
    shell:
        """
        awk '$1 ~ /_RagTag$/ {{print $1}}' {input.agp} | sort -u > results/01_scaffolding/main_scaffolds.txt 
        
        if [ ! -s results/01_scaffolding/main_scaffolds.txt ]; then
            echo "Error: No canonical scaffolds found in AGP." > {log}
            exit 1
        fi
        
        seqtk subseq {input.scaffold} results/01_scaffolding/main_scaffolds.txt > {output.main_fasta} 2>> {log}
        """

rule extract_and_fix_euchromatin:
    input:
        ref = config["ref_genome"],
        bed = config["euchromatin_bed"]
    output:
        fixed_fa = "results/02_preprocessing/ISO1-r6.58_euchromatin.fixed.fasta"
    conda: 
        "../envs/utilities.yaml"
    shell:
        """
        bedtools getfasta -fi {input.ref} -bed {input.bed} | \
        sed -e 's/X:277911-18930000/X/g' \
            -e 's/2L:82455-19570000/2L/g' \
            -e 's/2R:8860000-24684540/2R/g' \
            -e 's/3L:158639-18438500/3L/g' \
            -e 's/3R:9497000-31845060/3R/g' > {output.fixed_fa}
        """

rule repeat_masker:
    input:
        fa = "results/02_preprocessing/ISO1-r6.58_euchromatin.fixed.fasta",
        lib = config["repeat_lib"]
    output:
        masked = "results/02_preprocessing/ISO1-r6.58_euchromatin.fixed.fasta.masked"
    log:
        "logs/02_preprocessing/repeatmasker.log"
    conda: 
        "../envs/masking.yaml"
    threads: config["threads"]
    shell:
        "RepeatMasker -lib {input.lib} -pa {threads} -dir results/02_preprocessing/ {input.fa} > {log} 2>&1"

rule discover_pams:
    input:
        masked = "results/02_preprocessing/ISO1-r6.58_euchromatin.fixed.fasta.masked"
    output:
        fa = "results/03_pam_discovery/ISO1.masked.euchromatic.PAMs.fa",
        bed = "results/03_pam_discovery/ISO1.masked.euchromatic.PAMs.bed"
    conda: 
        "../envs/python_bio.yaml"
    shell:
        "python3 workflow/scripts/pam_discovery.py -i {input.masked} -o results/03_pam_discovery/ISO1.masked.euchromatic.PAMs"

rule sample_pams:
    input:
        all_pams = "results/03_pam_discovery/ISO1.masked.euchromatic.PAMs.fa"
    output:
        sample = "results/03_pam_discovery/ISO1.masked.euchromatic.PAMs.10pct.fa"
    params:
        percentage = "0.01",
        seed = "123"
    conda: 
        "../envs/utilities.yaml"
    threads: 1 
    resources:
        mem_mb = 8000, 
        runtime = 15
    shell:
        "seqtk sample -s {params.seed} {input.all_pams} {params.percentage} > {output.sample}"

rule nucmer_synteny:
    input:
        ref = "results/02_preprocessing/ISO1-r6.58_euchromatin.fixed.fasta",
        qry = "results/01_scaffolding/BL54591.canonical.scaffolded.fasta"
    output:
        delta = "results/04_synteny/ISO1_BL54591.delta"
    log:
        "logs/04_synteny/nucmer.log"
    conda: 
        "../envs/alignment.yaml"
    threads: 1  
    resources:
        mem_mb = 16000,   
        runtime = 120    
    shell:
        "nucmer -p results/04_synteny/ISO1_BL54591 {input.ref} {input.qry} 2> {log}"