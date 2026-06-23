
rule scaffold_assembly:
    input:
        ref = config["ref_genome"],
        qry = config["query_assembly"]
    output:
        scaffold = "results/01_scaffolding/ragtag.scaffold.fasta",
        agp = "results/01_scaffolding/ragtag.scaffold.agp"
    conda: 
        "../envs/scaffolding.yaml"
    shell:
        """
        ragtag.py scaffold {input.ref} {input.qry} -o results/01_scaffolding/
        """

rule extract_main_scaffolds:
    input:
        scaffold = "results/01_scaffolding/ragtag.scaffold.fasta",
        agp = "results/01_scaffolding/ragtag.scaffold.agp"
    output:
        main_fasta = "results/01_scaffolding/BL54591.canonical.scaffolded.fasta"
    conda: 
        "../envs/utilities.yaml"
    shell:
        """
        ## grab names of scaffolded components
        awk '$1 ~ /_RagTag$/ {{print $1}}' {input.agp} | sort -u > results/01_scaffolding/main_scaffolds.txt 
        
        if [ ! -s results/01_scaffolding/main_scaffolds.txt ]; then
            echo "Error: No canonical scaffolds found in AGP."
            exit 1
        fi
        
        seqtk subseq {input.scaffold} results/01_scaffolding/main_scaffolds.txt > {output.main_fasta}
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
    conda: 
        "../envs/masking.yaml"
    threads: config["threads"]
    shell:
        """
        RepeatMasker -lib {input.lib} -pa {threads} -dir results/02_preprocessing/ {input.fa}
        """

rule discover_pams:
    input:
        masked = "results/02_preprocessing/ISO1-r6.58_euchromatin.fixed.fasta.masked"
    output:
        fa = "results/03_pam_discovery/ISO1.masked.euchromatic.PAMs.fa",
        bed = "results/03_pam_discovery/ISO1.masked.euchromatic.PAMs.bed"
    conda: 
        "../envs/python_bio.yaml"
    shell:
        """
        python3 workflow/scripts/pam_discovery.py -i {input.masked} -o results/03_pam_discovery/ISO1.masked.euchromatic.PAMs
        """

rule sample_pams:
    input:
        all_pams = "results/03_pam_discovery/ISO1.masked.euchromatic.PAMs.fa"
    output:
        sample = "results/03_pam_discovery/ISO1.masked.euchromatic.PAMs.10pct.fa"
    params:
        percentage = "0.1",
        seed = "123"
    conda: 
        "../envs/utilities.yaml"
    threads: 1 
    resources:
        mem_mb = 8000, 
        runtime = 15
    shell:
        """
        seqtk sample -s {params.seed} {input.all_pams} {params.percentage} > {output.sample}
        """

rule nucmer_synteny:
    input:
        ref = "results/02_preprocessing/ISO1-r6.58_euchromatin.fixed.fasta",
        qry = "results/01_scaffolding/BL54591.canonical.scaffolded.fasta"
    output:
        delta = "results/04_synteny/ISO1_BL54591.delta"
    conda: 
        "../envs/alignment.yaml"
    threads: 1  
    resources:
        mem_mb = 16000,   
        runtime = 120    
    shell:
        """
        nucmer -p results/04_synteny/ISO1_BL54591 {input.ref} {input.qry}
        """

rule bowtie_index_ref:
    input:
        ref = config["ref_genome"]
    output:
        idx = "results/05_alignment/ISO1/ISO1.1.ebwt" 
    params:
        prefix = "results/05_alignment/ISO1/ISO1"
    conda: 
        "../envs/alignment.yaml"
    shell:
        "bowtie-build {input.ref} {params.prefix}"

rule bowtie_align_ref:
    input:
        idx = "results/05_alignment/ISO1/ISO1.1.ebwt",
        pams = get_pam_fasta
    output:
        sam = "results/05_alignment/ISO1_pams.sam"
    params:
        prefix = "results/05_alignment/ISO1/ISO1",
        mismatch = config["bowtie_mismatches"]
    conda: 
        "../envs/alignment.yaml"
    threads: config["threads"]
    resources:
        mem_mb = 32000,   
        runtime = 360    
    shell:
        """
        bowtie -f -v {params.mismatch} -a --best --sam -p {threads} {params.prefix} {input.pams} > {output.sam}
        """

rule bowtie_align_sample:
    input:
        idx = "results/05_alignment/ISO1/ISO1.1.ebwt",
        pams = get_pam_fasta
    output:
        sam = "results/05_alignment/ISO1_pams_sample.sam"
    params:
        prefix = "results/05_alignment/ISO1/ISO1",
        mismatch = config["bowtie_mismatches"]
    conda: 
        "../envs/alignment.yaml"
    threads: config["threads"]
    resources:
        mem_mb = 32000,   
        runtime = 360    
    shell:
        """
        bowtie -f -v {params.mismatch} -a --best --sam -p {threads} {params.prefix} {input.pams} > {output.sam}
        """


# ==========================================
# Missing Stage 2 Rules: BL54591 Alignment
# ==========================================
rule bowtie_index_qry:
    input:
        ref = "results/01_scaffolding/BL54591.canonical.scaffolded.fasta"
    output:
        idx = "results/05_alignment/BL54591/BL54591.1.ebwt" 
    params:
        prefix = "results/05_alignment/BL54591/BL54591"
    conda: 
        "../envs/alignment.yaml"
    shell:
        "bowtie-build {input.ref} {params.prefix}"

rule bowtie_align_qry:
    input:
        idx = "results/05_alignment/BL54591/BL54591.1.ebwt",
        pams = get_pam_fasta # Using the same input function to respect the 10% test mode!
    output:
        sam = "results/05_alignment/BL54591_pams.sam"
    params:
        prefix = "results/05_alignment/BL54591/BL54591",
        mismatch = config["bowtie_mismatches"]
    conda: 
        "../envs/alignment.yaml"
    threads: config["threads"]
    resources:
        mem_mb = 32000,
        runtime = 360
    shell:
        """
        bowtie -f -v {params.mismatch} -a --best --sam -p {threads} {params.prefix} {input.pams} > {output.sam}
        """

# ==========================================
# Stage 3: Compute CFD Scores
# ==========================================
rule compute_cfd_iso1:
    input:
        sam = "results/05_alignment/ISO1_pams.sam",
        ref = config["ref_genome"],
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

rule classify_pam_orthology:
    input:
        ref_cfd = "results/06_cfd_scores/iso1_raw_cfd.csv",
        query_cfd = "results/06_cfd_scores/bl_raw_cfd.csv",
        delta = "results/04_synteny/ISO1_BL54591.delta"
    output:
        summary = "results/07_summary/pam_orthology_summary.csv"
    conda:
        "../envs/python_bio.yaml"
    resources:
        mem_mb = 8000,
        runtime = 60
    shell:
        """
        python3 workflow/scripts/classify_pam_orthology.py \
            --ref-cfd {input.ref_cfd} \
            --query-cfd {input.query_cfd} \
            --delta {input.delta} \
            --out {output.summary} \
            --tol 1000
        """

# ==========================================
# OLD: Annotate Synteny Blocks -- set for deprecation
# ==========================================
rule annotate_synteny_iso1:
    input:
        cfd = "results/06_cfd_scores/iso1_raw_cfd.csv",
        delta = "results/04_synteny/ISO1_BL54591.delta"
    output:
        synt = "results/06_cfd_scores/iso1_synteny.csv"
    conda:
        "../envs/python_bio.yaml"
    resources:
        mem_mb = 8000,
        runtime = 30
    shell:
        """
        python3 workflow/scripts/annotate_synteny.py --cfd {input.cfd} --delta {input.delta} --out {output.synt}
        """

rule annotate_synteny_bl:
    input:
        cfd = "results/06_cfd_scores/bl_raw_cfd.csv",
        delta = "results/04_synteny/ISO1_BL54591.delta"
    output:
        synt = "results/06_cfd_scores/bl_synteny.csv"
    conda:
        "../envs/python_bio.yaml"
    resources:
        mem_mb = 8000,
        runtime = 30
    shell:
        """
        python3 workflow/scripts/annotate_synteny.py --cfd {input.cfd} --delta {input.delta} --out {output.synt}
        """

# ==========================================
# OLD: Merge Synteny-Corrected Tables -- set for deprecation
# ==========================================
rule merge_synteny_tables:
    input:
        iso1_synt = "results/06_cfd_scores/iso1_synteny.csv",
        bl_synt = "results/06_cfd_scores/bl_synteny.csv"
    output:
        final_csv = "results/07_summary/pam_summary.csv"
    conda:
        "../envs/python_bio.yaml"
    resources:
        mem_mb = 4000,
        runtime = 15
    shell:
        """
        python3 workflow/scripts/merge_synteny_corrected.py {input.iso1_synt} {input.bl_synt} {output.final_csv}
        """
