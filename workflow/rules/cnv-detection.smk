RESULTS_DIR = "results"

# Calculate GC content for genomic windows
rule calculate_gc_content:
    input:
        reference = CONTIG_REF_GENOME
    output:
        gc_content = "results/gc_content.tsv"
    log:
        "logs/cnv/gc.log"
    params:
        window_size = WINDOW_SIZE,
        interval = INTERVAL
    threads: 10
    conda:
        "../scripts/env.yaml" 
    script:
        "../scripts/gc_content.py"

# Extract read counts in windows
rule extract_counts_for_hmm:
    input:
        bam = "results/region_bams/{sample}.sorted.bam",
        bai = "results/region_bams/{sample}.sorted.bam.bai"
    output:
        counts = "results/cnv/{sample}/counts_for_hmm.tsv"
    log:
        "logs/cnv/counts/{sample}.log"
    conda:
        "../scripts/env.yaml"
    params:
        window_size = WINDOW_SIZE,
        interval = INTERVAL,
        min_mapq = MAPQ_THRESHOLD
    script:
        "../scripts/counts_for_hmm.py"

# Calculate the proportion of MAPQ0 reads in each window
rule calculate_mapq0_proportions:
    input:
        bam = "results/region_bams/{sample}.sorted.bam",
        bai = "results/region_bams/{sample}.sorted.bam.bai"
    output:
        mapq0_proportions = "results/cnv/{sample}/mapq0_proportions.tsv"
    log:
        "logs/cnv/mapq0/{sample}.log"
    conda:
        "../scripts/env.yaml"
    params:
        window_size = WINDOW_SIZE,
        interval = INTERVAL
    script:
        "../scripts/calculate_mapq0_proportions.py"

# Calculate windowed accessibility from GC content and MAPQ0 proportions
rule calculate_accessibility:
    input:
        gc_content = "results/gc_content.tsv",
        mapq0_proportions = "results/cnv/{sample}/mapq0_proportions.tsv"
    output:
        accessibility = "results/cnv/{sample}/accessibility.tsv"
    log:
        "logs/cnv/accessibility/{sample}.log"
    conda:
        "../scripts/env.yaml"
    params:
        mapq0_threshold = MAX_MAPQ0_PROPORTION,
        accessibility_threshold = ACCESSIBILITY_THRESHOLD
    script:
        "../scripts/calculate_accessibility.py"

# Calculate median coverage by GC content
rule calculate_median_coverage_by_gc:
    input:
        counts = expand("results/cnv/{sample}/counts_for_hmm.tsv", sample=SAMPLE_NAMES),
        gc_content = "results/gc_content.tsv",
        accessibility = expand("results/cnv/{sample}/accessibility.tsv", sample=SAMPLE_NAMES),
        mapq0= expand("results/cnv/{sample}/mapq0_proportions.tsv", sample=SAMPLE_NAMES)
    output:
        median_coverage = "results/cnv/median_coverage_by_gc.tsv",
        variance = "results/cnv/coverage_variance.tsv"
    params:
        accessibility_threshold = ACCESSIBILITY_THRESHOLD,
        mapq0_threshold=MAPQ_THRESHOLD
    conda:
        "../scripts/env.yaml"
    log:
        "logs/cnv/covgc.log"
    script:
        "../scripts/calculate_median_coverage_by_gc.py"

rule run_hmm:
    input:
        gc_content = "results/gc_content.tsv",
        median_coverage = "results/cnv/median_coverage_by_gc.tsv",
        variance = "results/cnv/coverage_variance.tsv",
        accessibility = expand("results/cnv/{sample}/accessibility.tsv", sample=SAMPLE_NAMES),
        mapq0_proportions = expand("results/cnv/{sample}/mapq0_proportions.tsv", sample=SAMPLE_NAMES)
    output:
        hmm_output = expand("results/cnv/{sample}/hmm_output.tsv", sample=SAMPLE_NAMES),
        cnv_calls = expand("results/cnv/{sample}/cnv_calls.csv", sample=SAMPLE_NAMES)
    log:
        "logs/cnv/hmm.log"
    conda:
        "../scripts/env.yaml"
    params:
        counts_dir = "results/cnv",
        hmm_output_dir = "results/cnv",
        transition_probability = TRANSITION_PROBABILITY,
        max_copy_number = MAX_COPY_NUMBER,
        max_mapq0_proportion = MAX_MAPQ0_PROPORTION,
        window_size = WINDOW_SIZE
    script:
        "../scripts/hmm_process.py"


# Detect breakpoints for CNVs
rule detect_breakpoints:
    input:
        bam = "results/region_bams/{sample}.sorted.bam",
        cnv_calls = "results/cnv/{sample}/cnv_calls.csv"
    output:
        breakpoints = "results/cnv/{sample}/breakpoints.csv",
        pre_clipping_fastq = "results/cnv/{sample}/pre_clipping.fastq",
        post_clipping_fastq = "results/cnv/{sample}/post_clipping.fastq"
    log:
        "logs/cnv/breakpoints/{sample}.log"
    conda:
        "../scripts/env.yaml"
    params:
        min_mapq = MAPQ_THRESHOLD
    script:
        "../scripts/breakpoint_detector.py"

# Calculate modal CNVs for all genes in GFF
rule calculate_modal_cnvs:
    input:
        cnv_files = expand("results/cnv/{sample}/cnv_calls.csv", sample=SAMPLE_NAMES),
        annotation = GFF
    output:
        modal_cnvs = "results/modal_cnvs.csv"
    log:
        "logs/cnv/modal_cnv.log"
    conda:
        "../scripts/env.yaml"
    params:
        samples = SAMPLE_NAMES
    script:
        "../scripts/modal_cnv.py"