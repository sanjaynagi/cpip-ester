RESULTS_DIR = "results"

# Calculate GC content for genomic windows
rule calculate_gc_content:
    input:
        reference = REF_GENOME
    output:
        gc_content = "results/gc_content.tsv"
    log:
        "logs/cnv/gc.log"
    params:
        window_size = WINDOW_SIZE,
        interval = INTERVAL
    threads: 10
    conda:
        "../scripts/env.yaml"  # Optional: create conda environment with dependencies
    script:
        "../scripts/gc_content.py"

# Extract read counts in windows
rule extract_counts_for_hmm:
    input:
        bam = "results/region_bams/{sample}.sorted.bam"
    output:
        counts = f"{RESULTS_DIR}/cnv/counts/{{sample}}/counts_for_hmm.tsv"
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
        bam = "results/region_bams/{sample}.sorted.bam"
    output:
        mapq0_proportions = f"{RESULTS_DIR}/cnv/counts/{{sample}}/mapq0_proportions.tsv"
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
        mapq0_proportions = f"{RESULTS_DIR}/cnv/counts/{{sample}}/mapq0_proportions.tsv"
    output:
        accessibility = f"{RESULTS_DIR}/cnv/counts/{{sample}}/accessibility.tsv"
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
        counts = f"{RESULTS_DIR}/cnv/counts/{{sample}}/counts_for_hmm.tsv",
        gc_content = "results/gc_content.tsv",
        accessibility = f"{RESULTS_DIR}/cnv/counts/{{sample}}/accessibility.tsv"
    output:
        median_coverage = f"{RESULTS_DIR}/cnv/counts/{{sample}}/median_coverage_by_gc.tsv",
        variance = f"{RESULTS_DIR}/cnv/counts/{{sample}}/coverage_variance.tsv"
    conda:
        "../scripts/env.yaml"
    log:
        "logs/cnv/covgc/{sample}.log"
    script:
        "../scripts/calculate_median_coverage_by_gc.py"

# Apply HMM for CNV calling
rule run_hmm:
    input:
        counts = f"{RESULTS_DIR}/cnv/counts/{{sample}}/counts_for_hmm.tsv",
        gc_content = "results/gc_content.tsv",
        median_coverage = f"{RESULTS_DIR}/cnv/counts/{{sample}}/median_coverage_by_gc.tsv",
        variance = f"{RESULTS_DIR}/cnv/counts/{{sample}}/coverage_variance.tsv",
        accessibility = f"{RESULTS_DIR}/cnv/counts/{{sample}}/accessibility.tsv",
        mapq0_proportions = f"{RESULTS_DIR}/cnv/counts/{{sample}}/mapq0_proportions.tsv"
    output:
        hmm_output = f"{RESULTS_DIR}/cnv/hmm/{{sample}}/hmm_output.tsv",
        cnv_calls = f"{RESULTS_DIR}/cnv/hmm/{{sample}}/cnv_calls.csv"
    log:
        "logs/cnv/hmm/{sample}.log"
    conda:
        "../scripts/env.yaml"
    params:
        transition_probability = TRANSITION_PROBABILITY,
        max_copy_number = MAX_COPY_NUMBER,
        max_mapq0_proportion = MAX_MAPQ0_PROPORTION,
        window_size = WINDOW_SIZE
    script:
        "../scripts/hmm_process.py"

# Detect breakpoints for CNVs
rule detect_breakpoints:
    input:
        bam = f"{RESULTS_DIR}/region_bams/{{sample}}.sorted.bam",
        cnv_calls = f"{RESULTS_DIR}/cnv/hmm/{{sample}}/cnv_calls.csv"
    output:
        breakpoints = f"{RESULTS_DIR}/cnv/breakpoints/{{sample}}/breakpoints.csv",
        pre_clipping_fastq = f"{RESULTS_DIR}/cnv/breakpoints/{{sample}}/pre_clipping.fastq",
        post_clipping_fastq = f"{RESULTS_DIR}/cnv/breakpoints/{{sample}}/post_clipping.fastq"
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
        cnv_files = expand(f"{RESULTS_DIR}/cnv/hmm/{{sample}}/cnv_calls.csv", sample=SAMPLE_NAMES),
        annotation = GFF
    output:
        modal_cnvs = f"{RESULTS_DIR}/modal_cnvs.csv"
    log:
        "logs/cnv/modal_cnv.log"
    conda:
        "../scripts/env.yaml"
    params:
        samples = SAMPLE_NAMES
    script:
        "../scripts/modal_cnv.py"