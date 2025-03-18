RESULTS_DIR = "results"

# Calculate GC content for genomic windows
rule calculate_gc_content:
    input:
        reference = REF_GENOME
    output:
        gc_content = "results/gc_content.tsv"
    log:
        "logs/cnv/gc.log"
    shell:
        """
        bedtools makewindows -g <(samtools faidx {input.reference}) -w {WINDOW_SIZE} -s {INTERVAL} | \
        bedtools nuc -fi {input.reference} -bed - | \
        awk 'BEGIN{{OFS="\\t"}} NR>1 {{print $1,$2,$5}}' > {output.gc_content}
        """

# Extract read counts in windows
rule extract_counts_for_hmm:
    input:
        bam = "results/region_bams/{sample}.sorted.bam"
    output:
        counts = f"{RESULTS_DIR}/{{sample}}/counts/counts_for_hmm.tsv"
    log:
        "logs/cnv/counts/{sample}.log"
    params:
        window_size = WINDOW_SIZE,
        interval = INTERVAL,
        min_mapq = MAPQ_THRESHOLD
    script:
        "workflow/scripts/counts_for_hmm.py"

# Calculate the proportion of MAPQ0 reads in each window
rule calculate_mapq0_proportions:
    input:
        bam = "results/region_bams/{sample}.sorted.bam"
    output:
        mapq0_proportions = f"{RESULTS_DIR}/{{sample}}/counts/mapq0_proportions.tsv"
    log:
        "logs/cnv/mapq0/{sample}.log"
    params:
        window_size = WINDOW_SIZE,
        interval = INTERVAL
    script:
        "workflow/scripts/calculate_mapq0_proportions.py"

# Calculate windowed accessibility from GC content and MAPQ0 proportions
rule calculate_accessibility:
    input:
        gc_content = "results/gc_content.tsv",
        mapq0_proportions = f"{RESULTS_DIR}/{{sample}}/counts/mapq0_proportions.tsv"
    output:
        accessibility = f"{RESULTS_DIR}/{{sample}}/counts/accessibility.tsv"
    log:
        "logs/cnv/accessibility/{sample}.log"
    params:
        mapq0_threshold = MAX_MAPQ0_PROPORTION,
        accessibility_threshold = ACCESSIBILITY_THRESHOLD
    script:
        "workflow/scripts/calculate_accessibility.py"

# Calculate median coverage by GC content
rule calculate_median_coverage_by_gc:
    input:
        counts = f"{RESULTS_DIR}/{{sample}}/counts/counts_for_hmm.tsv",
        gc_content = "results/gc_content.tsv",
        accessibility = f"{RESULTS_DIR}/{{sample}}/counts/accessibility.tsv"
    output:
        median_coverage = f"{RESULTS_DIR}/{{sample}}/counts/median_coverage_by_gc.tsv",
        variance = f"{RESULTS_DIR}/{{sample}}/counts/coverage_variance.tsv"
    log:
        "logs/cnv/covgc/{sample}.log"
    script:
        "workflow/scripts/calculate_median_coverage_by_gc.py"

# Apply HMM for CNV calling
rule run_hmm:
    input:
        counts = f"{RESULTS_DIR}/{{sample}}/counts/counts_for_hmm.tsv",
        gc_content = "results/gc_content.tsv",
        median_coverage = f"{RESULTS_DIR}/{{sample}}/counts/median_coverage_by_gc.tsv",
        variance = f"{RESULTS_DIR}/{{sample}}/counts/coverage_variance.tsv",
        accessibility = f"{RESULTS_DIR}/{{sample}}/counts/accessibility.tsv",
        mapq0_proportions = f"{RESULTS_DIR}/{{sample}}/counts/mapq0_proportions.tsv"
    output:
        hmm_output = f"{RESULTS_DIR}/{{sample}}/hmm/hmm_output.tsv",
        cnv_calls = f"{RESULTS_DIR}/{{sample}}/hmm/cnv_calls.csv"
    log:
        "logs/cnv/hmm/{sample}.log"
    params:
        transition_probability = TRANSITION_PROBABILITY,
        max_copy_number = MAX_COPY_NUMBER,
        max_mapq0_proportion = MAX_MAPQ0_PROPORTION
    script:
        "workflow/scripts/hmm_process.py"

# Detect breakpoints for CNVs
rule detect_breakpoints:
    input:
        bam = f"{RESULTS_DIR}/region_bams/{{sample}}.sorted.bam",
        cnv_calls = f"{RESULTS_DIR}/{{sample}}/hmm/cnv_calls.csv"
    output:
        breakpoints = f"{RESULTS_DIR}/{{sample}}/breakpoints/breakpoints.csv",
        pre_clipping_fastq = f"{RESULTS_DIR}/{{sample}}/breakpoints/pre_clipping.fastq",
        post_clipping_fastq = f"{RESULTS_DIR}/{{sample}}/breakpoints/post_clipping.fastq"
    log:
        "logs/cnv/breakpoints/{sample}.log"
    params:
        min_mapq = MAPQ_THRESHOLD
    script:
        "workflow/scripts/breakpoint_detector.py"

# Calculate modal CNVs for all genes in GFF
rule calculate_modal_cnvs:
    input:
        cnv_files = expand(f"{RESULTS_DIR}/{{sample}}/hmm/cnv_calls.csv", sample=SAMPLE_NAMES),
        annotation = GFF
    output:
        modal_cnvs = f"{RESULTS_DIR}/modal_cnvs.csv"
    log:
        "logs/cnv/modal_cnv.log"
    params:
        samples = SAMPLE_NAMES
    script:
        "workflow/scripts/modal_cnv.py"