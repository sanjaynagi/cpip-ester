

# Align reads to whole chromosome reference
rule align_reads:
    input:
        r1="results/temp_reads/{sample}_1.fastq.gz",
        r2="results/temp_reads/{sample}_2.fastq.gz",
        ref=REF_GENOME
    output:
        pipe("results/pipes/{sample}.sam")
    log:
        "logs/align/bwa/{sample}.log"
    params:
        rg=lambda wildcards: f"@RG\\tID:{sample_to_srr[wildcards.sample]}\\tSM:{wildcards.sample}"
    threads: 8
    shell:
        "resources/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t {threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2} > {output} 2> {log}"

# Convert SAM to sorted BAM for the whole chromosome, filtering out unmapped reads
rule sam_to_sorted_bam:
    input:
        sam="results/pipes/{sample}.sam"
    output:
        bam="results/region_bams/{sample}.sorted.bam"
    threads: 4  # Increased threads since we're doing more work
    log:
        "logs/align/sam_to_bam/{sample}.log"
    shell:
        """
        samtools view -b -h -F 4 -@ {threads} {input.sam} | \
        samtools sort -@ {threads} -o {output.bam} - 2> {log}
        """

# Index the region BAM files
rule index_region_bam:
    input:
        bam="results/region_bams/{sample}.sorted.bam"
    output:
        bai="results/region_bams/{sample}.sorted.bam.bai"
    log:
        "logs/align/bam_index_region/{sample}.log"
    shell:
        "samtools index {input.bam} 2> {log}"


# # Index sorted BAM
# rule index_temp_bam:
#     input:
#         bam="results/temp_bams/{sample}.sorted.bam"
#     output:
#         bai="results/temp_bams/{sample}.sorted.bam.bai"
#     log:
#         "logs/bam_index_temp/{sample}.log"
#     shell:
#         "samtools index {input.bam} 2> {log}"

# Calculate coverage in 300bp windows using mosdepth
rule calculate_coverage:
    input:
        bam="results/region_bams/{sample}.sorted.bam",
        bai="results/region_bams/{sample}.sorted.bam.bai"
    output:
        regions="results/coverage/{sample}.regions.bed.gz",
        global_dist="results/coverage/{sample}.mosdepth.global.dist.txt",
        summary="results/coverage/{sample}.mosdepth.summary.txt"
    params:
        prefix="results/coverage/{sample}",
        window=WINDOW_SIZE
    log:
        "logs/align/mosdepth/{sample}.log"
    priority: 10
    threads: 4
    shell:
        "mosdepth --no-per-base --by {params.window} --threads {threads} {params.prefix} {input.bam} 2> {log}"

# Extract the 2Mb region of interest and save as final BAM
# rule extract_region:
#     input:
#         bam="results/temp_bams/{sample}.sorted.bam",
#         bai="results/temp_bams/{sample}.sorted.bam.bai"
#     output:
#         bam="results/region_bams/{sample}.bam"
#     params:
#         region=REGION
#     threads: 4
#     log:
#         "logs/extract_region/{sample}.log"
#     shell:
#         "samtools view -b -h -@ {threads} -o {output.bam} {input.bam} {params.region} 2> {log}"