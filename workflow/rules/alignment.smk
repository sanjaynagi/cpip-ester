# Align reads to whole chromosome reference - updated for SRR processing
rule align_reads:
    input:
        r1="results/dl/temp_reads/{sample}_{srr}_1.fastq.gz",
        r2="results/dl/temp_reads/{sample}_{srr}_2.fastq.gz",
        ref=REF_GENOME
    output:
        pipe("results/pipes/{sample}_{srr}.sam")
    log:
        "logs/align/bwa/{sample}_{srr}.log"
    params:
        rg=lambda wildcards: f"@RG\\tID:{wildcards.srr}\\tSM:{wildcards.sample}"
    threads: 8
    shell:
        "resources/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t {threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2} > {output} 2> {log}"

# Convert SAM to sorted BAM - updated for SRR processing
rule sam_to_sorted_bam:
    input:
        sam="results/pipes/{sample}_{srr}.sam"
    output:
        bam="results/srr_bams/{sample}_{srr}.sorted.bam"
    threads: 4
    log:
        "logs/align/sam_to_bam/{sample}_{srr}.log"
    shell:
        """
        samtools view -b -h -F 4 -@ {threads} {input.sam} | \
        samtools sort -@ {threads} -o {output.bam} - 2> {log}
        """

# Index the SRR BAM files
rule index_srr_bam:
    input:
        bam="results/srr_bams/{sample}_{srr}.sorted.bam"
    output:
        bai="results/srr_bams/{sample}_{srr}.sorted.bam.bai"
    log:
        "logs/align/bam_index_srr/{sample}_{srr}.log"
    shell:
        "samtools index {input.bam} 2> {log}"

rule merge_sample_bams_and_region:
    input:
        bams=lambda wildcards: expand("results/srr_bams/{sample}_{srr}.sorted.bam", 
                                      sample=wildcards.sample, 
                                      srr=sample_to_srr[wildcards.sample]),
        bais=lambda wildcards: expand("results/srr_bams/{sample}_{srr}.sorted.bam.bai", 
                                      sample=wildcards.sample, 
                                      srr=sample_to_srr[wildcards.sample])
    output:
        bam="results/region_bams/{sample}.sorted.bam"
    threads: 4
    log:
        "logs/align/merge_bams/{sample}.log"
    params:
        region=REGION,
        # If only one BAM file, use direct samtools view instead of merge+view
        cmd=lambda wildcards, input: (
            "samtools view -b" if len(input.bams) == 1 
            else "samtools merge -f -"
        )
    shell:
        """
        if [ "{params.cmd}" = "samtools view -b" ]; then
            {params.cmd} {input.bams[0]} {params.region} > {output.bam} 2> {log}
        else
            {params.cmd} -@ {threads} {input.bams} 2> {log} | \
            samtools view -b -@ {threads} - {params.region} > {output.bam} 2>> {log}
        fi
        """

# Index sorted BAM
rule index_final_merged_bam:
    input:
        bam="results/region_bams/{sample}.sorted.bam"
    output:
        bai="results/region_bams/{sample}.sorted.bam.bai"
    log:
        "logs/bam_final/{sample}.log"
    shell:
        "samtools index {input.bam} 2> {log}"


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
