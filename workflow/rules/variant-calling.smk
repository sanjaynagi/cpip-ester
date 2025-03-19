# Create reference dictionary and index files
rule prepare_reference:
    input:
        ref=REF_GENOME
    output:
        dict=REF_GENOME.replace(".fa", ".dict"),
    log:
        "logs/call/prepare_reference/prepare.log"
    shell:
        """
        gatk CreateSequenceDictionary -R {input.ref} -O {output.dict} 2> {log}
        """

# Mark duplicates in the region BAM files
rule mark_duplicates:
    input:
        bam="results/region_bams/{sample}.sorted.bam"
    output:
        bam="results/region_bams_md/{sample}.md.bam",
        metrics="results/region_bams_md/{sample}.metrics.txt"
    log:
        "logs/mark_duplicates/{sample}.log"        # Create sequence dictionary

    resources:
        mem_mb=16000
    shell:
        """
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY LENIENT 2> {log}
        """

rule index_bam_md:
    input:
        bam="results/region_bams_md/{sample}.md.bam",
    output:
        bai="results/region_bams_md/{sample}.md.bam.bai",
    log:
        "logs/call/index_bam_md/{sample}.log"
    shell:
        """
        samtools index {input.bam} 2> {log}
        """
   

# Call variants with HaplotypeCaller directly on duplicate-marked BAMs
rule haplotype_caller:
    input:
        bam="results/region_bams_md/{sample}.md.bam",
        bai="results/region_bams_md/{sample}.md.bam.bai",
        ref=REF_GENOME,
        ref_dict=REF_GENOME.replace(".fa", ".dict"),
        ref_fai=REF_GENOME + ".fai"
    output:
        gvcf="results/vcf/gvcfs/{sample}.g.vcf.gz",
        idx="results/vcf/gvcfs/{sample}.g.vcf.gz.tbi"
    # params:
    #     region=REGION
    log:
        "logs/call/haplotype_caller/{sample}.log"
    resources:
        mem_mb=8000
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            --emit-ref-confidence GVCF \
            --create-output-variant-index true 2> {log}
        """

        #            -L {params.region} \


# Combine GVCFs from all samples
rule combine_gvcfs:
    input:
        gvcfs=expand("results/vcf/gvcfs/{sample}.g.vcf.gz", sample=SAMPLE_NAMES),
        ref=REF_GENOME,
        ref_dict=REF_GENOME.replace(".fa", ".dict"),
        ref_fai=REF_GENOME + ".fai"
    output:
        gvcf="results/vcf/combined/combined.g.vcf.gz",
        idx="results/vcf/combined/combined.g.vcf.gz.tbi",
    params:
        gvcfs=lambda wildcards, input: [f"-V {vcf}" for vcf in input.gvcfs],
        # region=REGION
    log:
        "logs/call/combine_gvcfs/combine.log"
    resources:
        mem_mb=16000
    shell:
        """
        gatk CombineGVCFs \
            -R {input.ref} \
            {params.gvcfs} \
            -O {output.gvcf} 2> {log}
        """

    # -L {params.region}

# Genotype the combined GVCF
rule genotype_gvcfs:
    input:
        gvcf="results/vcf/combined/combined.g.vcf.gz",
        idx="results/vcf/combined/combined.g.vcf.gz.tbi",
        ref=REF_GENOME,
        ref_dict=REF_GENOME.replace(".fa", ".dict"),
        ref_fai=REF_GENOME + ".fai"
    output:
        vcf="results/vcf/raw_variants.vcf.gz",
        idx="results/vcf/raw_variants.vcf.gz.tbi"
    # params:
    #     region=REGION
    log:
        "logs/call/genotype_gvcfs/genotype.log"
    resources:
        mem_mb=16000
    shell:
        """
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {input.gvcf} \
            -O {output.vcf} \
            --create-output-variant-index true 2> {log}
        """

#            -L {params.region} \




# Extract only SNPs
rule extract_snps:
    input:
        vcf="results/vcf/raw_variants.vcf.gz",
        ref=REF_GENOME
    output:
        vcf="results/vcf/raw_snps.vcf.gz",
        idx="results/vcf/raw_snps.vcf.gz.tbi"
    log:
        "logs/call/extract_snps/extract.log"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            --select-type-to-include SNP 2> {log}
        """

# Filter SNPs
rule filter_snps:
    input:
        vcf="results/vcf/raw_snps.vcf.gz",
        ref=REF_GENOME
    output:
        vcf="results/vcf/filtered_snps.vcf.gz",
        idx="results/vcf/filtered_snps.vcf.gz.tbi"
    log:
        "logs/call/filter_snps/filter.log"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \
            --filter-name "snp_filter" 2> {log}
        """