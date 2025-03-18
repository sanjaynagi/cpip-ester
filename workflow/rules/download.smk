# Get download links using ffq with rate limiting - updated for SRR processing
rule get_download_links:
    output:
        json="results/ffq/{sample}_{srr}.json"
    params:
        srr=lambda wildcards: wildcards.srr,
        sleep=lambda x: f"{1 + random.randint(0, 3)}"
    log:
        "logs/dl/dl_link/{sample}_{srr}.log"
    resources:
        ff1=1
    retries: 3
    shell:
        """
        # Sleep to avoid rate limiting
        sleep {params.sleep}
        ffq --ftp {params.srr} > {output.json} 2> {log}
        """

# Extract URLs from ffq output - updated for SRR processing
rule extract_urls:
    input:
        json="results/ffq/{sample}_{srr}.json"
    output:
        r1_url="results/urls/{sample}_{srr}_1.url",
        r2_url="results/urls/{sample}_{srr}_2.url"
    log:
        "logs/dl/extract_urls/{sample}_{srr}.log"
    shell:
        """
        jq -r '.[0].url' {input.json} > {output.r1_url} 2> {log}
        jq -r '.[1].url' {input.json} > {output.r2_url} 2>> {log}
        """

# Download and decompress R1 reads - updated for SRR processing
rule download_r1:
    input:
        url="results/urls/{sample}_{srr}_1.url"
    output:
        temp("results/temp_reads/{sample}_{srr}_1.fastq.gz")
    log:
        "logs/dl/dl_r1/{sample}_{srr}.log"
    resources:
        dl=1
    retries: 3
    shell:
        "curl $(cat {input.url}) > {output} 2> {log}"

# Download and decompress R2 reads - updated for SRR processing
rule download_r2:
    input:
        url="results/urls/{sample}_{srr}_2.url"
    output:
        temp("results/temp_reads/{sample}_{srr}_2.fastq.gz")
    log:
        "logs/dl/dl_r2/{sample}_{srr}.log"
    resources:
        dl=1
    retries: 3
    shell:
        "curl $(cat {input.url}) > {output} 2> {log}"