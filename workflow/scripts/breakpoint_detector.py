#!/usr/bin/env python3
"""
Detect CNV breakpoints using soft-clipped reads in the BAM file
"""

import pysam
import pandas as pd
import numpy as np

# Get inputs, outputs, and parameters from Snakemake
bam_file = snakemake.input.bam
cnv_calls_file = snakemake.input.cnv_calls
breakpoints_file = snakemake.output.breakpoints
pre_clipping_fastq = snakemake.output.pre_clipping_fastq
post_clipping_fastq = snakemake.output.post_clipping_fastq
min_mapq = int(snakemake.params.min_mapq)

print(f"Detecting breakpoints in BAM file: {bam_file}")
print(f"Using CNV calls from: {cnv_calls_file}")
print(f"Minimum mapping quality: {min_mapq}")

# Load CNV calls
cnv_calls = pd.read_csv(cnv_calls_file)

# Open the BAM file
alignment = pysam.AlignmentFile(bam_file, 'rb')

# Open output FASTQ files
pre_clipping_fastq_file = open(pre_clipping_fastq, 'w')
post_clipping_fastq_file = open(post_clipping_fastq, 'w')

# Function to collect soft-clipped reads in a genomic region
def collect_soft_clipped_reads(alignment, chromosome, start, end, min_clip_length=10):
    """
    Collect soft-clipped reads in a given genomic region
    
    Parameters:
    alignment: pysam.AlignmentFile object
    chromosome: Chromosome name
    start: Start position
    end: End position
    min_clip_length: Minimum length of soft-clipped sequence to consider
    
    Returns:
    Dictionary of breakpoints with their supporting reads
    """
    print(f"Searching for breakpoints in {chromosome}:{start}-{end}")
    
    # Fetch alignments in the region
    alignments = alignment.fetch(chromosome, start, end)
    
    # Initialize counters and containers
    soft_clipping_start_points = []  # Clipping at the end of read
    soft_clipping_end_points = []    # Clipping at the start of read
    
    for read in alignments:
        # Skip unmapped reads
        if read.is_unmapped:
            continue
            
        # Skip reads with empty positions vector (unlikely but possible)
        if len(read.positions) == 0:
            continue
            
        # Get CIGAR string
        cigar = read.cigar
        
        # Skip if there's only one CIGAR operation (can't be soft-clipped)
        if len(cigar) == 1:
            continue
            
        # Check for soft-clipping at the start of the read
        soft_clipping_end_point = None
        if (cigar[0][0] == 4) and (cigar[0][1] >= min_clip_length):
            soft_clipping_end_point = read.positions[0]
            
            # Write soft-clipped sequence to FASTQ
            pre_clipping_fastq_file.write(f'@{read.query_name}_left_of_pos_{soft_clipping_end_point}\n')
            pre_clipping_fastq_file.write(f'{read.query_sequence[:cigar[0][1]]}\n')
            pre_clipping_fastq_file.write('+\n')
            pre_clipping_fastq_file.write(f'{read.query_qualities[:cigar[0][1]].tobytes().decode("ascii")}\n')
        
        # Check for soft-clipping at the end of the read
        soft_clipping_start_point = None
        if (cigar[-1][0] == 4) and (cigar[-1][1] >= min_clip_length):
            # The last position in the reference + 1 (first position of soft-clipping)
            soft_clipping_start_point = read.positions[-1] + 1
            
            # Write soft-clipped sequence to FASTQ
            post_clipping_fastq_file.write(f'@{read.query_name}_right_of_pos_{soft_clipping_start_point}\n')
            clipped_seq_start = len(read.query_sequence) - cigar[-1][1]
            post_clipping_fastq_file.write(f'{read.query_sequence[clipped_seq_start:]}\n')
            post_clipping_fastq_file.write('+\n')
            post_clipping_fastq_file.write(f'{read.query_qualities[clipped_seq_start:].tobytes().decode("ascii")}\n')
        
        # Store breakpoint info if mapping quality is sufficient
        if read.mapping_quality >= min_mapq:
            if soft_clipping_start_point:
                soft_clipping_start_points.append({
                    'Position': soft_clipping_start_point,
                    'Clipped_sequence': read.query_sequence[clipped_seq_start:],
                    'Type': 'start'
                })
            if soft_clipping_end_point:
                soft_clipping_end_points.append({
                    'Position': soft_clipping_end_point,
                    'Clipped_sequence': read.query_sequence[:cigar[0][1]],
                    'Type': 'end'
                })
    
    return soft_clipping_start_points, soft_clipping_end_points

# Process each CNV region to find breakpoints
breakpoints = []

for _, cnv in cnv_calls.iterrows():
    # Skip if copy number is 2 (no CNV)
    if cnv['CopyNumber'] == 2:
        continue
        
    # Look for breakpoints in the region
    padding = 500  # Look 500bp on each side of the CNV region
    start_points, end_points = collect_soft_clipped_reads(
        alignment, 
        cnv['Chromosome'], 
        max(0, cnv['Start'] - padding), 
        cnv['End'] + padding
    )
    
    # Associate breakpoints with CNV
    for bp in start_points + end_points:
        breakpoints.append({
            'Chromosome': cnv['Chromosome'],
            'Position': bp['Position'],
            'Type': bp['Type'] + ' breakpoint',
            'Associated_CNV_start': cnv['Start'],
            'Associated_CNV_end': cnv['End'],
            'Copy_number': cnv['CopyNumber']
        })

# Close FASTQ files
pre_clipping_fastq_file.close()
post_clipping_fastq_file.close()

# Create and save breakpoints dataframe
breakpoints_df = pd.DataFrame(breakpoints)
if len(breakpoints_df) > 0:
    breakpoints_df.to_csv(breakpoints_file, index=False)
else:
    # Create empty dataframe with correct columns if no breakpoints found
    empty_df = pd.DataFrame(columns=[
        'Chromosome', 'Position', 'Type', 
        'Associated_CNV_start', 'Associated_CNV_end', 'Copy_number'
    ])
    empty_df.to_csv(breakpoints_file, index=False)

print(f"Breakpoint detection complete. Found {len(breakpoints)} potential breakpoints.")