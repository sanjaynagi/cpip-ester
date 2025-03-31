#!/usr/bin/env python3
"""
Calculate modal CNVs for all genes in the GFF file
"""

import pandas as pd
import numpy as np
import allel
import os
from collections import Counter

# Get inputs, outputs, and parameters from Snakemake
cnv_files = snakemake.input.cnv_files
annotation_file = snakemake.input.annotation
modal_cnvs_file = snakemake.output.modal_cnvs
samples = snakemake.params.samples

# Load GFF file using scikit-allel
print("Loading GFF file...")
gff_df = allel.gff3_to_dataframe(annotation_file, attributes=['ID']).sort_values(['seqid', 'start'])

# Filter for gene features
print("Extracting all genes from GFF file...")
genes_df = gff_df[gff_df['type'] == 'protein_coding_gene']
print(f"Found {len(genes_df)} genes in the GFF file")

def calculate_modal_copy_number(cnv_data, chrom, start, end):
    """
    Calculate the modal copy number for a region
    
    Parameters:
    cnv_data: DataFrame with CNV calls
    chrom: Chromosome name
    start: Start position
    end: End position
    
    Returns:
    The most common copy number in the region
    """
    # Filter CNV data for the region
    region_cnvs = cnv_data[
        (cnv_data['Chromosome'] == chrom) & 
        (cnv_data['End'] > start) & 
        (cnv_data['Start'] < end)
    ]
    
    if len(region_cnvs) == 0:
        return 2  # Assume diploid if no CNVs in the region
    
    # For overlapping CNVs, calculate weighted copy number by overlap size
    weighted_copy_numbers = []
    
    for _, cnv in region_cnvs.iterrows():
        # Calculate overlap
        overlap_start = max(start, cnv['Start'])
        overlap_end = min(end, cnv['End'])
        overlap_size = overlap_end - overlap_start
        
        # Only count if there's actual overlap
        if overlap_size > 0:
            weighted_copy_numbers.extend([cnv['CopyNumber']] * overlap_size)
    
    if not weighted_copy_numbers:
        return 2  # Assume diploid if no overlap
    
    # Get modal copy number
    copy_number_counts = Counter(weighted_copy_numbers)
    modal_cn = copy_number_counts.most_common(1)[0][0]
    
    return modal_cn

chrom, positions = 'CM027412.1:135390000-139390000'.split(':')
region_start = int(positions.split('-')[0])

# Process each CNV file and gene
results = []

for i, cnv_file in enumerate(cnv_files):
    sample = samples[i]
    print(f"Processing sample: {sample}")
    
    # Load CNV data
    cnv_data = pd.read_csv(cnv_file)

    cnv_data['Start'] = cnv_data['Start'] + region_start
    cnv_data['End'] = cnv_data['End'] + region_start
    cnv_data['Chromosome'] = chrom
    
    # Process each gene
    for _, gene in genes_df.iterrows():
        gene_name = gene['ID']
        gene_chrom = gene['seqid']
        # GFF is 1-based, so we don't need to adjust
        gene_start = gene['start']
        gene_end = gene['end']
        
        # Calculate modal copy number
        modal_cn = calculate_modal_copy_number(cnv_data, gene_chrom, gene_start, gene_end)
        
        # Store result
        results.append({
            'Sample': sample,
            'Region': gene_name,
            'Chromosome': gene_chrom,
            'Start': gene_start,
            'End': gene_end,
            'ModalCopyNumber': modal_cn,
            'Length': gene_end - gene_start
        })

# Create and save results
results_df = pd.DataFrame(results)
results_df.to_csv(modal_cnvs_file, index=False)

print(f"Modal CNV analysis complete. Results saved to {modal_cnvs_file}")
