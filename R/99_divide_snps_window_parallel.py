#!/bin/python3

import pandas as pd
from multiprocessing import Pool, cpu_count
import argparse

def process_chromosome(args):
    """Process a single chromosome to assign SNPs to windows."""
    chrom, chrom_snps, window_size = args
    chrom_snps = chrom_snps.copy()
    chrom_snps['window'] = chrom_snps['pos'] // window_size
    chrom_snps['window_start'] = chrom_snps['window'] * window_size
    chrom_snps['window_end'] = (chrom_snps['window'] + 1) * window_size - 1
    return chrom_snps

def divide_snps_parallel(snp_file, window_size, output_file):
    """Divide SNPs into equal-sized windows and process chromosomes in parallel."""
    # Load the SNP dataset
    snps = pd.read_csv(snp_file, sep='\t')
    
    # Ensure the SNP file contains the required columns
    if not {'chrom', 'pos'}.issubset(snps.columns):
        raise ValueError("SNP file must contain 'chrom' and 'pos' columns.")
    
    # Group SNPs by chromosome
    grouped = list(snps.groupby('chrom'))
    
    # Prepare arguments for multiprocessing
    args = [(chrom, group, window_size) for chrom, group in grouped]
    
    # Determine the number of processes to use
    num_processes = min(cpu_count(), len(grouped))
    
    # Process each chromosome in parallel
    with Pool(num_processes) as pool:
        result = pool.map(process_chromosome, args)
    
    # Concatenate results from all processes
    result_df = pd.concat(result, ignore_index=True)
    
    # Save to the output file
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"Results saved to {output_file}")

def main():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Divide SNPs into equal-sized windows by chromosome.")
    parser.add_argument("snp_file", help="Path to the input SNP file (TSV format).")
    parser.add_argument("window_size", type=int, help="Size of each window.")
    parser.add_argument("output_file", help="Path to the output file (TSV format).")
    args = parser.parse_args()
    
    # Call the function with the provided arguments
    divide_snps_parallel(args.snp_file, args.window_size, args.output_file)

if __name__ == "__main__":
    main()

