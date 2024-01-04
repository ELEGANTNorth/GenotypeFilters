#!/usr/bin/python3
# Last updated: December 28th, 2023 

# Load packages
import argparse
import os 
import sys
import pandas as pd 


######################################
#   CREATE THE COMMANDLINE OPTIONS
######################################
parser = argparse.ArgumentParser(description='Program to filter genotype data (vcf format files) for ACMG actionable variants.')
parser.add_argument('--input-vcf', dest='input_file', type=str, required=True, 
                    help='Path to input vcf format file containing variants to be filtered.')
parser.add_argument('--human-genome', dest='hg', required=False, choices=['37', '38', 37, 38], default='37', 
                    help='Specifies the human genome reference version to use. Accepts either 37 or 38. Default is 37.')
parser.add_argument('--acmg-list', dest='acmg_path', required=False, type=str, default='clinvar_acmg3.2_subset.txt.gz', 
                    help='Path to the file containing the ACMG v3.2 variants to be excluded. The file must be long format and contain positions in both GRCh37 and 38 coordinates (colname=Assembly). Furthermore it must contrain the columns Assembly, Chromosome and PositionVCF.')
parser.add_argument('--snps-only', dest='snps_only', action='store_true', default=False, 
                    help='If specified, only SNPs will be kept in the filtered file.')
parser.add_argument('--output-vcf', dest='output_file', type=str, required=True, 
                    help='Path to output vcf format file containing filtered variants.')

args = parser.parse_args()
print(f'Arguments:\n{args}')


############################
#       MAIN PROGRAM
############################
# The vcf file format includes a set of info lines starting with "##". 
# These info lines are also included in the new filtered output file. 
filesize = os.stat(args.input_file).st_size/10**6     # size in megabytes
print(f'Size of input vcf file is: {filesize:.2f} megabytes')

# Make a smaller filter for ACMG list to match the human genome reference 
acmg_filter = pd.read_table(args.acmg_path, engine='pyarrow', dtype_backend='pyarrow', dtype = {'PositionVCF': 'object', 'Chromosome': 'object'}) #dtype='object')
acmg_filter = acmg_filter.loc[acmg_filter.Assembly == f'GRCh{args.hg}']
print('\tColumns of ACMG filter data:', acmg_filter.columns, '\n\tWith shape:', acmg_filter.shape)

# Check if the file was gzipped 
# File can now be opened and read line by line 
flag = False
#kept = 0
removed = 0
totalLines = 0

print('\n# Data is not gzipped.')
infile = open(args.input_file, 'r')
outfile = open(args.output_file, 'w')

if args.snps_only: 
    # Iterate through the lines of the file and check for chr number to ensure smaller subset of the acmg set is needed 
    for line0 in infile: 
        # Now a line was found
        if flag: 
            totalLines += 1
            line = line0[:-1].split()

            # CHeck if new chromosome
            if CHR != line[idx_chr]:
                acmg_sub = acmg_filter.loc[acmg_filter.Chromosome == line[idx_chr]]
                print("New ACMG subset:", acmg_sub.shape)

                # Make variables ready for next line 
                CHR = line[idx_chr]

            # Check if position matches
            match = acmg_sub.loc[acmg_sub.PositionVCF == int(line[idx_pos])]
            # if line[idx_pos] in acmg_sub.PositionVCF: 
            if match.shape[0] > 0: 
                removed += 1
                continue
            elif len(line[idx_ref]) > 2 or len(line[idx_alt]) > 2: 
                removed += 1
                continue
            else: 
                # kept += 1
                outfile.write(line0)

        # Simply re-write the lines of info to the filtered out
        elif line0[:2] == "##": 
            outfile.write(line0)
        else: 
            # The header was found
            flag = True
            header = line0[:-1].split()
            print('\n### Number of columns in file:', len(header))

            # VCF line format: ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', ... samples]
            print('Info fields:', header[:10])

            # Get indexes to save 
            idx_chr = header.index('#CHROM')
            idx_pos = header.index('POS')
            idx_ref = header.index('REF')
            idx_alt = header.index('ALT')
            CHR = "0"

            # Write the header to out
            outfile.write(line0)
else: 
    # Iterate through the lines of the file and check for chr number to ensure smaller subset of the acmg set is needed 
    for line0 in infile: 
        # Now a line was found
        if flag: 
            totalLines += 1
            line = line0[:-1].split()

            # CHeck if new chromosome
            if CHR != line[idx_chr]:
                acmg_sub = acmg_filter.loc[acmg_filter.Chromosome == line[idx_chr]]
                print("New ACMG subset:", acmg_sub.shape)

                # Make variables ready for next line 
                CHR = line[idx_chr]

            # Check if position matches
            match = acmg_sub.loc[acmg_sub.PositionVCF == int(line[idx_pos])]
            # if line[idx_pos] in acmg_sub.PositionVCF: 
            if match.shape[0] > 0: 
                removed += 1
                continue
            else: 
                # kept += 1
                outfile.write(line0)

        # Simply re-write the lines of info to the filtered out
        elif line0[:2] == "##": 
            outfile.write(line0)
        else: 
            # The header was found
            flag = True
            header = line0[:-1].split()
            print('\n### Number of columns in file:', len(header))

            # VCF line format: ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', ... samples]
            print('Info fields:', header[:10])

            # Get indexes to save 
            idx_chr = header.index('#CHROM')
            idx_pos = header.index('POS')
            idx_ref = header.index('REF')
            idx_alt = header.index('ALT')
            CHR = "0"

            # Write the header to out
            outfile.write(line0)

# Close the file 
infile.close()
outfile.close()

# Print status and end 
kept = totalLines - removed 
print(f'\n### Status of filtering:\n\tNumber of variants: {totalLines}\n\tNumber of variants kept: {kept}\n\tNumber of variants removed: {removed}')
print("\n### DONE WITH FILTERS!")

