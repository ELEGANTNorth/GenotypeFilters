# GenotypeFilters

This repository contains scripts and reference data for filtering genotype data. 

## ACMG filter
Initially, the ACMG filter has been created. This will take in vcf format genotype data files and exclude all variants that match positions in the accompanying ```clinvar_acmg3.2_subset.txt.gz```. The positions are listed with both assembly GRCh37 and 38, but the python program ```match_filter_vcf.py``` has the flag ```--human-genome```, where one of the two can be specified. 

The list of ACMG variants was created by linking the ACMG variants to ClinVar. 
![acmg](ACMG_variants.png)

It is recommended to subset the input vcf file by chromosome and run the chromosomes in parallel. The input vcf file must have a header (```['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', ... samples]```) and optionally the vcf format info lines started by "##". 

```
usage: match_filter_vcf.py [-h] --input-vcf INPUT_FILE
                           [--human-genome {37,38,37,38}]
                           [--acmg-list ACMG_PATH] [--snps-only] --output-vcf
                           OUTPUT_FILE

Program to filter genotype data (vcf format files) for ACMG actionable
variants.

optional arguments:
  -h, --help            show this help message and exit
  --input-vcf INPUT_FILE
                        Path to input vcf format file containing variants to
                        be filtered.
  --human-genome {37,38,37,38}
                        Specifies the human genome reference version to use.
                        Accepts either 37 or 38. Default is 37.
  --acmg-list ACMG_PATH
                        Path to the file containing the ACMG v3.2 variants to
                        be excluded. The file must be long format and contain
                        positions in both GRCh37 and 38 coordinates
                        (colname=Assembly). Furthermore it must contrain the
                        columns Assembly, Chromosome and PositionVCF.
  --snps-only           If specified, only SNPs will be kept in the filtered
                        file.
  --output-vcf OUTPUT_FILE
                        Path to output vcf format file containing filtered
                        variants.
```
