# GenotypeFilters
This repository contains scripts and reference data for filtering genotype data. 

## SNPplatforms
This section describes a filtering script for selecting the SNPs matching a particular genotype array, e.g. Illumina's Infinium Omni 2.5Exome v1.5 array. The script is run using bash shell, R and the Genome Analysis Toolkit (GATK). The genotype array is given as a gzipped interval list, which is provided in the folder in .gz format. 

### Requirements and recommendations
- R >=4.3.0 (with tidyverse)
- GATK >=4.4.0.0 
- Java v17 is required to run GATK v4.4.0.0
- Interval lists downloaded to match your respective reference genome (b37 or b38)
- The input vcf files must follow the vcf format >=v4.1 and be accompanied by an index. This can be created e.g. by ```bcftools index -t <path to file>.vcf.gz```.

*Note:* GATK is run using the wrapper script as provided by GATK, meaning that the java options are passed through as described in this [GATK documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360035531892-GATK4-command-line-syntax). 

*NB!* ```elegant_filter_gatk.sh``` changed name to: ```wgs2snp_array.sh```

### Usage
```
usage: ./wgs2snp_array.sh [--help] --input-vcf <file-path> [OTHER OPTIONS] --out <file-path>

Usage:
	-i  | --input-vcf	- Variant call format input file (.vcf or .vcf.gz).
	-pf | --prefix-chr	- Prefix on the chromosome column of the input-vcf. Default: "".
	-hg | --human-genome	- Specifies the human genome reference version to use. Accepts either 37 or 38. Default: 37.
	-b  | --positions	- Path to interval list of positions to include in gzipped format.
				  The file contains columns snpID, rsID (if available) and position given as zero-indexed.
				  Be sure to download the file matching your reference genome. 
	-p  | --padding		- Number of bases to add padding around the included intervals. Default: 1.
	-w  | --no-wrapper	- Flag to deactivate the use of the GATK wrapper script, which is recommended way to run the toolkit. Default: Wrapper is used by default.
	-g  | --gatk-path	- Path to the gatk toolkit, if the wrapper is not used, ".jar" file ending is added automatically. Default: gatk.
	-m  | --memory		- Memory allocation in G to use for java-specific arguments. Default: 4.
	-o  | --out		- Output path for the filtered file.
	-h  | --help		- Print Usage message.


```

### Interval files 
The interval lists needed for the ```--positions``` flag are located in ```GenotypeFilters/SNPplatforms/interval_lists```. Here, interval lists are available for two SNP arrays, namely Illumina's Infinium Omni2.5Exome v1.5 and Global diversity array (GDA) v1.0. A combined version of unique positions across these two arrays is also available. 
The positions have been made available for two versions of the human genome reference: b37 \& b38. Make sure to use the positions list file that matches your version of the human genome used for calling variants. 

| Array name                   | Version | \# of sites | File name                                       |
|:-----------------------------|:--------|:------------|:------------------------------------------------|
| Infinium Omni2.5Exome        | 1.5     |  ~2.5 mio.  | ```ilm_omni2.5exome_v1.5_intervals_b*.txt.gz``` |
| Global diversity array (GDA) | 1.0     |  ~1.8 mio.  | ```gda_v1.0_intervals_b*.txt.gz```              |
| Combined O2.5Ex \& GDA       |         |  ~3.6 mio.  | ```combined_ilm_omni2.5exome_v1.5_gda_v1.0_intervals_b3*.txt.gz``` |

Please note, that the provided interval lists were filtered to exclude the ACMG variants classified to be reported back (ACMG v3.2, 2023). The list of ACMG variants was created by linking the ACMG genes and classification to ClinVar. 
![acmg](ACMG_variants.png)



## ACMG filter (deprecated)
Initially, the ACMG filter has been created. This will take in vcf format genotype data files and exclude all variants that match positions in the accompanying ```clinvar_acmg3.2_subset.txt.gz```. The positions are listed with both assembly GRCh37 and 38, but the python program ```match_filter_vcf.py``` has the flag ```--human-genome```, where one of the two can be specified. 


### Requirements and recommendations
The ACMGfilter script requires 
- python >=3.9
- pandas >=2
- pyarrow backend for pandas v2

It is recommended to subset the input vcf file by chromosome and run the chromosomes in parallel. The input vcf file must have a header (```['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', ... samples]```) and optionally the vcf format info lines started by "##". 


### Usage 
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
                        Path to input variant call format (.vcf) file containing variants to
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







