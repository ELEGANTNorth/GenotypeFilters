#!/usr/bin/bash

# Bash shell program to run the ELEGANT filters using GATK toolkit
# Last updated: June 12th, 2024

## SET DEFAULTS 
GATKPATH=gatk
PREFIXCHR=""
HG=37
POS=gatk_intervals_b37.txt.gz
PADDING=1
MEM=4

##################
#   FUNCTIONS
##################
# Usage function 
function usage()
{
    echo -e "\nUsage:"
    echo -e "\t-i  | --input-vcf\t\t- Variant call format input file (.vcf or .vcf.gz)."
    echo -e "\t-pf | --prefix-chr\t\t- Prefix on the chromosome column of the input-vcf. Default: \"\"."
    echo -e "\t-hg | --human-genome\t\t- Specifies the human genome reference version to use. Accepts either 37 or 38. Default: 37."
    echo -e "\t-b  | --positions\t\t- Path to interval list of positions to include in gzipped format. The file contains columns snpID, rsID (if available) and position given as zero-indexed. Be sure to download the file matching your reference genome. "
    echo -e "\t-p  | --padding\t\t- Number of bases to add padding around the included intervals. Default: 1."
    echo -e "\t-g  | --gatk-path\t\t- Path to the gatk toolkit. Default: gatk."
    echo -e "\t-m  | --memory\t\t- Memory allocation in G to use for java-specific arguments. Default: 4."
    echo -e "\t-o  | --out\t\t- Output path for the filtered file."
    echo -e "\t-h  | --help\t\t- Print Usage message. "
}

###########################
# Loop through arguments
for arg in "$@"
do 
    case $arg in 
        -h | --help)    # print help message and exit
            usage
            exit 0
            ;;
        -i | --input-vcf)
            INVCF="$2"
            shift   # Remove argument name from processing 
            shift   # Remove argument value from processing 
            ;;
        -pf | --prefix-chr)
            PREFIXCHR="$2"
            shift   # Remove argument name from processing 
            shift   # Remove argument value from processing 
            ;;
        -hg | --human-genome)
            HG="$2"
            shift   # Remove argument name from processing 
            shift   # Remove argument value from processing 
            ;;
        -b | --positions)
            POS="$2"
            shift   # Remove argument name from processing 
            shift   # Remove argument value from processing 
            ;;
        -p | --padding)
            PADDING="$2"
            shift   # Remove argument name from processing 
            shift   # Remove argument value from processing 
            ;;
        -g | --gatk-path)
            GATKPATH="$2"
            shift   # Remove argument name from processing 
            shift   # Remove argument value from processing 
            ;;
        -m | --memory)
            MEM="$2"
            shift   # Remove argument name from processing 
            shift   # Remove argument value from processing 
            ;;
        -o | --out-vcf)
            OUTVCF="$2"
            shift   # Remove argument name from processing 
            shift   # Remove argument value from processing 
            ;;
    esac
done

# Status on applied settings
echo -e "Settings applied: \n\tinput-vcf: $INVCF \n\tHuman genome: GRCh${HG} \n\tPositions filter file: $POS \n\tGATK path: $GATKPATH \n\tMemory allocation: ${MEM}G \n\tOut path: $OUTVCF \n\n"

# GET THE SUBSET OF POSITIONS from the first column 
gzip -dc $POS | cut -d" " -f1 | tail -n +2 > $POS.GRCh${HG}.list

# Get the contig lengths
if [[ "$INVCF" == *.gz ]]; then
    echo "this file is zipped"
    Rscript -e "con <- gzfile('$INVCF');count <- 1;flag <- TRUE;while (flag) {line <- scan(con, skip = count-1, nlines = 1, what = 'character', sep=NULL, quiet=TRUE)[1];if (!startsWith(line, '##')) { flag <- FALSE } else if (startsWith(line, '##contig=')) { cat(line, '\n') };count <- count + 1; }; close(con)" > contig_lengths.txt
else 
    echo "this file is not zipped"
    Rscript -e "con <- file('$INVCF');count <- 1;flag <- TRUE;while (flag) {line <- scan(con, skip = count-1, nlines = 1, what = 'character', sep=NULL, quiet=TRUE)[1];if (!startsWith(line, '##')) { flag <- FALSE } else if (startsWith(line, '##contig=')) { cat(line, '\n') }; count <- count + 1; }; close(con)" > contig_lengths.txt
fi

# Check the contig lengths can contain all the positions needed 
# options(scipen = 999); - Disables printing results in scientific notation.
Rscript -e "library(tidyverse); options(scipen = 999); data1 <- readLines('contig_lengths.txt', warn = FALSE); dict1 <- list(); for (entry in data1) { if (!startsWith(entry, '##')) { break } else if (startsWith(entry, '##contig=')) { key <- sub('.*ID=${PREFIXCHR}([0-9]+),.*', '\\\1', entry); value <- as.integer(sub('.*length=(\\\d+)>.*', '\\\1', entry)); dict1[[key]] <- value }; }; positions <- read_delim('$POS.GRCh${HG}.list', delim=':', col_names=c('CHR', 'POS')) %>% separate_wider_delim(POS, delim='-', names=c('start', 'stop')) %>% mutate(stop = as.integer(stop)); rows <- data.frame(); for (key in names(dict1)){ key_set <- dplyr::filter(positions, CHR == key); rows <- bind_rows(rows, dplyr::filter(key_set, stop < dict1[[key]]))}; rows %>% mutate(start = as.character(as.integer(start)), stop = as.character(stop)) %>% unite(start, stop, col = 'POS', sep = '-') %>% unite(CHR, POS, col = 'X', sep = ':') %>% write_delim('$POS.GRCh${HG}_filtered.list', delim=' ', col_names=FALSE)"

# clean-up 
rm contig_lengths.txt
rm $POS.GRCh${HG}.list

# RUN VARIANT SELECTION
$GATKPATH --java-options "-Xmx${MEM}G" SelectVariants --variant $INVCF --L $POS.GRCh${HG}_filtered.list --interval-padding $PADDING --output $OUTVCF
