#!/usr/bin/bash

# Bash shell program to run the ELEGANT filters using GATK toolkit
# Last updated: January 25th, 2024

## SET DEFAULTS 
GATKPATH=gatk
HG=37
POS=gatk_intervals_b37.txt.gz
PADDING=1

##################
#   FUNCTIONS
##################
# Usage function 
function usage()
{
    echo -e "\nUsage:"
    echo -e "\t-i | --input-vcf\t\t- Variant call format input file (.vcf or .vcf.gz)."
    echo -e "\t-hg | --human-genome\t\t- Specifies the human genome reference version to use. Accepts either 37 or 38. Default: 37."
    echo -e "\t-p | --positions\t\t- Path to interval list of positions to include in gzipped format. The file contains columns snpID, rsID (if available) and position given as zero-indexed. Be sure to download the file matching your reference genome. "
    echo -e "\t-g | --gatk-path\t\t- Path to the gatk toolkit. Default: gatk."
    echo -e "\t-p | --padding\t\t- Number of bases to add padding around the included intervals. Default: 1."
    echo -e "\t-o | --out\t\t- Output path for the filtered file."
    echo -e "\t-h | --help\t\t- Print Usage message. "
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
        -hg | --human-genome)
            HG="$2"
            shift   # Remove argument name from processing 
            shift   # Remove argument value from processing 
            ;;
        -p | --positions)
            POS="$2"
            shift   # Remove argument name from processing 
            shift   # Remove argument value from processing 
            ;;
        -g | --gatk-path)
            GATKPATH="$2"
            shift   # Remove argument name from processing 
            shift   # Remove argument value from processing 
            ;;
        -p | --padding)
            PADDING="$2"
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


echo -e "Settings applied: \n\tinput-vcf: $INVCF \n\tHuman genome: GRCh${HG} \n\tPositions filter file: $POS \n\tGATK path: $GATKPATH \n\tOut path: $OUTVCF \n\n"

# GET THE SUBSET OF POSITIONS
gzip -dc $POS | cut -d" " -f3 > $POS.GRCh${HG}.list

# Get the contig lengths
if [[ "$INVCF" == *.gz ]]; then
    echo "this file is zipped"
    Rscript -e "con <- gzfile('$INVCF'); data1 <- readLines(con, warn = FALSE); flag <- TRUE; for (line in data1) { if (!startsWith(line, '##')) { break } else if (startsWith(line, '##contig=')) {cat(line, '\n') } }; close(con)" > contig_lengths.txt
else 
    echo "this file is not zipped"
    Rscript -e "con <- file('$INVCF'); data1 <- readLines(con, warn = FALSE); flag <- TRUE; for (line in data1) { if (!startsWith(line, '##')) { break } else if (startsWith(line, '##contig=')) {cat(line, '\n') } }; close(con)" > contig_lengths.txt
fi

# Check the contig lengths can contain all the positions needed 
Rscript -e "data1 <- readLines('contig_lengths.txt', warn = FALSE); data2 <- readLines('$POS.GRCh${HG}.list', warn = FALSE); dict1 <- list(); flag <- TRUE; for (entry in data1) { if (!startsWith(entry, '##')) { break } else if (startsWith(entry, '##contig=')) { key <- sub('.*ID=([0-9]+),.*', '\\1', entry); value <- as.integer(sub('.*length=(\\d+)>.*', '\\1', entry)); dict1[[key]] <- value } };for (entry in data2) { key <- sub('^(.*):.*', '\\1', entry); value <- as.integer(sub('.*-(\\d+)$', '\\1', entry)); if (!is.na(value) && key %in% names(dict1) && dict1[[key]] > value) { cat(entry, '\n') } };" > $POS.GRCh${HG}_filtered.list

# clean-up 
rm contig_lengths.txt
rm $POS.GRCh${HG}.list

# RUN VARIANT SELECTION
$GATKPATH SelectVariants --variant $INVCF --L $POS.GRCh${HG}_filtered.list --interval-padding $PADDING --output $OUTVCF

