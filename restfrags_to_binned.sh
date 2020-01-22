#!/bin/bash

# Program to generate files for restriction enzyme density binned and smoothed
# at different levels

if [[ $# != 6 ]]; then
    echo "Program to generate a restriction enzyme density profile with binning and smoothing"
    echo "Usage: restfrags_to_binned.sh FRAGFILE CHROMSIZES CHROM BIN WINDOW OUTFILE"
    echo "    where FRAGFILE       is the restriction enzyme list used by capC-MAP"
    echo "    where CHROMSIZES     is the chromosome sizes file used by capC-MAP"
    echo "    where CHROM          is the chromosome of interested, e.g. chrZ"
    echo "    where BIN            is the bin size"
    echo "    where WINDOW         is the smoothing window size"
    echo "    where OUTFILE        is a file name for the output bedGraph"
    exit 0
fi


fragfile=$1 
chromsizes=$2 
chrom=$3 
bin=$4 
window=$5
outfile=$6 

if [[ ! -f $fragfile ]]; then
    echo "ERROR : cannot find file $fragfile"
    exit 0
fi

if [[ ! -f $chromsizes ]]; then
    echo "ERROR : cannot find file $chromsizes"
    exit 0
fi

if [[ -f $outfile ]]; then
    echo "ERROR : file $outfile already exists, will not overwrite"
    exit 0
fi

##########################################################################################

awk -v c=$chrom '{OFS="\t"
if ($1==c) print $1,$2,$2+2,1
}' $fragfile > temp

capCpileup2binned -i temp -c $chromsizes -o $outfile -b $bin $window -t restfrags


rm temp
