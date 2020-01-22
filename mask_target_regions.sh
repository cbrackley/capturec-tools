#!/bin/bash

# Program to mask target regions in capC-MAP binned and smoothed profiles

if [[ $# != 3 ]]; then
    echo "Program to mask target regions in capC-MAP binned and smoothed profiles"
    echo "Usage: mask_target_regions.sh CONFIGFILE TARGETFILE DATADIR"
    exit 0
fi


configfile=$1
targetfile=$2
datadir=$3

if [[ ! -f $configfile ]]; then
    echo "ERROR : cannot find file $configfile"
    exit 0
fi

if [[ ! -f $targetfile ]]; then
    echo "ERROR : cannot find file $targetfile"
    exit 0
fi

if [[ ! -d $datadir ]]; then
    echo "ERROR : cannot find directory $datadir"
    exit 0
fi


##############################################################################


# get bin sizes
mapfile -t binwinds < <(awk '{
if ($1=="BIN") print $2,$3
}' $configfile )

# get exclusion zone (capC-MAP default is 500)
exc=$( awk 'BEGIN{E=-1}
{if ($1==EXCLUDE) {E=$2}}
END{if (E!=-1) {print E} else {print 500}}' $configfile )


# loop through bin sizes
for (( bw=0 ; bw<${#binwinds[@]} ; bw++ ))
do

    bin_window=(${binwinds[$bw]})
    bin=${bin_window[0]}
    window=${bin_window[1]}


    # get mask regions
    if [[ "$exc" -gt "$window" ]]; then
	maskwidth=$exc
    else
	maskwidth=$window
    fi    
    awk -v mw=$maskwidth '{OFS="\t"
        print $1,$2-mw,$3+mw,$4
         }' $targetfile > temp_mask.bed


    # loop through targets
    while read line
    do
	line=($line)
	target=${line[3]} 
	
	infile=${datadir}/captured_bin_${bin}_${window}_RPM_${target}.bdg
	outfile=${datadir}/masked_bin_${bin}_${window}_RPM_${target}.bdg

	if [[ -f $infile ]]; then
	    if [[ ! -f $outfile ]]; then
	    # now generate masked file
	    head -n 1 $infile > $outfile
	    sort -k1,1 -k2,2n \
		<( bedtools subtract -a $infile -b temp_mask.bed ) \
		<( bedtools intersect -wa -a $infile -b temp_mask.bed \
		   | awk '{OFS="\t"; print $1,$2,$3,-100}' ) \
		>> $outfile    
	    else
		echo "Warning, file $outfile already exists - did not overwrite"
	    fi
	else
	    echo "Warning, could not find file $infile so skipping this"
	fi

    done < $targetfile


done


rm temp_mask.bed

