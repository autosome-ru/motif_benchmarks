#!/usr/bin/env sh
# We suppose that FASTA is already in single-line format
FLANK_5=$2
FLANK_3=$3
awk "{if ( NR%2==1 ){print \$0} else { print \"${FLANK_5}\" \$0 \"${FLANK_3}\" }}" $1
