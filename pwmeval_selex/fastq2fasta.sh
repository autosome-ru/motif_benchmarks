#!/usr/bin/env sh
awk '{if ( NR%4==1 || NR%4==2){print $0}}' $1 | sed 's/^@/>/'
