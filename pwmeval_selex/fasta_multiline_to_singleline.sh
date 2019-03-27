#!/usr/bin/env sh
awk -e '($0 ~ /^>/){ if(NR!=1){ printf("\n") }; print $0 }; ($0 !~ /^>/) { printf("%s",$0) }' $1
