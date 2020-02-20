#!/usr/bin/env sh
awk -e '{print $1 "\t" int(($2+$3)/2) "\t" int(1 + ($2+$3)/2) "\t" "." "\t" $5 }' -- "$1"
