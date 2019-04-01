#!/usr/bin/env sh
set -e -u -o pipefail
ASSEMBLY_NAME=$1
mkdir -p /tmp/assembly_${ASSEMBLY_NAME}/
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/${ASSEMBLY_NAME}/chromosomes/  /tmp/assembly_${ASSEMBLY_NAME}/
zcat /tmp/assembly_${ASSEMBLY_NAME}/*.fa.gz
