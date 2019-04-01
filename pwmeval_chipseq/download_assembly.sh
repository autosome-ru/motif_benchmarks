#!/usr/bin/env sh
set -e -u -o pipefail
ASSEMBLY_NAME=$1
OUTPUT_FILE=$2
mkdir -p "/tmp/assembly_${ASSEMBLY_NAME}/"
rsync -a -P "rsync://hgdownload.soe.ucsc.edu/goldenPath/${ASSEMBLY_NAME}/chromosomes/*.fa.gz"  "/tmp/assembly_${ASSEMBLY_NAME}/"
find "/tmp/assembly_${ASSEMBLY_NAME}/" -type f -name '*.fa.gz' -print0 | xargs -0 zcat > "/tmp/assembly_${ASSEMBLY_NAME}.fa"
mv "/tmp/assembly_${ASSEMBLY_NAME}.fa" "$OUTPUT_FILE"
