#!/usr/bin/env bash
/app/print_matrix_from_json.R /workdir/config.json pcm > /motif.pcm
/app/print_matrix_from_json.R /workdir/config.json pfm > /motif.pfm
ln -s /data/common_data/selex_sequences.fa /seq.fa
ln -s /results /workdir/persistent
Rscript /app/calculate_roc.R --json --seed 13 > /workdir/persistent/result.json
