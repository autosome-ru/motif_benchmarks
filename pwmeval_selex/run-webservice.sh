#!/usr/bin/env bash
/app/print_matrix_from_json.R /workdir/config.json pcm > /motif.pcm
/app/print_matrix_from_json.R /workdir/config.json pfm > /motif.pfm
Rscript /app/calculate_roc.R --json
