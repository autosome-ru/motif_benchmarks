#!/usr/bin/env Rscript
source("/app/pcm2pfm_utils.R")

args = commandArgs(trailingOnly = TRUE)
input_filename = if (length(args) > 0) args[1] else "stdin"
output_filename = if (length(args) > 1) args[2] else "stdout"
pcm2pfm_files(input_filename, output_filename)
