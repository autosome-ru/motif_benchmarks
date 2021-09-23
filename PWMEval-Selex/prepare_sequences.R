#!/usr/bin/env Rscript
library(optparse)
source('/app/utils.R')
source('/app/seq_preprocessing.R')
source('/app/arglist_options.R')

option_list = arglist_sequence_options
usage = paste("\n",
              "docker run --rm --entrypoint /app/prepare_sequences.R -v {results}:/sequences -v {Selex FASTA}:/seq[.fa|.fq][.gz]  pwmeval_selex [options]\n",
              " or\n",
              "docker run --rm --entrypoint /app/prepare_sequences.R -v {results}:/sequences  pwmeval_selex --seq-url {Selex FASTA URL} [options]\n")
description = paste("\n",
                    "Note!\n",
                    "  All local paths (for FASTA file and results folder) should be absolute.\n",
                    "  Sequences format can be derived from extension.\n",
                    "  You can use fa/fasta extensions for FASTA files and fq/fastq for FASTQ files.\n",
                    "  Also you can use gz extension for gzipped sequences.\n",
                    "  So that /seq.fastq.gz is a correct way to pass a gzipped FASTQ file.\n",
                    "  Options like --fa/--fq, --gz/--not-compressed override derived format,\n",
                    "  what is especially useful for passing data via url.\n",
                    "  In case when format is specified via options, `/seq` with extension omitted can be used.\n")
opt_parser <- OptionParser(option_list=option_list, usage = usage, description=description);
opts_and_args <- parse_args(opt_parser, positional_arguments=TRUE);
opts <- opts_and_args[[1]]
args <- opts_and_args[[2]]

pos_seq_fn = obtain_and_preprocess_sequences(opts)
neg_seq_fn = tempfile()

if (is.na(opts$seed)) {
  system(paste("/app/seqshuffle", shQuote(pos_seq_fn), ">", shQuote(neg_seq_fn)))
} else {
  system(paste("/app/seqshuffle -s", opts$seed, shQuote(pos_seq_fn), ">", shQuote(neg_seq_fn)))
}

pos_seq_fn = append_flanks(pos_seq_fn, opts)
neg_seq_fn = append_flanks(neg_seq_fn, opts)

if (endsWith(opts$positive_fn, '.gz')) {
  pos_seq_fn = compress_file(pos_seq_fn, "gz")
}
if (endsWith(opts$negative_fn, '.gz')) {
  neg_seq_fn = compress_file(neg_seq_fn, "gz")
}

dummy <- file.copy(pos_seq_fn, opts$positive_fn)
dummy <- file.copy(neg_seq_fn, opts$negative_fn)
