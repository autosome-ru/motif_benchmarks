#!/usr/bin/env Rscript
library(optparse)
source('/app/utils.R')
source('/app/seq_preprocessing.R')
source('/app/arglist_options.R')

option_list = c(
  arglist_sequence_options,
  make_option(c("--positive-file"), dest='positive_fn', type='character', default='/sequences/positive.fa.gz', help="Resulting positive sequences filename"),
  make_option(c("--negative-file"), dest='negative_fn', type='character', default='/sequences/negative.fa.gz', help="Resulting negative sequences filename")
)
usage = paste(
  "\n",
  "docker run --rm \\",
  "           --entrypoint /app/prepare_sequences.R \\",
  "           --volume {Selex FASTA/FASTQ}:/seq[.fa|.fq][.gz] \\",
  "           --volume {results}:/sequences \\",
  "           pwmeval_selex \\",
  "             --seq /seq.fa.gz \\",
  "             [options]\n"
)
description = paste("\n",
                    "Note!\n",
                    "  All local paths (for FASTA file and results folder) should be absolute.\n",
                    "  Sequences format can be derived from extension.\n",
                    "  You can use fa/fasta extensions for FASTA files and fq/fastq for FASTQ files.\n",
                    "  Also you can use gz extension for gzipped sequences.\n",
                    "  So that `--seq /seq.fastq.gz` is a correct way to pass a gzipped FASTQ file.\n",
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

dummy <- file.copy(pos_seq_fn, opts$positive_fn, overwrite=TRUE)
dummy <- file.copy(neg_seq_fn, opts$negative_fn, overwrite=TRUE)
