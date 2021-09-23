#!/usr/bin/env Rscript
library(optparse)
source('/app/utils.R')
source('/app/seq_preprocessing.R')

option_list = list(
  make_option(c("--seq-url"), dest= 'seq_url', type='character', default=NA, help="Use FASTA file located at some URL"),

  make_option(c("--gz"), dest="compression_gz", default=FALSE, action="store_true", help="Force un-gzipping sequences"),
  make_option(c("--not-compressed"), dest="compression_no", default=FALSE, action="store_true", help="Prevent un-gzipping sequences"),

  make_option(c("--fastq"), dest='seq_format_fastq', default=FALSE, action="store_true", help="Use FASTQ"),
  make_option(c("--fasta"), dest='seq_format_fasta', default=FALSE, action="store_true", help="Use FASTA"),
  
  make_option(c("--seq-length"), dest="seq_length", type='integer', default=NA, action="store", metavar="LENGTH", help="Specify length of sequences. All sequences of different length will be rejected."),
  make_option(c("--allow-iupac"), dest="allow_iupac", default=FALSE, action="store_true", help="Allow IUPAC sequences (by default only ACGT are valid)."),
  make_option(c("--non-redundant"), dest="non_redundant", default=FALSE, action="store_true", help="Retain only unique sequences."),
  make_option(c("--flank-5"), dest="flank_5", type='character', default='', help="Append 5'-flanking sequence (adapter+barcode) to sequences"),
  make_option(c("--flank-3"), dest="flank_3", type='character', default='', help="Append 3'-flanking sequence (adapter+barcode) to sequences"),

  make_option(c("--seed"), type="integer", default=NA, help="Set a seed for generation of random negative control"),
  make_option(c("--maxnum-reads"), dest="maxnum_reads", type="integer", default=NA, help="Set a maximal number of reads to subsample"),
  
  make_option(c("--positive-file"), dest='positive_fn', type='character', default='/sequences/positive.fa.gz', help="Resulting positive sequences filename"),
  make_option(c("--negative-file"), dest='negative_fn', type='character', default='/sequences/negative.fa.gz', help="Resulting negative sequences filename")
)
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
