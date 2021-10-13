arglist_sequence_options = list(
  make_option(c("--seq"), dest='seq_fn', type='character', default=NA, help="Specify FASTA file"),
  make_option(c("--seq-url"), dest='seq_url', type='character', default=NA, help="Use FASTA file located at some URL"),

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
  make_option(c("--maxnum-reads"), dest="maxnum_reads", type="integer", default=NA, help="Set a maximal number of reads to subsample")
)

arglist_motif_options = list(
  make_option(c("--motif"), dest='motif_fn', type='character', default=NA, help="Specify motif file"),
  make_option(c("--motif-url"), dest='motif_url', type='character', default=NA, help="Use PFM file located at some URL"),
  make_option(c("--pfm"), default=FALSE, action="store_true", help="Force use of PFM matrix"),
  make_option(c("--pcm"), default=FALSE, action="store_true", help="Force use of PCM matrix"),
  make_option(c("--pseudo-weight"), dest="pseudo_weight", type="double", default=0.0001, help="Set a pseudo-weight to re-normalize the frequencies of the positional-probability matrix (PFM) [default=%default]")
)
