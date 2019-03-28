#!/usr/bin/env Rscript
library(rjson)
library(optparse)

# Define Approximate AUC computating function
auc_approx <- function(pos, neg, top_fraction, n_bins) {
  N = round(top_fraction*length(POS))
  pos_sorted = sort(pos, decreasing=T)[1:N]
  neg_sorted = sort(neg, decreasing=T)[1:N]

  # Define common break points for score historgrams (n_bins, the more the better)
  max = max(pos_sorted, neg_sorted)
  min = min(pos_sorted, neg_sorted)
  inc = (max - min) / n_bins
  breaks = seq(min, max, inc)

  # Compute density distributions for scores of positives and negative examples
  # Note: *(max-min)/n_bins) ensures that densities sum to 1!
  dpos = rev(hist(pos_sorted, breaks=breaks, plot=F)$density * (max-min)/n_bins)
  dneg = rev(hist(neg_sorted, breaks=breaks, plot=F)$density * (max-min)/n_bins)

  # AUROC computation for binned score values (fast)
  sens = rep(0, n_bins)
  fpr  = rep(0, n_bins)
  sens[1] = dpos[1]
  fpr[1]  = dneg[1]
  for(i in 2:n_bins) {
    sens[i] = sens[i-1] + dpos[i]
    fpr[i] = fpr[i-1] + dneg[i]
  }
  AUC = dneg[1] * 0.5 * (sens[1] - 0)
  for(i in 2:n_bins) {
    AUC = AUC + dneg[i] * 0.5 * (sens[i] + sens[i-1])
  }
  return( list(auc=AUC, fpr=fpr, tpr=sens) )
}

store_roc <- function(roc_data, output_filename) {
  write.table(list(tpr=roc_data$tpr, fpr=roc_data$fpr), row.names=FALSE, quote=FALSE, sep="\t", file=output_filename)
}

plot_roc <- function(roc_data, output_filename) {
  png(output_filename, width = width, height = height, pointsize = pointsize)
  plot(1 - roc_data$fpr, roc_data$tpr, t="s", col="red", xlab="1-Specificity", ylab="Sensitivity", main="ROC curve of TF motif predictor", 
      xlim=c(1,0), lwd = 4, cex.axis=1.4, cex.main=1.4, cex.lab=1.4, cex.sub=1.4, family='Liberation')
  txt <- paste("AUC = ", signif(roc_data$auc, digits = 6), sept="")
  legend(x="bottomright", legend=txt, col=c("red"), lty=1, lwd=4, cex=1.2)
  dummy <- dev.off() # store to a variable to suppress a message
}

roc_curve_as_points_list <- function(tpr, fpr) {
  n_bins <- length(tpr)
  roc_curve <- list()
  for(i in 1:n_bins) {
    roc_curve[[i]] <- list(tpr=tpr[[i]], fpr=fpr[[i]])
  }
  return(roc_curve)
}

guess_seq_format <- function(filename) {
  if (endsWith(filename, ".gz")) {
    compression = "gz"
  } else {
    compression = "no"
  }

  filename = sub("\\.gz$", "", filename)

  if (endsWith(filename, ".fa") || endsWith(filename, ".fasta")) {
    seq_format = "fasta"
  } else if (endsWith(filename, ".fq") || endsWith(filename, ".fastq")) {
    seq_format = "fastq"
  } else { # default
    seq_format = "fasta"
  }
  return(list(compression=compression, seq_format=seq_format))
}

# override sequences format/compression based on command-line options
refine_seq_format_guess <- function(guessed_format, opts) {
  seq_format = guessed_format$seq_format
  compression = guessed_format$compression

  if (opts$seq_format_fasta) {
    seq_format = 'fasta'
  }
  if (opts$seq_format_fastq) {
    seq_format = 'fastq'
  }
  if (opts$compression_no) {
    compression = 'no'
  }
  if (opts$compression_gz) {
    compression = 'gz'
  }
  return(list(compression=compression, seq_format=seq_format))
}

download_file <- function(url) {
  # We want to find out original name of the downloaded file.
  # So we create a folder and download the only file to it
  # and can get the name of that file
  dirname = tempfile()
  dir.create(dirname)
  system(paste("wget -P", shQuote(dirname), shQuote(url)))
  original_fn = list.files(dirname)[1]
  return(file.path(dirname, original_fn))
}

decompress_file <- function(filename, compression) {
  if (compression == "no") {
    return(filename)
  } else if (compression == "gz") {
    tmp_fn = tempfile()
    system(paste("gzip -cd", filename, " > ", shQuote(tmp_fn)))
    return(tmp_fn)
  } else {
    simpleError("Unknown compression format")
  }
}

fastq2fasta <- function(seq_filename) {
  tmp_fn = tempfile()
  system(paste("/app/fastq2fasta.sh", seq_filename, " > ", shQuote(tmp_fn)))
  return(tmp_fn)
}

convert2fasta <- function(seq_filename, seq_format) {
  if (seq_format == 'fasta') {
    return(seq_filename)
  } else if (seq_format == 'fastq') {
    return(fastq2fasta(seq_filename))
  } else {
    simpleError("Incorrect format")
  }
}

# in addition to filtering it also joins FASTA spreaded on multiple lines into single-string format
filter_fasta <- function(seq_filename, opts) {
  only_acgt = 'yes'
  if (opts$allow_iupac) {
    only_acgt = 'no'
  }

  seq_length = 'no'
  if (!is.na(opts$seq_length)) {
    seq_length = opts$seq_length
  }

  tmp_fn = tempfile()
  system(paste("/app/filter_fasta", shQuote(seq_filename), seq_length,  only_acgt, " > ", shQuote(tmp_fn)))
  return(tmp_fn)
}

filter_redundant <- function(seq_filename, opts) {
  if (opts$non_redundant) {
    tmp_fn = tempfile()
    system(paste("/app/filter_redundant.sh", shQuote(seq_filename), " > ", shQuote(tmp_fn)))
    return(tmp_fn)
  } else {
    return(seq_filename)
  }
}

find_mounted_sequences_file <- function() {
  fasta_files = c('/seq.fasta', '/seq.fa', '/seq.fasta.gz', '/seq.fa.gz')
  fastq_files = c('/seq.fastq', '/seq.fq', '/seq.fastq.gz', '/seq.fq.gz')
  no_format_files = c('/seq', '/seq.gz')
  acceptable_seq_files = c(no_format_files, fasta_files, fastq_files)
  existing_seq_files = file.exists(acceptable_seq_files)

  if (sum(existing_seq_files) == 0) {
    simpleError("Provide a file with SELEX sequences. Either mount to /seq or its counterparts, or pass it via URL.")
  } else if (sum(existing_seq_files) > 1) {
    simpleError("Provide the only file with SELEX sequences.")
  }

  seq_filename = acceptable_seq_files[existing_seq_files][1]
  return(seq_filename)
}


obtain_and_preprocess_sequences <- function(opts) {
  if (!is.na(opts$seq_url)) {
    seq_filename = download_file(opts$seq_url)
  } else {
    seq_filename = find_mounted_sequences_file()
  }

  seq_format_info = refine_seq_format_guess(guess_seq_format(seq_filename), opts)

  # process sequences file into uncompressed FASTA file
  seq_filename = decompress_file(seq_filename, seq_format_info$compression)
  seq_filename = convert2fasta(seq_filename, seq_format_info$seq_format)
  seq_filename = filter_fasta(seq_filename, opts)
  seq_filename = filter_redundant(seq_filename, opts)
  file.copy(seq_filename, "/workdir/positive.fa")
}

obtain_and_preprocess_motif <- function(opts) {
  if (!is.na(opts$motif_url)) {
    motif_filename = download_file(opts$motif_url)
  } else {
    motif_filename = find_mounted_motif_file()
  }

  motif_format = refine_motif_format_guess(guess_motif_format(motif_filename), opts)
  motif_filename = get_ppm(motif_filename, motif_format)
  file.copy(motif_filename, "/workdir/motif.ppm")
}

find_mounted_motif_file <- function() {
  ppm_files = c('/motif.ppm', '/motif.pfm', '/matrix.ppm', '/matrix.pfm')
  pcm_files = c('/motif.pcm', '/matrix.pcm')
  no_format_files = c('/motif', '/matrix')
  acceptable_motif_files = c(no_format_files, ppm_files, pcm_files)
  existing_motif_files = file.exists(acceptable_motif_files)

  if (sum(existing_motif_files) == 0) {
    simpleError("Provide a file with positional frequencies/counts matrix. Either mount to /motif or its counterparts, or pass it via URL.")
  } else if (sum(existing_motif_files) > 1) {
    simpleError("Provide the only file with positional frequencies/counts matrix")
  }

  motif_filename = acceptable_motif_files[existing_motif_files][1]
  return(motif_filename)
}

guess_motif_format <- function(filename) {
  if (endsWith(filename, ".pcm")) {
    seq_format = "pcm"
  } else if (endsWith(filename, ".pfm") || endsWith(filename, ".ppm")) {
    seq_format = "ppm"
  } else { # default
    seq_format = "pcm"
  }
  return(seq_format)
}

# override motif format based on command-line options
refine_motif_format_guess <- function(guessed_format, opts) {
  format = guessed_format
  if (opts$ppm) {
    format = 'ppm'
  } else if (opts$pcm) {
    format = 'pcm'
  }
  return(format)
}

get_ppm <- function(filename, format) {
  if (format == 'pcm') {
    tmp_fn = tempfile()
    system(paste("/app/pcm2ppm.R", shQuote(filename), " > ", shQuote(tmp_fn)))
    return(tmp_fn)
  } else if (format == 'ppm') {
    return(filename)
  } else {
    simpleError("Unknown motif format")
  }
}


width = 800
height = 800
quality = 100
pointsize = 20

option_list = list(
  make_option(c("--seq-url"), dest= 'seq_url', type='character', default=NA, help="Use FASTA file located at some URL"),
  make_option(c("--motif-url"), dest= 'motif_url', type='character', default=NA, help="Use PPM file located at some URL"),
  make_option(c("--plot"), dest="plot_image", default=FALSE, action="store_true", help="Plot ROC curve"),
  make_option(c("--plot-filename"), dest="image_filename", type="character", default="roc_curve.png", metavar='FILENAME', help="Specify plot filename [default=%default]"),
  make_option(c("--roc"), dest="store_roc", default=FALSE, action="store_true", help="Store ROC curve point"),
  make_option(c("--roc-filename"), dest="roc_filename",type="character", default="roc_curve.tsv", help="Specify ROC curve points filename [default=%default]"),
  make_option(c("--json"), dest="jsonify_results", default=FALSE, action="store_true", help="Print results as a json file"),

  make_option(c("--gz"), dest="compression_gz", default=FALSE, action="store_true", help="Force un-gzipping sequences"),
  make_option(c("--not-compressed"), dest="compression_no", default=FALSE, action="store_true", help="Prevent un-gzipping sequences"),

  make_option(c("--fastq"), dest='seq_format_fastq', default=FALSE, action="store_true", help="Use FASTQ"),
  make_option(c("--fasta"), dest='seq_format_fasta', default=FALSE, action="store_true", help="Use FASTA"),
  
  make_option(c("--ppm"), default=FALSE, action="store_true", help="Force use of PPM matrix"),
  make_option(c("--pcm"), default=FALSE, action="store_true", help="Force use of PCM matrix"),

  make_option(c("--seq-length"), dest="seq_length", type='integer', default=NA, action="store", metavar="LENGTH", help="Specify length of sequences. All sequences of different length will be rejected."),
  make_option(c("--allow-iupac"), dest="allow_iupac", default=FALSE, action="store_true", help="Allow IUPAC sequences (by default only ACGT are valid)."),
  make_option(c("--non-redundant"), dest="non_redundant", default=FALSE, action="store_true", help="Retain only unique sequences."),
  make_option(c("--top"), dest="top_fraction", type="double", default=0.1, help="Fraction of top sequences to take [default=%default]"),
  make_option(c("--bins"), dest="num_bins", type="integer", default=1000, help="Number of bins for ROC computations [default=%default]"),
  make_option(c("--pseudo-weight"), dest="pseudo_weight", type="double", default=0.0001, help="Set a pseudo-weight to re-normalize the frequencies of the positional-probability matrix (PPM) [default=%default]")
)
usage = paste("\n",
              "docker run --rm  -v {PPM}:/motif.ppm  -v {Selex FASTA}:/seq[.fa|.fq][.gz]  pwmeval_selex [options]\n",
              "  or\n",
              "docker run --rm  -v {PPM}:/motif.ppm  -v {Selex FASTA}:/seq[.fa|.fq][.gz]  -v {results}:/results  pwmeval_selex [options]\n",
              " or\n",
              "docker run --rm  pwmeval_selex --motif-url {Motif URL} --seq-url {Selex FASTA URL} [options]\n")
description = paste("\n",
                    "Note!\n",
                    "  All local paths (for FASTA file, PPM file and results folder) should be absolute.\n",
                    "  Sequences format can be derived from extension.\n",
                    "  You can use fa/fasta extensions for FASTA files and fq/fastq for FASTQ files.\n",
                    "  Also you can use gz extension for gzipped sequences.\n",
                    "  So that /seq.fastq.gz is a correct way to pass a gzipped FASTQ file.\n",
                    "  Options like --fasta/--fastq, --gz/--not-compressed override derived format,\n",
                    "  what is especially useful for passing sequences via url.\n",
                    "  In case when format is specified via options `/seq` with extension omitted can be used.\n")
opt_parser <- OptionParser(option_list=option_list, usage = usage, description=description);
opts_and_args <- parse_args(opt_parser, positional_arguments=TRUE);
opts <- opts_and_args[[1]]
args <- opts_and_args[[2]]

dummy = obtain_and_preprocess_sequences(opts)
dummy = obtain_and_preprocess_motif(opts)

system(paste("/app/seqshuffle /workdir/positive.fa > /workdir/negative.fa"))
system(paste("/app/pwm_scoring -r -w", opts$pseudo_weight, "-m motif.ppm /workdir/positive.fa  > /workdir/PPM_scores_positive.txt"))
system(paste("/app/pwm_scoring -r -w", opts$pseudo_weight, "-m motif.ppm /workdir/negative.fa  > /workdir/PPM_scores_negative.txt"))

POS <- log10(as.numeric(read.table("/workdir/PPM_scores_positive.txt", header=F)[,1]))
NEG <- log10(as.numeric(read.table("/workdir/PPM_scores_negative.txt", header=F)[,1]))

roc_data <- auc_approx(POS, NEG, opts$top_fraction, opts$num_bins)
if (opts$plot_image) {
  plot_roc(roc_data, file.path('/results', opts$image_filename))
}
if (opts$store_roc) {
  store_roc(roc_data, file.path('/results', opts$roc_filename))
}

if (opts$jsonify_results) {
  metrics <- list(roc_auc=roc_data$auc)
  supplementary <- list(roc_curve=roc_curve_as_points_list(roc_data$tpr, roc_data$fpr))
  results <- list(metrics=metrics, supplementary=supplementary)
  writeLines(toJSON(results))
} else{
  writeLines(as.character(roc_data$auc))
}
