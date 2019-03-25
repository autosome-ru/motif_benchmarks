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

width = 800
height = 800
quality = 100
pointsize = 20

option_list = list(
  make_option(c("--url"), type='character', default=NA, help="Use FASTA file located at some URL"),
  make_option(c("--plot"), dest="plot_image", default=FALSE, action="store_true", help="Plot ROC curve"),
  make_option(c("--plot-filename"), dest="image_filename", type="character", default="roc_curve.png", metavar='FILENAME', help="Specify plot filename [default=%default]"),
  make_option(c("--roc"), dest="store_roc", default=FALSE, action="store_true", help="Store ROC curve point"),
  make_option(c("--roc-filename"), dest="roc_filename",type="character", default="roc_curve.tsv", help="Specify ROC curve points filename [default=%default]"),
  make_option(c("--json"), dest="jsonify_results", default=FALSE, action="store_true", help="Print results as a json file"),
  make_option(c("--fastq"), default=FALSE, action="store_true", help="Use FASTQ file instead of FASTA"),
  make_option(c("--top"), dest="top_fraction", type="double", default=0.1, help="Fraction of top sequences to take [default=%default]"),
  make_option(c("--bins"), dest="num_bins", type="integer", default=1000, help="Number of bins for ROC computations [default=%default]"),
  make_option(c("--pseudo-weight"), dest="pseudo_weight", type="double", default=0.0001, help="Set a pseudo-weight to re-normalize the frequencies of the letter-probability matrix (LPM) [default=%default]")
)
usage = paste("\n",
              "docker run --rm  -v {PPM}:/motif.ppm  -v {Selex FASTA}:/seq.fa  pwmeval_selex [options]\n",
              "  or\n",
              "docker run --rm  -v {PPM}:/motif.ppm  -v {Selex FASTA}:/seq.fa  -v {results}:/results  pwmeval_selex [options]\n",
              " or\n",
              "docker run --rm  -v {PPM}:/motif.ppm  pwmeval_selex --url {Selex FASTA URL} [options]\n")
description = paste("\n",
                    "Note!\n",
                    "  All local paths (for FASTA file, PPM file and results folder) should be absolute.\n",
                    "  FASTA and FASTQ are interchangeable. Gzipped files with .gz extension can be used")
opt_parser <- OptionParser(option_list=option_list, usage = usage, description=description);
opts_and_args <- parse_args(opt_parser, positional_arguments=TRUE);
opts <- opts_and_args[[1]]
args <- opts_and_args[[2]]

positive_sequences_fn = args[1]
matrix_fn = args[2]

if (is.na(opts$url)) {
  system(paste("cp /seq.fa /workdir/positive.seq"))
} else {
  url <- opts$url
  if (endsWith(url, '.gz')) {
    system(paste("wget -O /workdir/positive.seq.gz", shQuote(url)))
    system(paste("gzip -d /workdir/positive.seq.gz"))
  } else {
    system(paste("wget -O /workdir/positive.seq", shQuote(url)))
  }
}

if (opts$fastq) {
  system("/app/fastq2fasta.sh /workdir/positive.seq  > /workdir/positive.fa")
} else {
  system("cp /workdir/positive.seq /workdir/positive.fa")
}

system(paste("ln -s /motif.ppm /workdir/motif.ppm"))
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
