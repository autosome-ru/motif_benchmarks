#!/usr/bin/env Rscript
library('rjson')
library('optparse')
library('MASS')
library('caTools')
library('PRROC')
source('/app/utils.R')
source('/app/motif_preprocessing.R')
source('/app/peak_preprocessing.R')
source('/app/assembly_preprocessing.R')

plot_curve <- function(curve_data, image_filename, width = 800, height = 800) {
  png(image_filename, width = width, height = height)
  plot(curve_data, lwd = 4, cex.axis=1.4, cex.main=1.4, cex.lab=1.4, cex.sub=1.4)
  # txt <- paste("AUC = ", roc_data$auc, sept="")
  # legend(x="bottomright", legend=txt, col=c("red"), lty=1, lwd=4, cex=1.2)
  dummy <- dev.off()
}

roc_tpr_fpr <- function(roc_curve) {
  n_points <- length(roc_curve[,1])
  tpr <- list()
  fpr <- list()
  for(i in 1:n_points) {
    fpr[[i]] <- roc_curve[i, 1]
    tpr[[i]] <- roc_curve[i, 2]
  }
  return(list(tpr=unlist(tpr), fpr=unlist(fpr)))
}

pr_precision_recall <- function(pr_curve) {
  n_points <- length(pr_curve[,1])
  recall <- list()
  precision <- list()
  for(i in 1:n_points) {
    recall[[i]] <- pr_curve[i, 1]
    precision[[i]] <- pr_curve[i, 2]
  }
  return(list(precision=unlist(precision), recall=unlist(recall)))
}

roc_curve_as_points_list <- function(tpr, fpr) {
  n_bins <- length(tpr)
  roc_curve <- list()
  for(i in 1:n_bins) {
    roc_curve[[i]] <- list(tpr=tpr[[i]], fpr=fpr[[i]])
  }
  return(roc_curve)
}

pr_curve_as_points_list <- function(precision, recall) {
  n_bins <- length(precision)
  pr_curve <- list()
  for(i in 1:n_bins) {
    pr_curve[[i]] <- list(precision=precision[[i]], recall=recall[[i]])
  }
  return(pr_curve)
}

store_roc <- function(roc_data, output_filename) {
  write.table(list(fpr=roc_data$fpr, tpr=roc_data$tpr), row.names=FALSE, quote=FALSE, sep="\t", file=output_filename)
}

store_pr <- function(pr_data, output_filename) {
  write.table(list(recall=pr_data$recall, precision=pr_data$precision), row.names=FALSE, quote=FALSE, sep="\t", file=output_filename)
}

width = 800
height = 800

option_list = list(
  make_option(c("--assembly-name"), dest='assembly_name', action="store", type='character', default=NA, help="Choose assembly by name"),
  make_option(c("--motif-url"), dest= 'motif_url', type='character', default=NA, help="Use PFM file located at some URL"),
  make_option(c("--peaks-url"), dest= 'peaks_url', type='character', default=NA, help="Use peaks file located at some URL"),

  make_option(c("--pfm"), default=FALSE, action="store_true", help="Force use of PFM matrix"),
  make_option(c("--pcm"), default=FALSE, action="store_true", help="Force use of PCM matrix"),

  make_option(c("--narrowPeak"), dest='peak_format_narrowPeak', default=FALSE, action="store_true", help="Peaks are formatted in narrowPeak (peaks are reshaped into constant-size peaks around summit of a peak)"),
  make_option(c("--bed"), dest='peak_format_bed', default=FALSE, action="store_true", help="Peaks are formatted in bed (peaks are reshaped into constant-size peaks around center of a peak)"),

  make_option(c("--gz"), dest="compression_gz_peaks", default=FALSE, action="store_true", help="Force un-gzipping peaks"),
  make_option(c("--not-compressed"), dest="compression_no_peaks", default=FALSE, action="store_true", help="Prevent un-gzipping peaks"),

  make_option(c("--plot-roc"), dest="plot_roc_image", default=FALSE, action="store_true", help="Plot ROC curve"),
  make_option(c("--plot-roc-filename"), dest="roc_image_filename", type="character", default="/results/roc_curve.png", help="Specify plot filename [default=%default]"),
  make_option(c("--plot-pr"), dest="plot_pr_image", default=FALSE, action="store_true", help="Plot ROC curve"),
  make_option(c("--plot-pr-filename"), dest="pr_image_filename", type="character", default="/results/pr_curve.png", help="Specify plot filename [default=%default]"),
  make_option(c("--roc"), dest="store_roc", default=FALSE, action="store_true", help="Store ROC curve point"),
  make_option(c("--roc-filename"), dest="roc_filename",type="character", default="roc_curve.tsv", help="Specify ROC curve points filename [default=%default]"),
  make_option(c("--pr"), dest="store_pr", default=FALSE, action="store_true", help="Store PR curve point"),
  make_option(c("--pr-filename"), dest="pr_filename",type="character", default="pr_curve.tsv", help="Specify PR curve points filename [default=%default]"),
  make_option(c("--json"), dest="jsonify_results", default=FALSE, action="store_true", help="Print results as a json file"),

  make_option(c("--top"), type="integer", dest="num_top_peaks", default=500, help="Number of top peaks to take [default=%default]")
)
usage = paste("\n",
              "docker run --rm  -v {PFM}:/motif.pfm  -v {Peaks}:/peaks[.bed|.narrowPeak][.gz] -v {assembly folder}:/assembly/ pwmeval_chipseq --assembly-name {UCSC assembly name} [options]\n",
              "  or\n",
              "docker run --rm  -v {PFM}:/motif.pfm  -v {Peaks}:/peaks[.bed|.narrowPeak][.gz] -v {assembly folder}:/assembly/ -v {results}:/results  pwmeval_chipseq --plot --roc [options]\n",
              "  or\n",
              "docker run --rm  -v {assembly FASTA}:/assembly.fa -v {assembly chrom-sizes}:/assembly.chrom.sizes -v {assembly FASTA index}:/assembly.fa.fai  pwmeval_chipseq --motif-url {Motif URL} --peaks-url {Peaks URL} [options]\n",
              " or\n",
              "docker run --rm pwmeval_chipseq --motif-url {Motif URL} --peaks-url {Peaks URL} --assembly-name {UCSC assembly name} [options]\n")

description = paste("\n",
                    "Note!\n",
                    "  All local paths (for peaks file, PFM file, assembly folder and results folder) should be absolute.\n",
                    "  Peaks and motif format can be derived from extension.\n",
                    "  You can use bed/narrowPeak extensions for peaks files.\n",
                    "  Also you can use gz extension for gzipped peaks.\n",
                    "  So that /peaks.narrowPeak.gz is a correct way to pass a gzipped narrowPeak file.\n",
                    "  Options like --bed/--narrowPeak, --gz/--not-compressed, --pcm/--pfm override derived format,\n",
                    "  what is especially useful for passing data via url.\n",
                    "  In case when format is specified via options, `/peaks` and `/motif` with extension omitted can be used.\n",
                    "\n",
                    "  Assembly can be passed as separate files (FASTA, its index and chromosome sizes),\n",
                    "  or as a folder /assembly with necessary files and assembly name specified. If some of necessary files\n",
                    "  don't exists, they are calculated on the fly and are stored into /assembly folder.\n",
                    "  If the folder is mounted to a file system, an assembly and all supplementary files\n",
                    "  will be stored for reuse in the following runs (first time it'll be slow).\n",
                    "  --assembly-name [hg38/mm9/...] allows one to choose a necessary genome assembly amongst several ones.\n",
                    "  If specified assembly doesn't exist in a specified folder, it will be automatically downloaded from UCSC (not always possible).\n",
                    "  But one should note that if /assembly folder not mounted, genome and supplementary data will live only during container existence")
opt_parser <- OptionParser(option_list=option_list, usage = usage, description=description);
opts_and_args <- parse_args(opt_parser, positional_arguments=TRUE);
opts <- opts_and_args[[1]]
args <- opts_and_args[[2]]

dummy = obtain_and_preprocess_motif(opts)
dummy = obtain_and_preprocess_peaks(opts)
assembly = obtain_and_preprocess_assembly(opts)

system(paste("sort -k5,5nr /workdir/peak_centers_scored.bed | head -n", opts$num_top_peaks, " > /workdir/top_peaks.bed"))
system(paste("/app/bedtools slop -i /workdir/top_peaks.bed -g ", shQuote(assembly$sizes_fn), " -l 124 -r 125  > /workdir/positive_peaks.bed"))
system(paste("/app/bedtools slop -i /workdir/top_peaks.bed -g ", shQuote(assembly$sizes_fn), " -l -301 -r 550  > /workdir/negative_peaks.bed"))
system(paste("/app/bedtools getfasta -bed /workdir/positive_peaks.bed -fi ", shQuote(assembly$fasta_fn), "  > /workdir/positive.seq"))
system(paste("/app/bedtools getfasta -bed /workdir/negative_peaks.bed -fi ", shQuote(assembly$fasta_fn), "  > /workdir/negative.seq"))

system("/app/pwm_scoring -r -u -m /workdir/motif.pfm /workdir/positive.seq  > /workdir/positive_PWM.out")
system("/app/pwm_scoring -r -u -m /workdir/motif.pfm /workdir/negative.seq  > /workdir/negative_PWM.out")

pos <- as.matrix(read.table("/workdir/positive_PWM.out"))
neg <- as.matrix(read.table("/workdir/negative_PWM.out"))
roc_infos <- roc.curve(pos, neg, curve=TRUE)
pr_infos <- pr.curve(pos, neg, curve=TRUE)
roc_auc = roc_infos$auc
pr_auc = pr_infos$auc.integral
pr_auc_dg = pr_infos$auc.davis.goadrich
roc_data <- roc_tpr_fpr(roc_infos$curve)
pr_data <- pr_precision_recall(pr_infos$curve)

if (opts$store_roc) {
  store_roc(roc_data, file.path('/results', opts$roc_filename))
}
if (opts$store_pr) {
  store_pr(pr_data, file.path('/results', opts$pr_filename))
}

if (opts$plot_roc_image) {
  plot_curve(roc_infos, opts$roc_image_filename, width = width, height = height)
}
if (opts$plot_pr_image) {
  plot_curve(pr_infos, opts$pr_image_filename, width = width, height = height)
}

if (opts$jsonify_results) {
  metrics <- list(roc_auc=roc_auc, pr_auc=pr_auc, pr_auc_davis_goadrich=pr_auc_dg)
  supplementary <- list(
    roc_curve=roc_curve_as_points_list(roc_data$tpr, roc_data$fpr),
    pr_curve=pr_curve_as_points_list(pr_data$precision, pr_data$recall)
  )
  results <- list(metrics=metrics, supplementary=supplementary)
  writeLines(toJSON(results))
} else{
  writeLines(as.character(roc_auc))
  writeLines(as.character(pr_auc))
}
