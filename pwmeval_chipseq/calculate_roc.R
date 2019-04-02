#!/usr/bin/env Rscript
library('rjson')
library('optparse')
library('MASS')
library('caTools')
suppressPackageStartupMessages(library('pROC'))
source('/app/utils.R')
source('/app/motif_preprocessing.R')
source('/app/peak_preprocessing.R')
source('/app/assembly_preprocessing.R')

plot_roc <- function(roc_data, image_filename, width = 800, height = 800, pointsize = 20) {
  png(image_filename, width = width, height = height, pointsize = pointsize)
  plot(roc_data, main="ROC curve of motif predictor", col = "red", lwd = 4, cex.axis=1.4, cex.main=1.4, cex.lab=1.4, cex.sub=1.4)
  txt <- paste("AUC = ", roc_data$auc, sept="")
  legend(x="bottomright", legend=txt, col=c("red"), lty=1, lwd=4, cex=1.2)
  dummy <- dev.off()
}

roc_tpr_fpr <- function(roc_infos) {
  roc_coords <- coords(roc_infos, 'all', as.list=TRUE, ret=c("sensitivity", "1-specificity"))
  n_points <- length(roc_coords)
  tpr <- list()
  fpr <- list()
  for(i in 1:n_points) {
    tpr[[i]] <- roc_coords[[i]]$sensitivity
    fpr[[i]] <- roc_coords[[i]]$"1-specificity"
  }
  return(list(tpr=unlist(tpr), fpr=unlist(fpr)))
}

roc_curve_as_points_list <- function(tpr, fpr) {
  n_bins <- length(tpr)
  roc_curve <- list()
  for(i in 1:n_bins) {
    roc_curve[[i]] <- list(tpr=tpr[[i]], fpr=fpr[[i]])
  }
  return(roc_curve)
}

store_roc <- function(roc_data, output_filename) {
  print(names(roc_data))
  print(class(roc_data$tpr))
  print(class(roc_data$tpr))
  print(class(roc_data$fpr))
  write.table(list(tpr=roc_data$tpr, fpr=roc_data$fpr), row.names=FALSE, quote=FALSE, sep="\t", file=output_filename)
}

width = 800
height = 800
pointsize = 20

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

  make_option(c("--plot"), dest="plot_image", default=FALSE, action="store_true", help="Plot ROC curve"),
  make_option(c("--plot-filename"), dest="image_filename", type="character", default="/results/roc_curve.png", help="Specify plot filename [default=%default]"),
  make_option(c("--roc"), dest="store_roc", default=FALSE, action="store_true", help="Store ROC curve point"),
  make_option(c("--roc-filename"), dest="roc_filename",type="character", default="roc_curve.tsv", help="Specify ROC curve points filename [default=%default]"),
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
W = wilcox.test(pos, neg, alternative ="g")$statistic
AUC = W/length(pos)/length(neg)
pos_labs <- rep(1, 500)
neg_labs <- rep(0, 500)
pos <- cbind(pos, pos_labs)
neg <- cbind(neg, neg_labs)
comb_sets <- rbind(pos, neg)
roc_infos <- roc(response = comb_sets[, 2], predictor = comb_sets[, 1])
auc = as.numeric(roc_infos$auc) # by default auc is a complex object which is displayed not as a number but number + text
roc_data <- roc_tpr_fpr(roc_infos)

if (opts$store_roc) {
  store_roc(roc_data, file.path('/results', opts$roc_filename))
}

if (opts$plot_image) {
  plot_roc(roc_infos, opts$image_filename, width = width, height = height, pointsize = pointsize)
}

if (opts$jsonify_results) {
  metrics <- list(roc_auc=auc)
  supplementary <- list(roc_curve=roc_curve_as_points_list(roc_data$tpr, roc_data$fpr))
  results <- list(metrics=metrics, supplementary=supplementary)
  writeLines(toJSON(results))
} else{
  writeLines(as.character(auc))
}
