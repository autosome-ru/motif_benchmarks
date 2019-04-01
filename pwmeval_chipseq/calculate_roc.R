#!/usr/bin/env Rscript
library('optparse')
library('MASS')
library('caTools')
suppressPackageStartupMessages(library('pROC'))
source('/app/utils.R')
source('/app/motif_preprocessing.R')
source('/app/peak_preprocessing.R')

plot_roc <- function(roc_data, image_filename, width = 800, height = 800, pointsize = 20) {
  png(image_filename, width = width, height = height, pointsize = pointsize)
  plot(roc_data, main="ROC curve of motif predictor", col = "red", lwd = 4, cex.axis=1.4, cex.main=1.4, cex.lab=1.4, cex.sub=1.4)
  txt <- paste("AUC = ", roc_data$auc, sept="")
  legend(x="bottomright", legend=txt, col=c("red"), lty=1, lwd=4, cex=1.2)
  dev.off()
}

width = 800
height = 800
pointsize = 20

option_list = list(
  make_option(c("--assembly-name"), dest='assembly_name', action="store", type='character', default=NA, help="Choose assembly by name"),
  make_option(c("--motif-url"), dest= 'motif_url', type='character', default=NA, help="Use PPM file located at some URL"),
  make_option(c("--peaks-url"), dest= 'peaks_url', type='character', default=NA, help="Use peaks file located at some URL"),

  make_option(c("--ppm"), default=FALSE, action="store_true", help="Force use of PPM matrix"),
  make_option(c("--pcm"), default=FALSE, action="store_true", help="Force use of PCM matrix"),

  make_option(c("--narrowPeak"), dest='peak_format_narrowPeak', default=FALSE, action="store_true", help="Peaks are formatted in narrowPeak (peaks are reshaped into constant-size peaks around summit of a peak)"),
  make_option(c("--bed"), dest='peak_format_bed', default=FALSE, action="store_true", help="Peaks are formatted in bed (peaks are reshaped into constant-size peaks around center of a peak)"),

  make_option(c("--gz"), dest="compression_gz_peaks", default=FALSE, action="store_true", help="Force un-gzipping peaks"),
  make_option(c("--not-compressed"), dest="compression_no_peaks", default=FALSE, action="store_true", help="Prevent un-gzipping peaks"),

  make_option(c("--plot"), dest="plot_image", default=FALSE, action="store_true", help="Plot ROC curve"),
  make_option(c("--plot-filename"), dest="image_filename", type="character", default="/results/roc_curve.png", help="Specify plot filename [default=%default]"),
  make_option(c("--top"), type="integer", dest="num_top_peaks", default=500, help="Number of top peaks to take [default=%default]")
)

opt_parser <- OptionParser(option_list=option_list);
opts_and_args <- parse_args(opt_parser, positional_arguments=TRUE);
opts <- opts_and_args[[1]]
args <- opts_and_args[[2]]

dummy = obtain_and_preprocess_motif(opts)
dummy = obtain_and_preprocess_peaks(opts)

if (dir.exists("/assembly")) {
  if (is.na(opts$assembly_name)) {
    stop("Error! Specify assembly name")
  }
  assembly_fasta_fn = file.path("/assembly", paste0(opts$assembly_name, ".fa"))
  assembly_sizes_fn = file.path("/assembly", paste0(opts$assembly_name, ".sizes"))
  if (!file.exists(assembly_fasta_fn)) {
    system(paste("/app/download_assembly.sh", opts$assembly_name, shQuote(assembly_fasta_fn)))
  }
} else{
  if (file.exists("/assembly.fa")) {
    assembly_fasta_fn = "/assembly.fa"
    assembly_sizes_fn = "/assembly.chrom.sizes"
  } else {
    stop("Mount /assembly.fa file (also /assembly.chrom.sizes and /assembly.fa.fai not to recalculate them)")
  }
}

if (!file.exists(assembly_sizes_fn)) {
  writeLines(paste("Chromosome sizes file", assembly_sizes_fn, "not found, generating..."), con=stderr())
  system(paste("/app/chrom_sizes", shQuote(assembly_fasta_fn), " > ", shQuote(assembly_sizes_fn)))
}

system(paste("sort -k5,5nr /workdir/peak_centers_scored.bed | head -n", opts$num_top_peaks, " > /workdir/top_peaks.bed"))
system(paste("/app/bedtools slop -i /workdir/top_peaks.bed -g ", shQuote(assembly_sizes_fn), " -l 124 -r 125  > /workdir/positive_peaks.bed"))
system(paste("/app/bedtools slop -i /workdir/top_peaks.bed -g ", shQuote(assembly_sizes_fn), " -l -301 -r 550  > /workdir/negative_peaks.bed"))
system(paste("/app/bedtools getfasta -bed /workdir/positive_peaks.bed -fi ", shQuote(assembly_fasta_fn), "  > /workdir/positive.seq"))
system(paste("/app/bedtools getfasta -bed /workdir/negative_peaks.bed -fi ", shQuote(assembly_fasta_fn), "  > /workdir/negative.seq"))

system("/app/pwm_scoring -r -u -m /workdir/motif.ppm /workdir/positive.seq  > /workdir/positive_PWM.out")
system("/app/pwm_scoring -r -u -m /workdir/motif.ppm /workdir/negative.seq  > /workdir/negative_PWM.out")

pos <- as.matrix(read.table("/workdir/positive_PWM.out"))
neg <- as.matrix(read.table("/workdir/negative_PWM.out"))
W = wilcox.test(pos, neg, alternative ="g")$statistic
AUC = W/length(pos)/length(neg)
pos_labs <- rep(1, 500)
neg_labs <- rep(0, 500)
pos <- cbind(pos, pos_labs)
neg <- cbind(neg, neg_labs)
comb_sets <- rbind(pos, neg)
roc_data <- roc(response = comb_sets[, 2], predictor = comb_sets[, 1])
auc = as.numeric(roc_data$auc) # by default auc is a complex object which is displayed not as a number but number + text
writeLines(as.character(auc))

if (opts$plot_image) {
  plot_roc(roc_data, opts$image_filename, width = width, height = height, pointsize = pointsize)
}
