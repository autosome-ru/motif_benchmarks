#!/usr/bin/env Rscript
library('optparse')
library('MASS')
library('caTools')
library('pROC')

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
  make_option(c("--plot"), dest="plot_image", default=FALSE, action="store_true", help="Plot ROC curve"),
  make_option(c("--plot-filename"), dest="image_filename", type="character", default="/results/roc_curve.png", help="Specify plot filename [default=%default]"),
  make_option(c("--top"), type="integer", dest="num_top_peaks", default=500, help="Number of top peaks to take [default=%default]")
)

opt_parser <- OptionParser(option_list=option_list);
opts_and_args <- parse_args(opt_parser, positional_arguments=TRUE);
opts <- opts_and_args[[1]]
args <- opts_and_args[[2]]

system(paste("cp /peaks.narrowPeak /workdir/peaks.narrowPeak"))
system(paste("cp /motif.ppm /workdir/matrix.ppm"))

if (dir.exists("/assembly")) {
  assembly_fasta_fn = file.path("/assembly", opts$assembly_name)
  assembly_sizes_fn = file.path("/assembly", paste(opts$assembly_name, ".sizes", sep=""))
  if (!file.exists(assembly_fasta_fn)) {
    # wget it
  }
} else{
  if (file.exists("/assembly.fa")) {
    assembly_fasta_fn = "/assembly.fa"
    assembly_sizes_fn = "/assembly.chrom.sizes"
  } else {
    simpleError("Mount /assembly.fa file (also /assembly.chrom.sizes and /assembly.fa.fai not to recalculate them)")
  }
}

if (!file.exists(assembly_sizes_fn)) {
  writeLines(paste("Chromosome sizes file", assembly_sizes_fn, "not found, generating..."), con=stderr())
  system(paste("/app/chrom_sizes", shQuote(assembly_fasta_fn), " > ", shQuote(assembly_sizes_fn)))
}

system("/app/narrowpeak2bed /workdir/peaks.narrowPeak > /workdir/scored_peaks.bed")
system(paste("sort -k5,5nr /workdir/scored_peaks.bed | head -n", opts$num_top_peaks, " > /workdir/top_peaks.bed"))
system("/app/bedtools slop -i /workdir/top_peaks.bed -g /assembly.chrom.sizes -l 124 -r 125  > /workdir/positive_peaks.bed")
system("/app/bedtools slop -i /workdir/top_peaks.bed -g /assembly.chrom.sizes -l -301 -r 550  > /workdir/negative_peaks.bed")
system("/app/bedtools getfasta -bed /workdir/positive_peaks.bed -fi /assembly.fa  > /workdir/positive.seq")
system("/app/bedtools getfasta -bed /workdir/negative_peaks.bed -fi /assembly.fa  > /workdir/negative.seq")

system("/app/pwm_scoring -r -u -m /workdir/matrix.ppm /workdir/positive.seq  > /workdir/positive_PWM.out")
system("/app/pwm_scoring -r -u -m /workdir/matrix.ppm /workdir/negative.seq  > /workdir/negative_PWM.out")

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
