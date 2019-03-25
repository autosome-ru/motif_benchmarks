#!/usr/bin/env Rscript

library('MASS')
library('caTools')
library('pROC')

plot_roc <- function(roc_data, image_filename) {
  png(image_filename, width = 800, height = 800, pointsize=20)
  plot(roc_data, main="ROC curve of motif predictor", col = "red", lwd = 4, cex.axis=1.4, cex.main=1.4, cex.lab=1.4, cex.sub=1.4)
  txt <- paste("AUC = ", roc_data$auc, sept="")
  legend(x="bottomright", legend=txt, col=c("red"), lty=1, lwd=4, cex=1.2)
  dev.off()
}

option_list = list(
  make_option(c("--plot"), dest="plot_image", default=FALSE, action="store_true", help="Plot ROC curve"),
  make_option(c("--plot-filename"), dest="image_filename", type="character", default="/results/roc_curve.png", help="Specify plot filename [default=%default]"),
  make_option(c("--top"), type="integer", default=500, help="Number of top peaks to take [default=%default]"),
)

opt_parser <- OptionParser(option_list=option_list);
opts_and_args <- parse_args(opt_parser, positional_arguments=TRUE);
opts <- opts_and_args[[1]]
args <- opts_and_args[[2]]

plot_image = opts$plot_image
image_filename = opts$image_filename
# pseudo_weight = opts$pseudo_weight
num_top_peaks = opts$top
# n_bins = opts$bins

peaks_fn = args[1]
matrix_fn = args[2]

system(paste("ln -s", peaks_fn, "/workdir/peaks.narrowPeak"))
system(paste("ln -s", matrix_fn, "/workdir/matrix.ppm"))

system("/app/narrowpeak2bed /workdir/peaks.narrowPeak > /workdir/scored_peaks.bed")
system("sort -k5,5nr /workdir/scored_peaks.bed | head -n", num_top_peaks, " > /workdir/top_peaks.bed")
system("/app/bedtools slop -i /workdir/top_peaks.bed -g /data/genome.sizes -l 124 -r 125  > /workdir/positive_peaks.bed")
system("/app/bedtools slop -i /workdir/top_peaks.bed -g /data/genome.sizes -l -301 -r 550  > /workdir/negative_peaks.bed")
system("/app/bedtools getfasta -bed /workdir/positive_peaks.bed -fi /data/genome.fa  > /workdir/positive.seq")
system("/app/bedtools getfasta -bed /workdir/negative_peaks.bed -fi /data/genome.fa  > /workdir/negative.seq")

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
# write.table(roc_data$auc, col.names=F, row.names=T, file="$workDir/motif_hg19_35623_AUC.out", quote=F, sep=" ")
writeLines(roc_data$auc)

if (plot_image) {
  plot_roc(roc_data, image_filename)
}
