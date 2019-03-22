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

plot_roc <- function(roc_data, output_filename) {
  png(output_filename, width = width, height = height, pointsize = pointsize)
  plot(1 - roc_data$fpr, roc_data$tpr, t="s", col="red", xlab="1-Specificity", ylab="Sensitivity", main="ROC curve of TF motif predictor", xlim=c(1,0), lwd = 4, cex.axis=1.4, cex.main=1.4, cex.lab=1.4, cex.sub=1.4)
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
  make_option(c("--plot"), dest="plot_image", action="store_true", help="Plot ROC curve"),
  make_option(c("--plot-filename"), dest="image_filename", type="character", default="/results/roc_curve.png", help="Specify plot filename [default=%default]"),
  make_option(c("--top"), type="double", default=0.1, help="Fraction of top sequences to take [default=%default]"),
  make_option(c("--bins"), type="integer", default=1000, help="Number of bins for ROC computations [default=%default]"),
  make_option(c("--pseudo-weight"), dest="pseudo_weight", type="double", default=0.0001, help="Set a pseudo-weight to re-normalize the frequencies of the letter-probability matrix (LPM) [default=%default]")
)

opt_parser <- OptionParser(option_list=option_list);
opts_and_args <- parse_args(opt_parser, positional_arguments=TRUE);
opts <- opts_and_args[[1]]
args <- opts_and_args[[2]]

plot_image = opts$plot_image
image_filename = opts$image_filename
pseudo_weight = opts$pseudo_weight
top_fraction = opts$top
n_bins = opts$bins

positive_sequences_fn = args[1]
matrix_fn = args[2]

system(paste("cp", positive_sequences_fn, "/workdir/positive.seq"))
system(paste("cp", matrix_fn, "/workdir/matrix.ppm"))
system(paste("/app/seqshuffle positive.seq > /workdir/negative.seq"))
system(paste("/app/pwm_scoring -r -w", pseudo_weight, "-m matrix.ppm /workdir/positive.seq  > /workdir/positive_PWM.out"))
system(paste("/app/pwm_scoring -r -w", pseudo_weight, "-m matrix.ppm /workdir/negative.seq  > /workdir/negative_PWM.out"))

POS <- log10(as.numeric(read.table("/workdir/positive_PWM.out", header=F)[,1]))
NEG <- log10(as.numeric(read.table("/workdir/negative_PWM.out", header=F)[,1]))

roc_data <- auc_approx(POS, NEG, top_fraction, n_bins)
if (plot_image) {
  plot_roc(roc_data, image_filename)
}

metrics <- list(roc_auc=roc_data$auc)
supplementary <- list(roc_curve=roc_curve_as_points_list(roc_data$tpr, roc_data$fpr))
results <- list(metrics=metrics, supplementary=supplementary)

fileConn <- file("/results/result.json")
writeLines(toJSON(results), fileConn)
close(fileConn)
