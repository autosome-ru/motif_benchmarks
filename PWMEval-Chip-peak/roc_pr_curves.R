plot_curve <- function(curve_data, image_filename, width = 800, height = 800) {
  dir.create(dirname(image_filename), recursive=TRUE, showWarnings=FALSE)
  png(image_filename, width = width, height = height)
  plot(curve_data, lwd = 4, cex.axis=1.4, cex.main=1.4, cex.lab=1.4, cex.sub=1.4)
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
  dir.create(dirname(output_filename), recursive=TRUE, showWarnings=FALSE)
  write.table(list(fpr=roc_data$fpr, tpr=roc_data$tpr), row.names=FALSE, quote=FALSE, sep="\t", file=output_filename)
}

store_pr <- function(pr_data, output_filename) {
  dir.create(dirname(output_filename), recursive=TRUE, showWarnings=FALSE)
  write.table(list(recall=pr_data$recall, precision=pr_data$precision), row.names=FALSE, quote=FALSE, sep="\t", file=output_filename)
}
