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


# wgEncodeAwgTfbsSydhGm12891NfkbTnfaIggrabUniPk.narrowPeak | awk -e '{print $1 "\t" ($2+$10) "\t" ($2+$10+1) "\t" "." "\t" $7 }' | head

  sort -k2,2nr $MGAdir/hg19/remap/ENCSR000EAI.RELA.GM12891.narrowPeak.fps | sed -n 1,500p> $workDir/$pwmevalFile.fps;
  EPDtoFA -l-124 -r125 < $workDir/$pwmevalFile.fps 2>/dev/null |perl -pe 's/EP /EP_/' > $workDir/motif_hg19_35623_pos.seq;
  EPDtoFA -l301 -r550 < $workDir/$pwmevalFile.fps 2>/dev/null |perl -pe 's/EP /EP_/' > $workDir/motif_hg19_35623_neg.seq


  pwm_scoring -r -u -m $workDir/pwmFile $workDir/motif_hg19_35623_pos.seq  > $workDir/motif_hg19_35623_pos_PWM.out;
  pwm_scoring -r -u -m $workDir/ $workDir/motif_hg19_35623_neg.seq  > $workDir/motif_hg19_35623_neg_PWM.out


pos <- as.matrix(read.table("$workDir/motif_hg19_35623_pos_PWM.out"))
neg <- as.matrix(read.table("$workDir/motif_hg19_35623_neg_PWM.out"))
W = wilcox.test(pos, neg, alternative ="g")$statistic
AUC = W/length(pos)/length(neg)
pos_labs <- rep(1, 500)
neg_labs <- rep(0, 500)
pos <- cbind(pos, pos_labs)
neg <- cbind(neg, neg_labs)
comb_sets <- rbind(pos, neg)
roc_data <- roc(response = comb_sets[, 2], predictor = comb_sets[, 1])
write.table(roc_data$auc, col.names=F, row.names=T, file="$workDir/motif_hg19_35623_AUC.out", quote=F, sep=" ")
plot_roc(roc_data, "/results/roc_curve.png")

