#!/usr/bin/env bash
PSEUDO_WEIGHT=0.0001
cat /workdir/config.json | /app/jq -r '.ppm[] | join("\t")' > /workdir/matrix.ppm
/app/pwm_scoring -w $PSEUDO_WEIGHT -r -m matrix.ppm /data/positive.seq  > positive_PWM.out
/app/pwm_scoring -w $PSEUDO_WEIGHT -r -m matrix.ppm /data/negative.seq  > negative_PWM.out
Rscript /app/calculate_roc.R "$@"
