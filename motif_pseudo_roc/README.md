# motif_pseudo_roc
**Container name:** `vorontsovie/motif_pseudo_roc`

Benchmark for evaluation of weight matrices (PWMs) based on ChIP-seq data (sequences). Control data is simulated based on PWM score distribution.

Original benchmark is developed by Ivan Kulakovskiy.

## Instructions
PWM should be supplied in `/workdir/config.json`. 

`config.json` format:
```{"motif": ">motif_name\n1.2\t34\t-5.0\t6.78\n4\t3\t2\t1\n..."}```

FASTA file with peak sequences should be mounted as `/benchmark_specific_data/control.formatted.mfa`.

Results are to be stored in `/workdir/persistent/result.json`:
```
{
  "metrics": {
    "roc_auc": roc_auc,
    "logroc_auc": logroc_auc,
  },
  "supplementary": {
    "roc": [{"fpr": , "tpr": }, {"fpr": , "tpr": },...],
    "logroc": [{"logfpr": , "tpr": },...],
  },
}
```
