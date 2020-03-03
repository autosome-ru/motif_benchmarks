# predictions_roc
**Container name:** `vorontsovie/predictions_roc`

Benchmark to compute basic metrics of binary classification problem given scores of objects in positive and negative classes. Origin of data doesn't matter, but scores should be supplied from outside.

## Instructions
Actual labels should be mounted as `/benchmark_specific_data/ground_truth.txt`:
```
+
-
-
...
```

Prediction scores should be supplied in `/workdir/config.json` with the same objects ordering as labels have. 
`config.json` format:
```{"predictions": [1.23, -3.4, 0.56, ...]}```

Results are to be stored in `/workdir/persistent/result.json`:
```
{
  "metrics": {
    "roc_auc": roc_auc,
  },
  "supplementary": {
    "roc": [{"fpr": , "tpr": }, {"fpr": , "tpr": },...],
  },
}
```
