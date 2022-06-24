# motif_pseudo_roc
**Container name:** `vorontsovie/motif_pseudo_roc`

Benchmark for evaluation of weight matrices (PWMs) based on ChIP-seq data (sequences). Control data is simulated based on PWM score distribution.

Original benchmark is developed by Ivan Kulakovskiy.

## Instructions
### Console version
Motif should be supplied either as a mounted file `/motif.pwm` (or `/motif.pcm`, or `/motif.ppm`, or `/motif.pfm`) or via `--motif-url <URL>`. Conversion to PWM is performed internally.

Positive set of sequences can be passed via `/peaks.fa`. Also one can pass set of peak coordinates and an assembly. Peaks can be passed in different formats: `/peaks.bed`, `/peaks.narrowPeak` or just `/peaks`. In the last case one should also specify `--peak-format` which describes which columns contain chromosome, start and end coordinates and how should such a peak be treated. There are several options to treat peaks: they can be taken as is using format `entire`. Alternatively one can take fixed length intervals around peak center/summit; length of interval can be changed used `--peak-flank-size` option (default: 150nt in each direction). It is recommended to use fixed length intervals around peak summit.

For `.bed` format peaks are treated as `--peak-format 1,2,3,center`, i.e. center of each peak is taken.

For `.narrowPeak` format peak are treated as `--peak-format 1,2,3,summit:rel:10`, i.e. summit is taken. Summit coordinate is specified in 10-th column relative to peak start position.

Peak treatment can be redefined. E.g. content of `peaks.interval`:
```
#chromosome	from	to	some_info	summit_position	any_other_infos...
chr1	12	112	.	62	...
chrX	100500	100735	.	100598	...
...
```
then we can process it using `--peak-format 1,2,3,summit:abs:4` option.


Background model is crucial for this benchmark. We recommend to use dinucleotide background models. Default model `infer:di` is obtained by calculating dinculeotide frequencies of positive set. An alternative is to use predefined background, e.g. whole genome frequencies for hg38:
`--background di:0.09815922831551911,0.05049139142847783,0.0700008469320879,0.07700899311854521,0.07264462371481141,0.0516855514432905,0.01000850296695631,0.0700008469320879,0.059702607121927695,0.042460006040716786,0.0516855514432905,0.05049139142847783,0.06515399996155279,0.059702607121927695,0.07264462371481141,0.09815922831551911`.

If you test multiple motifs using prepare/evaluate stage separation, it's reasonable to store inferred background using `--store-background /bg.txt` during prepare stage and to load it from file during evaluation stage using `--background file:/bg.txt`.


### Invocation example:
```
docker run --rm  \
  --volume /path/to/genomes/:/assembly  \
  --volume $(pwd):/data:ro  \
  vorontsovie/motif_pseudo_roc:2.1.0  \
    evaluate \
    --assembly-name hg38 \
    --peaks /data/TF_peaks.bed \
    --motif /data/TF_motif.pwm
```
### Invocation using prepare/evaluate stages
```
docker run --rm -it  \
  --volume /path/to/genomes/:/assembly  \
  --volume $(pwd):/data:ro  \
  vorontsovie/motif_pseudo_roc:2.1.0  \
    /bin/sh

prepare --assembly-name hg38  --peaks /data/TF_peaks.bed \
  --positive-file /pos.fa  --store-background /bg.txt

evaluate --motif /data/TF_motif.pwm \
  --positive-file /pos.fa  --background file:/bg.txt
```

### Version configurable via `config.json` (outdated)
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
