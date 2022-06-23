# PWMEval-Chip-peak
**Container name:** `vorontsovie/pwmeval_chipseq`

ChIP-seq benchmark needs three types of files to run: motif itself (positional frequency matrix), list of ChIP-seq peaks with scores in BED or narrowPeak format, and a genome assembly in a single FASTA file.

To run these benchmarks one should supply data files either by mounting local files into container predefined path, or by specifying URL to download these files. As a result it reports ROC AUC at stdout. There are configuration options to return resulting metrics in json format instead of single number (also printed to stdout) and to produce ROC curve plot and coordinates of ROC curve points (FPR and TPR at different thresholds) in separate files.

Example:
```
docker run -v /path/to/genomes/:/assembly/ -v /path/to/data:/data vorontsovie/pwmeval_chipseq:1.1.2 evaluate --assembly-name hg38 --motif /data/motif.pfm --peaks /data/peaks.narrowPeak [benchmark options]
```

Motif is specified exactly the same way as it's [specified](https://github.com/autosome-ru/motif_benchmarks/blob/master/PWMEval-Selex/README.md) in PWMEval-Selex benchmark.

## Peaks format

Peaks can be passed by mounting a file with peaks and passing it via option `--peaks /path/to/peaks[.bed|.narrowPeak][.gz]` or by specifying a link `--peaks-url URL`. When format can't be derived by extension, it defaults to BED.

All peaks are transformed to have a standard length, but BED peaks and narrowPeaks are transformed in a different way:
* Peaks in BED format produce a symmetrical interval of constant length around peak center. The 5-th column is used as a peak score.
* Peaks in narrowPeak format produce a symmetrical interval of constant length around peak summit (relative position of summit in an interval is specified in the 10-th column). The 7-th column is used as a peak score.

One can specify peaks format explicitly using `--narrowPeak` or `--bed` options in combination with `--gz` and `--not-compressed`.

## Passing genome assembly

To pass an assembly, one should specify an assembly name using `--assembly-name` option. Assembly name normally is a short UCSC identifier such as hg38. Not to download an assembly each time, it's highly recommended to mount a folder `/assembly` with files named smth like `hg38.fa`, `hg38.chrom.sizes`, `hg38.fa.fai`. Missing files will be downloaded and/or generated automatically — and stored in this folder — but that can be a time consuming process during the first run.

It's not guaranted that a tool will successfully obtain an assembly from UCSC, in some cases you'd do this manually. Chromosome sizes and index always can be reconstructed from FASTA file, so you don't have to download them too.

If you don't want to specify an assembly name, you can mount your assembly smth like `/assembly/genome.fa` and a family (accordingly named files with an index and chromosome sizes) and to pass an option `--assembly-fasta /assembly/genome.fa`. But note that supplementary files are only stored between runs if `/assembly` folder is mounted in a local file system and names of supplementary files are derived from the name of assembly FASTA file (so `genome.fa` will result in not very meaningful names like `genome.fa.fai`  and `genome.chrom.sizes`).

## Customization

Options `--plot-roc` / `--plot-roc-filename`, `--roc` / `--roc-filename` (and their Precision-Recall `pr` counterparts) and `--json` work the same way as they do for `PWMEval-Selex` benchmark.

Option `--top` behaves *differently* from option with same name in `PWMEval-Selex`. Here it indicates a number of peaks which should be taken as a positive dataset, only best peaks (score is specified in bed/narrowPeak file, see details above) are taken into account.

## Prepare/Evaluate

It's possible to separate sequence preparation stage from motif evaluation stage. We don't need motif on the preparation stage and don't need peaks/assembly on the evaluation stage. These stages are glued together with `--positive-file` and `--negative-file` options. Such separation can be useful when one want to test dozens of motifs against the same dataset.
```
docker run \
    -v /path/to/genomes/:/assembly/  -v /path/to/data:/data  -v /path/to/temporary/storage:/sequences \
    vorontsovie/pwmeval_chipseq:1.1.2  prepare \
        --assembly-name hg38  --peaks /data/peaks.narrowPeak \
        --positive-file /sequences/positive.fa --negative-file /sequences/negative.fa
docker run \
    -v /path/to/data:/data  -v /path/to/temporary/storage:/sequences \
    vorontsovie/pwmeval_chipseq:1.1.2  evaluate \
        --motif /data/motif.pfm \
        --positive-file /sequences/positive.fa --negative-file /sequences/negative.fa
```
