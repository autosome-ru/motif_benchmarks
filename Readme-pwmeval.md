Here I present two containerized benchmarks (vorontsovie/pwmeval_selex and vorontsovie/pwmeval_chipseq) which can be used to benchmark positional frequency matrices on SELEX datasets and on ChIP-seq datasets. In order to run a container one need only a docker installed.

To get the last version of benchmarks one should invoke
```
docker pull vorontsovie/pwmeval_selex
docker pull vorontsovie/pwmeval_chipseq
```

To run these benchmarks one should supply data files either by mounting local files into container predefined path, or by specifying URL to download these files. As a result they report ROC AUC at stdout. There're configuration options to return resulting metrics in json format instead of single number (also printed to stdout) and to produce ROC curve plot and coordinates of ROC curve points (FPR and TPR at different thresholds) in separate files.

## Benchmark using SELEX data

In order to benchmark a motif on SELEX one has to specify a motif itself (positional frequency matrix) and a file with sequences obtained from a SELEX experiment (positive samples).
A motif can be given as a positional count matrix (with .pcm extension) or as a frequency matrix (with .pfm or .ppm extension). By default tool tries to guess matrix type from its extension, if an extension doesn't tell us meaningful information about a matrix type, frequency matrix is supposed. It's also possible to force usage of certain type by specifying options `--pcm/--pfm`.
Note that name and extension of a file in a local file system are not taken into account, only name of a mount point does matter.
There're several mount points at which benchmark will try to find a matrix: `/motif[.pfm|.ppm|.pcm]`. Square brackets mean that an extension can be omitted.
For example lets pass a file into a container with an option `-v` (mount a volume):
```
docker run -v /absolute/path/to/CTCF_HUMAN.H11MO.0.A.pcm:/motif.pcm [other docker options] vorontsovie/pwmeval_selex [benchmark options]
```

Another way to pass a matrix is via `--motif-url URL` option. In this case we don't have to mount a file
```
docker run [docker options] vorontsovie/pwmeval_selex --motif-url 'http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/pcm/CTCF_HUMAN.H11MO.0.A.pcm' [benchmark options]
```

Matrix format is a very simple plain text 4-column matrix (ACGT-ordered):
```
	[> optional header]
	1st position counts or probabilities in order A, C, G, T (as columns)
	2nd position
	...
```

SELEX sequences also can be passed by mounting them to one of `/seq[.fasta|.fa|.fastq|.fq][.gz]`
As one can guess, a tool accepts FASTA (.fasta and .fa extensions) and FASTQ files (.fastq and .fq extensions). Default format is FASTA. Sequences can be gzipped (it's derived from .gz extension). To override guessed format one can use `--fasta/--fastq` options in combination with `--gz/--not-compressed` options. To pass sequences via URL use `--seq-url` option.

Let's now run a complete example with passing everything by URL:

```
docker run --rm vorontsovie/pwmeval_selex --seq-url 'ftp.sra.ebi.ac.uk/vol1/run/ERR194/ERR194820/NFKB2_TTCAAT20NGA_R_2.fastq.gz' --motif-url 'http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/pcm/NFKB2_HUMAN.H11MO.0.B.pcm'
```

and by passing everything via mounting point:
```
wget 'ftp.sra.ebi.ac.uk/vol1/run/ERR194/ERR194820/NFKB2_TTCAAT20NGA_R_2.fastq.gz'
wget 'http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/pcm/NFKB2_HUMAN.H11MO.0.B.pcm'
docker run --rm  -v $(pwd)/NFKB2_TTCAAT20NGA_R_2.fastq.gz:/seq.fq.gz  -v $(pwd)/NFKB2_HUMAN.H11MO.0.B.pcm:/motif.pcm vorontsovie/pwmeval_selex
```

A docker option `--rm` means that we want to automatically remove a container after it finished its work. Note that these containers are not designed for reuse because names of files in different runs will clash! So there are no reasons not to use `--rm` option.

`-v $(pwd)/filename:/filename/in/container` matra means that we want to take a file `filename` from current folder — which is returned by `$(pwd)` — and mount it to a `/filename/in/container` file inside a container. Note that we should specify absolute pathes both in local file system and inside a container.

Also note that we mounted our file as `/seq.fq.gz` so that a benchmark know that we have a gzipped FASTQ file, and thus we don't need to specify it explicitly. Instead we could mount sequences to a `/seq` file and to specify `--gz --fastq` options.

Now lets discuss some options to change benchmarking procedure and to provide more information.

When option `--json` is specified, the benchmark's output is not just an ROC AUC value but a JSON object. This format allows one to extend benchmark with several other metrics such as AUPRC or recall@0.1. Also in supplementry section of this JSON object there are a ROC curve points.

The next pair of options requires user to have some folder to be mounted to `/results`. E.g. one can add smth like `-v $(pwd):/results` - in this case resulting files will be stored in current folder.
When `--plot` option is specified, benchmark generates a ROC curve plot `roc_curve.png` in that folder. Filename can be overrided using `--plot-filename` option.
When `--roc` option is specified, benchmark generates a table `roc_curve.tsv` in a tab-separated values format. This file writes out points of ROC curve in (FPR, TPR) coordinates. `--roc-filename` can change a name of resulting file.

Option `--seq-length L` tells that all sequences of different length should be rejected from a file with sequences.
Option `--allow-iupac` allows sequences to have N-nucleotide, by default all such sequences are rejected.
Option `--non-redundant` takes only unique sequences in benchmark.

Option `--top FRACTION` specifies a fraction of top-scoring sequences which should be taken into account. FRACTION is a number in [0,1] range with default value of 0.1
Option `--bins N` specifies number of bins in use for ROC-curve construction. Default value is 1000.
Option `--pseudo-weight W` specifies a pseudoweight to be added to a PFM. Default value is 0.0001.

As negative control data is obtained by random shuffling of a list of positive samples, metrics can differ a bit from run to run. In order to get randomness out, one can specify random number generator seeding value with a `--seed INT` option.

## Benchmark using ChIP-seq data
ChIP-seq benchmark needs three types of files to run: motif itself (also positional frequency matrix), list of ChIP-seq peaks with scores in BED or narrowPeak format, and a genome assembly in a single FASTA file.

Motif is specified exactly the same way as it's specified in a SELEX benchmark.

Peaks can be passed by mounting a file with peaks as `/peaks[.bed|.narrowPeak][.gz]` or by specifying a link `--peaks-url URL`. When format can't be derived by extension, it defaults to BED.
All peaks are transformed to have a standard length, but BED peaks and narrowPeaks are transformed in a different way:
* Peaks in BED format produce a symmetrical interval of constant length around peak center. The 5-th column is used as a peak score.
* Peaks in narrowPeak format produce a symmetrical interval of constant length around peak summit (relative position of summit in an interval is specified in the 10-th column). The 7-th column is used as a peak score.
One can specify peaks format explicitly using `--narrowPeak` or `--bed` options in combination with `--gz` and `--not-compressed`.

To pass an assembly, one should specify an assembly name using `--assembly-name` option. Assembly name normally is a short UCSC identifier such as hg38. Not to download an assembly each time, it's highly recommended to mount a folder `/assembly` with files named smth like `hg38.fa`, `hg38.chrom.sizes`, `hg38.fa.fai`. Missing files will be downloaded and/or generated automatically — and stored in this folder — but that can be a time consuming process during the first run.
It's not guaranted that a tool will successfully obtain an assembly from UCSC, in some cases you'd do this manually. Chromosome sizes and index always can be reconstructed from FASTA file, so you don't have to download them too.

If you don't want to specify an assembly name, you can mount your assembly as `/assembly.fa` or `/assembly/assembly.fa` and a family (accordingly named files with an index and chromosome sizes). But note that supplementary files are only stored between runs if `/assembly` folder is mounted in a local file system? and names of supplementary files are derived from the name of assembly FASTA file (so assembly.fa will result in not very meaningful names like `assembly.fa.fai`  and `assembly.chrom.sizes`).

Options `--plot` / `--plot-filename`, `--roc` / `--roc-filename` and `--json` work the same way as they do for `pwmeval_selex` benchmark.

Option `--top` behaves differently from option with same name in `pwmeval_selex`. Here it indicates number of peaks which should be taken as a positive dataset, only best peaks (score is specified in bed/narrowPeak file, see details above) are taken into account.
