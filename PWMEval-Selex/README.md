# PWMEval-Selex
**Container name:** `vorontsovie/pwmeval_selex`

In order to benchmark a motif against SELEX dataset one has to specify a motif itself (positional frequency matrix) and a file with sequences obtained from a SELEX experiment (positive samples).

To run these benchmarks one should supply data files either by mounting local files into container predefined path, or by specifying URL to download these files. As a result it reports ROC AUC at stdout. There are configuration options to return resulting metrics in json format instead of single number (also printed to stdout) and to produce ROC curve plot and coordinates of ROC curve points (FPR and TPR at different thresholds) in separate files.


## Motif format

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

## Sequences format

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

Also note that we mounted our file as `/seq.fq.gz` so that a benchmark know that we have a gzipped FASTQ file, and thus we don't need to specify it explicitly. Instead we could mount sequences to a `/seq` file and to specify `--gz --fastq` options.

## Customization

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

Options `--flank-5` and `--flank-3` concatenate 5' and 3' flanking sequences to provided sequences. It's reasonable to specify adapter and barcode sequences (e.g. 20bp from each side) as they can impact binding. These flanks are concatenated to control dataset after shuffling is done.

## Sequence set precalculation
(version 1.1.0 and later)

When you need to test a bunch of motifs against a single dataset, you can optimize the process by precalculating positive and negative sets of sequences once, and not recreating them again and again for each benchmark invocation.

Let's first prepare positive and negative datasets and then we will use it.
A script to prepare sequences live at entrypoint `/app/prepare_sequences.R`. This script results into two files, by default `/sequences/positive.fa.gz` and `/sequences/negative.fa.gz`. Don't forget to mount container's `/sequences` folder to your host machine, instead you will leave and lose all the results inside a container.
It's reasonable to specify a pair of options: `--positive-file FILENAME` and `--negative-file FILENAME` to make names of resulting files unique along different datasets.

This script accepts all the sequence-oriented options, which accepts the main entrypoint. And it doesn't deal with any motif, so doesn't use any motif-related options.

Usage:
```
docker run --rm \
    --entrypoint /app/prepare_sequences.R  \
    --volume $(pwd)/JUN_dataset.fastq.gz:/seq.fastq.gz:ro \
    --volume $(pwd)/prepared_sequences/:/sequences/  \
    vorontsovie/pwmeval_selex:1.1.0 \
        --positive-file /sequences/JUN_pos.fa.gz \
        --negative-file /sequences/JUN_neg.fa.gz \
        [options]...
```

To get an advantage of these precalculations, the main stage should get these precalculated files. In order to do it, mount these files (you don't need to mount original sequences) and specify options `--positive-file FILENAME` and `--negative-file FILENAME`. The precalculation script writes into these files, the main script reads from them.

Usage:
```
docker run --rm \
    --volume $(pwd)/prepared_sequences/:/sequences/:ro \
    --volume $(pwd)/motif.pcm:/motif.pcm:ro \
    vorontsovie/pwmeval_selex:1.1.0 \
        --positive-file /sequences/JUN_pos.fa.gz \
        --negative-file /sequences/JUN_neg.fa.gz \
        [options]...
```
