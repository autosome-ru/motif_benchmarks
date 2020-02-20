# TFBS-motifs benchmarks
Here we present several containerized benchmarks which can be used to assess performance of positional frequency matrices on ChIP-seq, SELEX and PBM datasets.

More info can be found in dedicated sections.

## Available benchmarks
* [PWMEval-Selex](https://github.com/autosome-ru/motif_benchmarks/blob/master/PWMEval-Selex/README.md) -- benchmark for evaluation of frequency matrices (PFMs) based on SELEX/SELEX-like data. Control data is obtained by shuffling sequences. Container implements almost the same approach as used in a [web-service](https://ccg.epfl.ch/pwmtools/pwmeval_selex.php). Original bechmark developed by Giovanna Ambrosini. Image: *vorontsovie/pwmeval_selex*.
* [PWMEval-Chip-peak](https://github.com/autosome-ru/motif_benchmarks/blob/master/PWMEval-Chip-peak/README.md) -- benchmark for- evaluation of frequency matrices (PFMs) based on ChIP-seq data (peak sets). Control data is obtained by taking 'shades' of peaks at some distance from actual peak. Container implements the same approach as used in a [web-service](https://ccg.epfl.ch/pwmtools/pwmeval_chippeak.php).  Original bechmark developed by Giovanna Ambrosini. Image: *vorontsovie/pwmeval_chipseq.
* [PWMBench-PBM](https://github.com/autosome-ru/motif_benchmarks/blob/master/PWMBench-PBM/README.md) -- benchmark for evaluation of frequency matrices (PFMs) based on PBM data. Original bechmark developed by Jan Grau.  Image: *vorontsovie/pwmbench_pbm*.
* [motif_pseudo_roc](https://github.com/autosome-ru/motif_benchmarks/blob/master/motif_pseudo_roc/README.md) -- benchmark for evaluation of weight matrices (PWMs) based on ChIP-seq data (sequences). Control data is simulated based on PWM score distribution. Original benchmark is developed by Ivan Kulakovskiy. Image: *vorontsovie/motif_pseudo_roc*.
* [predictions_roc](https://github.com/autosome-ru/motif_benchmarks/blob/master/predictions_roc/README.md) -- benchmark to computate basic metrics of binary classification problem given scores of objects in positive and negative classes. Origin of data doesn't matter. Image: *vorontsovie/predictions_roc*.

## General instructions

To run a container it's necessary to have docker installed (and user should have permissions to run a container). Each benchmark is containerized into a separate docker image.

To obtain the latest version of benchmark, use `docker pull vorontsovie/{benchmark_name}:latest` command.

To run a benchmark one should run smth like this command:
`docker run --rm  -v $(pwd)/NFKB2_TTCAAT20NGA_R_2.fastq.gz:/seq.fq.gz  -v $(pwd)/NFKB2_HUMAN.H11MO.0.B.pcm:/motif.pcm vorontsovie/pwmeval_selex --top 0.1`

Let's sort out command structure. General pattern is as follows:

`docker run [docker options] <image name> [benchmark options]`

Image name (`vorontsovie/pwmeval_selex` in this case) specifies type of benchmark.

There're several docker options you should use for these benchmarks. A docker option `--rm` means that we want to automatically remove a container after it finished its work. Note that these containers are not designed for reuse because names of files in different runs will clash! So there are no reasons not to use `--rm` option. The other options are typically used to pass data into container (see below).

### How to supply data

It's critical to understand that a container doesn't have access to an external file system unless you manually mount some external directories as internal ones. So `-v $(pwd)/filename:/filename/in/container` mantra means that we want to take a file `filename` from current folder — which is returned by `$(pwd)` — and mount it to a `/filename/in/container` file inside a container. Note that we should specify absolute pathes both in local file system and inside a container.

In this example we supplied a file with sequences from SELEX experiment and a motif to benchmark as `/seq.fq.gz` and `/motif.pcm` inside a container. These pathnames are standardized, so in each run you will typically mount data at the same paths (or choose among several alternative names for different formats). But for each benchmark necessary files are different, so there is no standardization between runs yet.

Some benchmarks (`PWMEval-Selex` and `PWMEval-Chip-peak`) provide options to pass data via url instead of mount point. In this case, data will be downloaded from specified URLs.

Standardization of input formates is one of directions for further improvements. We are going to develop several containers following convention that input data is supplied via `/workdir/config.json`, yielding results into `/workdir/results.json` and loading data from several specific folders. Take a look at `motif_pseudo_roc` and `predictions_roc` benchmarks) for some examples.

## Development notes

Benchmarks are developed using alpine base image to be as slim as possible. Github repository contains all the necessary files to rebuild these images.

The goal of this initiative - is to build a reusable and standardized set of benchmarks, so after development stage is over, we will use docker image release tags to specify certain versions of a protocol.

Feel free to experiment with your own benchmarking protocols and contribute.
