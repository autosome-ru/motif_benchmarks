# PWMBench-PBM
**Container name:** `vorontsovie/pwmbench_pbm`

Benchmark of PBM data developed by Jan Grau. Based on [jstacs](http://www.jstacs.de/index.php/Main_Page) project.

## Instructions

How to run benchmark:
```
docker run --rm --mount type=bind,src=$(pwd),dst=/data vorontsovie/pwmbench_pbm <METRIC> /data/pbm_data.txt /data/motif.mat
```

METRIC corresponds to different metrics.
* ASIS - correlation of intensity (taken as is) and log-sum-occupancy.
* EXP - correlation of intensity and sum-occupancy.
* LOG - correlation of log-intensity and log-sum-occupancy.
* ROC - AUC ROC
* PR - AUC PR
* ROCLOG - AUC ROC in logarithmed intensities
* PRLOG - AUC PR in logarithmed intensities
* all - calculate every metrics ({metrics_name: value, ...} in JSON format is printed)

Motif score is defined as sum-occupancy.

Motif have to be positional frequency matrix or positional count matrix. It's internally corrected with pseudocount of 0.0001.

Motif format (`motif.mat`):
```
>AHR_HUMAN.H11MO.0.B
0.26314477735061276	0.11859745819448689	0.36641340309888065	0.25184436135601973
0.07065186803239253	0.07710428433167413	0.2251456131484061	0.6270982344875272
0.14105450572642986	0.2850314719772257	0.13449544731512686	0.4394185749812176
0.0165402342545919	0.008555293579961341	0.9474208905975482	0.027483581567898402
0.0	0.9766155308814389	0.009695999390622859	0.013688469727938113
0.022350405419921632	0.005133176147976808	0.9702350068107785	0.0022814116213230235
0.0	0.02235040541992163	0.004562823242646056	0.9730867713374324
0.0	0.0	1.0	0.0
0.279810224889837	0.43437346977265195	0.10496230952256	0.180853995814951
```

PBM data format (`pbm_data.txt`):
```
370241.822966	GAACAATGTAAATTATTGAAAGGGCTAATTCAATTAGTCTGTGTTCCGTTGTCCGTGCTG
365545.496066	GATAACCGACGCCCATTAATTATATTAGCATTGAGCGTCTGTGTTCCGTTGTCCGTGCTG
353259.608349	ATTGATTGATGGCTAACTAAATTAAGCGCATGGAGGGTCTGTGTTCCGTTGTCCGTGCTG
347894.963074	CGTCTATTTTCGGGTAATTATCTCATAATGAGGTGGGTCTGTGTTCCGTTGTCCGTGCTG
334081.174145	GATGCTCCGGATTATTAAGTAATTAAATGAGTTTCCGTCTGTGTTCCGTTGTCCGTGCTG
...
```
