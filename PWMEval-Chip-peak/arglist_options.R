arglist_assembly_options = list(
  make_option(c("--assembly-name"), dest='assembly_name', action="store", type='character', default=NA, help="Choose assembly by name"),
  make_option(c("--assembly-fasta"), dest='assembly_fasta_fn', action="store", type='character', default=NA, help="Choose assembly FASTA file"),
  make_option(c("--assembly-sizes"), dest='assembly_sizes_fn', action="store", type='character', default=NA, help="Choose assembly chromosome sizes file")
)

arglist_peaks_options = list(
  make_option(c("--peaks"), dest= 'peaks_fn', type='character', default=NA, help="Specify peaks file"),
  make_option(c("--peaks-url"), dest= 'peaks_url', type='character', default=NA, help="Use peaks file located at some URL"),

  make_option(c("--narrowPeak"), dest='peak_format_narrowPeak', default=FALSE, action="store_true", help="Peaks are formatted in narrowPeak (peaks are reshaped into constant-size peaks around summit of a peak)"),
  make_option(c("--bed"), dest='peak_format_bed', default=FALSE, action="store_true", help="Peaks are formatted in bed (peaks are reshaped into constant-size peaks around center of a peak)"),

  make_option(c("--gz"), dest="compression_gz_peaks", default=FALSE, action="store_true", help="Force un-gzipping peaks"),
  make_option(c("--not-compressed"), dest="compression_no_peaks", default=FALSE, action="store_true", help="Prevent un-gzipping peaks"),

  make_option(c("--top"), type="integer", dest="num_top_peaks", default=500, help="Number of top peaks to take [default=%default]")
)

arglist_motif_options = list(
  make_option(c("--motif"), dest='motif_fn', type='character', default=NA, help="Specify motif file"),
  make_option(c("--motif-url"), dest= 'motif_url', type='character', default=NA, help="Use PFM file located at some URL"),

  make_option(c("--pfm"), default=FALSE, action="store_true", help="Force use of PFM matrix"),
  make_option(c("--pcm"), default=FALSE, action="store_true", help="Force use of PCM matrix")
)
