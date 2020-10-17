require 'json'
require 'optparse'

require_relative 'utils'
require_relative 'motif_preprocessing'
require_relative 'peak_preprocessing'
require_relative 'assembly_preprocessing'
# require 'bioinform'


options = {
#   image_filename: "/results/roc_curve.png",
#   roc_filename: "roc_curve.tsv",
#   num_top_peaks: 500,
#   plot_curve: false,
  flank_size: 150,
  pseudocount: :log,
  word_count: 100,
  curve_points: false,
}

option_parser = OptionParser.new{|opts|
  opts.on('--motif-url URL', 'Use PFM file located at some URL'){|url| options[:motif_url] = url }

  opts.on('--peaks-url URL', 'Use peaks file located at some URL'){|url| options[:peaks_url] = url }
  opts.on('--assembly-name NAME', 'Choose assembly by name') {|value| options[:assembly_name] = value }

  opts.on('--pfm', 'Force use of PFM matrix'){ options[:motif_format] = :pfm }
  opts.on('--pcm', 'Force use of PCM matrix'){ options[:motif_format] = :pcm }

  opts.on('--narrowPeak', 'Peaks are formatted in narrowPeak (peaks are reshaped into constant-size peaks around summit of a peak)'){
    options[:peaks_format] = :narrowPeak
  }
  opts.on('--bed', 'Peaks are formatted in bed (peaks are reshaped into constant-size peaks around center of a peak)'){
    options[:peaks_format] = :bed
  }

  opts.on('--gz', 'Force un-gzipping peaks'){ options[:peaks_compression] = :gz }
  opts.on('--not-compressed', 'Prevent un-gzipping peaks'){ options[:peaks_compression] = false }

  opts.on('--curve-points', 'ROC curve points'){ options[:curve_points] = true }
#   opts.on('--plot', 'Plot ROC curve'){ options[:plot_image] = true }
#   opts.on('--plot-filename', "Specify plot filename [default=#{options[:image_filename]}]") {|filename|
#     arg[:image_filename] = filename
#   }
#   opts.on('--roc', 'Store ROC curve point'){ options[:store_roc] = true }
#   opts.on('--roc-filename', "Specify ROC curve points filename [default=#{options[:roc_filename]}]"){|filename|
#     options[:roc_filename] = filename
#   }
#   opts.on('--json', 'Print results as a json file'){ options[:jsonify_results] = true }
#   opts.on('--top', "Number of top peaks to take [default=#{options[:num_top_peaks]}]"){ options[:num_top_peaks] = Integer(value) }
}

option_parser.parse!(ARGV)

assembly_infos = obtain_and_preprocess_assembly!(options)
obtain_and_preprocess_motif!(options, necessary_motif_type: :pwm, default_motif_type: :no_default)
obtain_and_preprocess_peak_sequences!(options, assembly_infos)

motif_length = read_matrix('/workdir/motif.pwm', num_columns: 4)[:matrix].length
ape_class = 'ru.autosome.ape.PrecalculateThresholds'
sarus_class = 'ru.autosome.SARUS'

system("java -cp /app/ape.jar #{ape_class} /workdir/motif.pwm --single-motif --background uniform > /workdir/motif.thr")

auc_opts = []
auc_opts << '--curve-points'  if options[:curve_points]
auc_opts = auc_opts.join(' ')
system("java -cp /app/sarus.jar #{sarus_class}  /workdir/positive.fa  /workdir/motif.pwm  besthit " + 
    " --output-scoring-mode pvalue --pvalues-file /workdir/motif.thr  --add-flanks" +
    " | ruby /app/calculate_auc.rb #{motif_length} - #{auc_opts} ")
