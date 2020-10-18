require 'json'
require 'optparse'

require_relative 'utils'
require_relative 'motif_preprocessing'
require_relative 'peak_preprocessing'
require_relative 'assembly_preprocessing'
# require 'bioinform'

require_relative 'sequence_dataset'

options = {
#   image_filename: "/results/roc_curve.png",
#   roc_filename: "roc_curve.tsv",
#   num_top_peaks: 500,
#   plot_curve: false,
  flank_size: 150,
  pseudocount: :log,
  word_count: 100,
  curve_points: false,
  summit_column: 10,
  background_type: :di, background: :infer,
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
  opts.on('--peak-format FORMAT', 'Peaks are formatted in a custom format. ' +
                                  'FORMAT is a `<chr column>,<start column>,<end column>,<value column><mode>` string. ' +
                                  'Column indiced are 1-based. Mode can be either `entire`, or `center`, or ' +
                                  '`summit:(abs|rel):<summit column>` for different modes of interval clipping. ' +
                                  'Summit types `abs`, `rel` are for absolute summit coordinates vs relative (from start) summit position'){|format|
  chr_column, start_column, end_column, value_column, mode = format.split(',')
    options[:peaks_format] = :custom
    options[:peaks_format_config] = {chr_column: Integer(chr_column), start_column: Integer(start_column), end_column: Integer(end_column), value_column: Integer(value_column)}
    if mode == 'entire'
      options[:peaks_format_config][:mode] = :entire
    elsif mode == 'center'
      options[:peaks_format_config][:mode] = :center
    elsif mode.split(':')[0] == 'summit'
      options[:peaks_format_config][:mode] = :summit
      options[:peaks_format_config][:summit_column] = Integer(mode.split(':')[2])

      summit_type = mode.split(':')[1]
      if summit_type == 'abs'
        options[:peaks_format_config][:summit_type] = :absolute
      elsif summit_type == 'rel'
        options[:peaks_format_config][:summit_type] = :relative
      else
        raise "Unknown summit type `#{summit_type}`. Should be abs or rel"
      end
    else
      raise 'Unknown mode of format'
    end
  }

  opts.on('--gz', 'Force un-gzipping peaks'){ options[:peaks_compression] = :gz }
  opts.on('--not-compressed', 'Prevent un-gzipping peaks'){ options[:peaks_compression] = false }

  opts.on('--background VALUE', "One can specify mono/dinucleotide background inferred from (di)nucleotide frequencies in dataset or can use custom background. " + 
                                "Available values: `infer:mono`, `infer:di`, `mono:<pA>,<pC>,<pG>,<pT>`, `di:<pAA>,<pAC>,...,<pTT>`, " +
                                "`gc:<GC content in [0,1] range>` or `uniform`. (Default: `infer:di`)"){|value|
    if value == 'uniform'
      options[:background_type] = :mono
      options[:background] = 'uniform'
    elsif value.start_with?('gc:')
      options[:background_type] = :mono
      options[:background] = Float(value.split(':')[1])
    elsif value.start_with?('mono:')
      options[:background_type] = :mono
      options[:background] = value.split(':')[1] # "A,C,G,T"
      raise  unless options[:background].split(',').map{|x| Float(x) }.size == 4
    elsif value.start_with?('di:')
      options[:background_type] = :di
      options[:background] = value.split(':')[1]  # "AA,AC,AG,AT,CA,CC,...TG,TT"
      raise  unless options[:background].split(',').map{|x| Float(x) }.size == 16
    elsif value.start_with?('infer:mono')
      options[:background_type] = :mono
      options[:background] = :infer
    elsif value.start_with?('infer:di')
      options[:background_type] = :di
      options[:background] = :infer
    else
      raise "Unknown background mode `#{value}`."
    end
  }

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

if options[:background] == :infer
  dataset = SequenceDataset.new('/workdir/positive.fa')
  if options[:background_type] == :mono
    background = dataset.local_mono_background.join(',')
  elsif options[:background_type] == :di
    background = dataset.local_di_background.join(',')
  else
    raise "Should not be here"
  end
else
  background = options[:background]
end


motif_length = read_matrix('/workdir/motif.pwm', num_columns: 4)[:matrix].length
sarus_class = 'ru.autosome.SARUS'

if options[:background_type] == :mono
  ape_class = 'ru.autosome.ape.PrecalculateThresholds'
  system("java -cp /app/ape.jar #{ape_class} /workdir/motif.pwm --single-motif --background #{background} > /workdir/motif.thr")
elsif options[:background_type] == :di
  ape_class = 'ru.autosome.ape.di.PrecalculateThresholds'
  system("java -cp /app/ape.jar #{ape_class} /workdir/motif.pwm --single-motif --background #{background} --from-mono > /workdir/motif.thr")
else
  raise "Should not be here"
end

auc_opts = []
auc_opts << '--curve-points'  if options[:curve_points]
auc_opts = auc_opts.join(' ')
system("java -cp /app/sarus.jar #{sarus_class}  /workdir/positive.fa  /workdir/motif.pwm  besthit " +
    " --output-scoring-mode pvalue --pvalues-file /workdir/motif.thr  --add-flanks" +
    " | ruby /app/calculate_auc.rb #{motif_length} - #{auc_opts} ")
