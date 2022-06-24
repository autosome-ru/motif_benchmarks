require 'json'
require 'optparse'
require_relative 'cli_args'

require_relative 'utils'
require_relative 'motif_preprocessing'
require_relative 'peak_preprocessing'
require_relative 'assembly_preprocessing'

require_relative 'background'
require_relative 'sequence_dataset'

options = {
  # num_top_peaks: 500,
  flank_size: 150,
  pseudocount: :log,
  word_count: 100,
  curve_points: false,
  summit_column: 10,
  background_type: :di, background: :infer,
}

option_parser = OptionParser.new{|opts|
  optparse_add_motif_opts(opts, options)
  optparse_add_peaks_opts(opts, options)
  opts.on('--background VALUE', "One can specify mono/dinucleotide background inferred from (di)nucleotide frequencies in dataset or can use custom background. " + 
                                "Available values: `infer:mono`, `infer:di`, `mono:<pA>,<pC>,<pG>,<pT>`, `di:<pAA>,<pAC>,...,<pTT>`, " +
                                "`gc:<GC content in [0,1] range>` or `uniform`. (Default: `infer:di`)"){|value|
    configure_background!(options, value)
  }

  opts.on('--curve-points', 'ROC curve points'){ options[:curve_points] = true }
  # opts.on('--top', "Number of top peaks to take [default=#{options[:num_top_peaks]}]"){ options[:num_top_peaks] = Integer(value) }
}

option_parser.parse!(ARGV)

# We don't need assembly when peaks are provided in FASTA format, thus we can ignore
# missing `--assembly-name` or not provided `/assembly.fa`.
# An exception will be triggered later if assembly is needed but absent.
assembly_infos = obtain_and_preprocess_assembly!(options, necessary: false)

# negative fasta is not used in this benchmark and not provided!
positive_fasta_fn = obtain_and_preprocess_peak_sequences!(options, assembly_infos)

motif_fn = obtain_and_preprocess_motif!(options, necessary_motif_type: :pwm, default_motif_type: :no_default)

if options[:background] == :infer
  dataset = SequenceDataset.new(positive_fasta_fn)
  if options[:background_type] == :mono
    background = local_mono_background(dataset).join(',')
  elsif options[:background_type] == :di
    background = local_di_background(dataset).join(',')
  else
    raise "Should not be here"
  end
else
  background = options[:background]
end

thresholds_fn = get_motif_thresholds(motif_fn, background_type: options[:background_type], background: background)

sarus_class = 'ru.autosome.SARUS'
motif_length = read_matrix(motif_fn, num_columns: 4)[:matrix].length

auc_opts = []
auc_opts << '--curve-points'  if options[:curve_points]
auc_opts = auc_opts.join(' ')
system("java -cp /app/sarus.jar #{sarus_class}  #{positive_fasta_fn}  #{motif_fn}  besthit " +
    " --output-scoring-mode pvalue --pvalues-file #{thresholds_fn}  --add-flanks" +
    " | ruby /app/calculate_auc.rb #{motif_length} - #{auc_opts} ")
