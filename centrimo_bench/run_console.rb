require 'json'
require 'optparse'

require_relative 'utils'
require_relative 'motif_preprocessing'
require_relative 'peak_preprocessing'
require_relative 'assembly_preprocessing'


# Resulting values are intentionally left as strings, not numbers. That's because some of values can be
# e.g. "1.5e-436" which can't be converted to float without loss of precision
def read_centrimo_results(filename)
  lines = File.readlines(filename).map(&:chomp)
  row = lines[1].split("\t").map(&:strip)
  db_index, motif_id, motif_alt_id, consensus, evalue, adj_pvalue, log_adj_pvalue, \
    bin_location, bin_width, total_width, sites_in_bin, total_sites, p_success, pvalue, mult_tests = *row
  {
    motif_id: motif_id, consensus: consensus,
    evalue: evalue, adj_pvalue: adj_pvalue, log_adj_pvalue: log_adj_pvalue,
    bin_location: bin_location, bin_width: bin_width, total_width: total_width, sites_in_bin: sites_in_bin,
    total_sites: total_sites, p_success: p_success, pvalue: pvalue, mult_tests: mult_tests,
  }
end

# aka `concentration`
def central_probability(sites, motif_length:, sequence_length:, total_sites:, window_size: 20)
  start = (0.5 * (sequence_length - motif_length - 1) - 0.5 * window_size).ceil
  stop = start + window_size
  sites[start .. stop].sum(0.0) / total_sites
end

options = {
  flank_size: 250,
  pseudocount: :log,
  word_count: 100,
  curve_points: false,
  summit_column: 10,
  background_type: :di, background: :infer,
  jsonify_results: false,
}

option_parser = OptionParser.new{|opts|
  opts.on('--motif-url URL', 'Use PFM file located at some URL'){|url| options[:motif_url] = url }

  opts.on('--peaks-url URL', 'Use peaks file located at some URL'){|url| options[:peaks_url] = url }
  opts.on('--assembly-name NAME', 'Choose assembly by name') {|value| options[:assembly_name] = value }

  opts.on('--pfm', 'Force use of PFM matrix'){ options[:motif_format] = :pfm }
  opts.on('--pcm', 'Force use of PCM matrix'){ options[:motif_format] = :pcm }

  opts.on('--word-count', "Specify word count (default count for PFM is #{options[:word_count]},\nfor PCM default is calculated as sum of counts)"){|value|
    options[:word_count] = Float(value)
  }
  opts.on('--pseudocount', 'Specify pseudocount (float value/`sqrt`/`log`. default: log, i.e. log(count))'){|value|
    options[:pseudocount] = Float(value) rescue value
  }

  opts.on('--narrowPeak', 'Peaks are formatted in narrowPeak (peaks are reshaped into constant-size peaks around summit of a peak)'){
    options[:peaks_format] = :narrowPeak
  }
  opts.on('--bed', 'Peaks are formatted in bed (peaks are reshaped into constant-size peaks around center of a peak)'){
    options[:peaks_format] = :bed
  }
  opts.on('--peak-flank-size VALUE', 'In center/summit peak preprocessing modes, peaks are transformed: ' +
                                     'peak center/summit point is taken and extended with flanks of fixed size in both directions. ' +
                                     "Default flank size is #{options[:flank_size]}nt.") {|value|
    opts[:flank_size] = Float(value)
  }
  opts.on('--peak-format FORMAT', 'Peaks are formatted in a custom format. ' +
                                  'FORMAT is a `<chr column>,<start column>,<end column>,<mode>` string. ' +
                                  'Column indiced are 1-based. Mode can be either `entire`, or `center`, or ' +
                                  '`summit:(abs|rel):<summit column>` for different modes of interval clipping. ' +
                                  'Summit types `abs`, `rel` are for absolute summit coordinates vs relative (from start) summit position. ' +
                                  'See also --peak-flank-size option.'){|format|
  chr_column, start_column, end_column, mode = format.split(',')
    options[:peaks_format] = :custom
    options[:peaks_format_config] = {chr_column: Integer(chr_column), start_column: Integer(start_column), end_column: Integer(end_column)}
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

  opts.on('--window-size', '...')

  opts.on('--gz', 'Force un-gzipping peaks'){ options[:peaks_compression] = :gz }
  opts.on('--not-compressed', 'Prevent un-gzipping peaks'){ options[:peaks_compression] = false }

  opts.on('--json', 'Print results as a json file'){ options[:jsonify_results] = true }
}

option_parser.parse!(ARGV)

# We don't need assembly when peaks are provided in FASTA format, thus we can ignore
# missing `--assembly-name` or not provided `/assembly.fa`.
# An exception will be triggered later if assembly is needed but absent.
assembly_infos = obtain_and_preprocess_assembly!(options, necessary: false)

obtain_and_preprocess_motif!(options, necessary_motif_type: :pfm, default_motif_type: :no_default)
obtain_and_preprocess_peak_sequences!(options, assembly_infos)

word_count = get_meme_word_count('/workdir/motif.pfm')
motif_length = get_meme_motif_length('/workdir/motif.pfm')
pseudocount = calculate_pseudocount(word_count, pseudocount: options[:pseudocount])

system("centrimo /workdir/positive.fa /workdir/motif.pfm --oc /results --verbosity 1 --motif-pseudo #{pseudocount}")
info = read_centrimo_results('/results/centrimo.tsv')
if !info[:motif_id]
  $stderr.puts "Fallback: all sequences were filtered out so we lower threshold"
  system("centrimo /workdir/positive.fa /workdir/motif.pfm --oc /results --verbosity 1 --motif-pseudo #{pseudocount} --score 0")
  info = read_centrimo_results('/results/centrimo.tsv')
end

info = read_centrimo_results('/results/centrimo.tsv')

# calculate concentration for different window sizes
sites = File.readlines('/results/site_counts.txt').drop(1).map(&:strip).map{|l| Float(l.split("\t").last) }

concentrations = (5..100).step(5).map{|window_size|
  concentration = central_probability(sites,
                                      motif_length: motif_length,
                                      sequence_length: 2 * options[:flank_size],
                                      total_sites: Integer(info[:total_sites]),
                                      window_size: window_size)
  {window_size: window_size, concentration: concentration}
}

info[:concentrations] = concentrations

if options[:jsonify_results]
  puts info.to_json
else
  puts info[:evalue]
end
