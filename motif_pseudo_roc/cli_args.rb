require 'optparse'

def configure_peaks_format!(options, format)
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
end

def configure_background!(options, value)
  return configure_background!(options, File.read(value.sub(/^file:/, ''))) if value.start_with?('file:')

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
end

def configure_num_top_peaks!(options, value)
  if value.downcase == 'all'
    options[:top_peaks] = {num_peaks: 'all'}
  elsif value.match?(/^\d+$/) # take original order (by line number)
    options[:top_peaks] = {num_peaks: Integer(value)}
  elsif value.match?(/^\d+:by:\d+:(max|min)$/) # top X peaks by Y column, take min/max values
    top_peaks_config = value.match(/^(?<amount>\d+):by:(?<column>\d+):(?<order>(max|min))$/).named_captures
    options[:top_peaks] = {
      num_peaks: Integer(top_peaks_config['amount']),
      order_by_column: Integer(top_peaks_config['column']),
      order: top_peaks_config['order'],
    }
  else
    raise "Unknown value for number of top peaks"
  end
end

def optparse_add_motif_opts(opts, options)
  opts.on('--motif FILE', 'Motif PFM file'){|filename| options[:motif_fn] = filename }
  opts.on('--motif-url URL', 'Use PFM file located at some URL'){|url| options[:motif_url] = url }
  opts.on('--pfm', 'Force use of PFM matrix'){ options[:motif_format] = :pfm }
  opts.on('--pcm', 'Force use of PCM matrix'){ options[:motif_format] = :pcm }  
end

def optparse_add_peaks_opts(opts, options)
  opts.on('--peaks FILE', 'Peaks file'){|filename| options[:peaks_fn] = filename }
  opts.on('--peaks-url URL', 'Use peaks file located at some URL'){|url| options[:peaks_url] = url }
  opts.on('--assembly-name NAME', 'Choose assembly by name') {|value| options[:assembly_name] = value }
  opts.on('--top VALUE', "Number of top peaks to take. Possible values: `all` or `<number>` or `<number>:by:<column>:<max|min>` [default=#{options[:num_top_peaks]}]"){|value|
    configure_num_top_peaks!(options, value)
  }

  opts.on('--narrowPeak', 'Peaks are formatted in narrowPeak (peaks are reshaped into constant-size peaks around summit of a peak)'){
    options[:peaks_format] = :narrowPeak
  }
  opts.on('--bed', 'Peaks are formatted in bed (peaks are reshaped into constant-size peaks around center of a peak)'){
    options[:peaks_format] = :bed
  }
  opts.on('--fasta', 'Peaks are provided as sequences in fasta format') {
    options[:peaks_format] = :fasta
  }
  opts.on('--peak-flank-size VALUE', 'In center/summit peak preprocessing modes, peaks are transformed: ' +
                                     'peak center/summit point is taken and extended with flanks of fixed size in both directions. ' +
                                     'Default flank size is 150nt.') {|value|
    opts[:flank_size] = Float(value)
  }
  opts.on('--peak-format FORMAT', 'Peaks are formatted in a custom format. ' +
                                  'FORMAT is a `<chr column>,<start column>,<end column>,<mode>` string. ' +
                                  'Column indiced are 1-based. Mode can be either `entire`, or `center`, or ' +
                                  '`summit:(abs|rel):<summit column>` for different modes of interval clipping. ' +
                                  'Summit types `abs`, `rel` are for absolute summit coordinates vs relative (from start) summit position. ' +
                                  'See also --peak-flank-size option.'){|format|
    configure_peaks_format!(options, format)
  }

  opts.on('--gz', 'Force un-gzipping peaks'){ options[:peaks_compression] = :gz }
  opts.on('--not-compressed', 'Prevent un-gzipping peaks'){ options[:peaks_compression] = false }
end
