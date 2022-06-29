require_relative 'utils'

def guess_peaks_format(filename)
  peaks_compression = (File.extname(filename) == '.gz') ? :gz : false
  basename = File.basename(filename, '.gz')

  case File.extname(basename)
  when '.fa'
    peaks_format = :fasta
  when '.bed'
    peaks_format = :bed
  when '.narrowPeak'
    peaks_format = :narrowPeak
  else # default
    peaks_format = :fasta
  end
  {peaks_compression: peaks_compression, peaks_format: peaks_format}
end

# override sequences format/compression based on command-line options
def refine_peaks_format_guess(guessed_format, opts)
  {
    peaks_compression: opts.fetch(:peaks_compression, guessed_format[:peaks_compression]),
    peaks_format: opts.fetch(:peaks_format, guessed_format[:peaks_format]),
  }
end

def peak_summit(peaks_filename, peaks_format_config)
  chr_column = peaks_format_config[:chr_column] - 1
  start_column = peaks_format_config[:start_column] - 1
  summit_column = peaks_format_config[:summit_column] - 1
  summit_type = peaks_format_config[:summit_type]
  tmp_file = register_new_tempfile('narrowPeak_summit.bed')
  File.open(peaks_filename){|f|
    f.each_line.each{|line|
      next  if line.start_with?('#')
      row = line.chomp.split
      chr = row[chr_column]
      start = Integer(row[start_column])
      if summit_type == :relative
        peak_offset = Integer(row[summit_column])
        peak_summit = start + peak_offset
      elsif summit_type == :absolute
        peak_summit = Integer(row[summit_column])
      else
        raise "Unknown type of summit `#{summit_type}`. Should be :relative or :absolute."
      end
      tmp_file.puts [chr, peak_summit, peak_summit + 1,].join("\t")
    }
  }
  tmp_file.close
  tmp_file.path
end

def peak_center(peaks_filename, peaks_format_config)
  chr_column = peaks_format_config[:chr_column] - 1
  start_column = peaks_format_config[:start_column] - 1
  end_column = peaks_format_config[:end_column] - 1
  tmp_file = register_new_tempfile('peak_center.bed')
  File.open(peaks_filename){|f|
    f.each_line.each{|line|
      next  if line.start_with?('#')
      row = line.chomp.split
      chr = row[chr_column]
      start = Integer(row[start_column])
      stop = Integer(row[end_column])
      center = (start + stop) / 2
      tmp_file.puts [chr, center, center + 1].join("\t")
    }
  }
  tmp_file.close
  tmp_file.path
end

def entire_peak(peaks_filename, peaks_format_config)
  chr_column = peaks_format_config[:chr_column] - 1
  start_column = peaks_format_config[:start_column] - 1
  end_column = peaks_format_config[:end_column] - 1
  tmp_file = register_new_tempfile('entire_peak.bed')
  File.open(peaks_filename){|f|
    f.each_line.each{|line|
      next  if line.start_with?('#')
      row = line.chomp.split
      chr = row[chr_column]
      start = Integer(row[start_column])
      stop = Integer(row[end_column])
      tmp_file.puts [chr, start, stop].join("\t")
    }
  }
  tmp_file.close
  tmp_file.path
end

def slop_peaks(peaks_filename, assembly_sizes_fn, flank_size)
  tmp_file = register_new_tempfile('slop_peaks.bed').tap(&:close)
  system("/app/bedtools slop -i #{peaks_filename.shellescape} -g #{assembly_sizes_fn.shellescape} -l #{flank_size - 1} -r #{flank_size} > #{tmp_file.path.shellescape}")
  tmp_file.path
end

def fasta_by_bed_peaks(peaks_bed_filename, assembly_fn)
  tmp_file = register_new_tempfile('peaks.fa').tap(&:close)
  system("/app/bedtools getfasta -bed #{peaks_bed_filename.shellescape} -fi #{assembly_fn.shellescape} > #{tmp_file.path.shellescape}")
  tmp_file.path
end

# control should be formatted FASTA (i.e. with seq length specified in header line)
def fasta_with_lengths(fasta_fn)
  tmp_file = register_new_tempfile('positive_peaks_formatted.fa')
  File.open(fasta_fn){|f|
    f.each_line.slice_before(/^>/).each{|hdr, *lines|
      seq = lines.map(&:strip).join
      tmp_file.puts "#{hdr.strip}:#{seq.size}"
      tmp_file.puts seq
    }
  }
  tmp_file.close
  tmp_file.path
end

def infer_peaks_format_config(peaks_format, opts)
  case peaks_format
  when :custom
    opts[:peaks_format_config]
  when :bed
    {chr_column: 1, start_column: 2, end_column: 3, mode: :center}
  when :narrowPeak
    {chr_column: 1, start_column: 2, end_column: 3, mode: :summit, summit_column: 10, summit_type: :relative}
  when :fasta
    {mode: :nop}
  else
  end
end

def get_top_peaks(peaks_filename, top_peaks_config)
  num_peaks = top_peaks_config[:num_peaks]
  return peaks_filename  if num_peaks == 'all'
  num_peaks = Integer(num_peaks)

  tmp_file = register_new_tempfile("top_peaks.#{File.extname(peaks_filename)}")

  peaks = File.readlines(peaks_filename).reject{|l|
    l.start_with?('#')
  }
  if top_peaks_config[:order]
    column = Integer(top_peaks_config[:order_by_column]) - 1  # 1-based column indices
    peaks = peaks.sort_by{|l|
      row = l.split("\t")
      Float(row[column])
    }
    case top_peaks_config[:order]
    when 'min'
      peaks = peaks.first(num_peaks)
    when 'max'
      peaks = peaks.reverse.first(num_peaks)
    else
      raise "Shouldn't be here"
    end
  else # order not specified — just take first lines
    peaks = peaks.first(num_peaks)
  end

  peaks.each{|l|
    tmp_file.puts(l)
  }
  tmp_file.close
  tmp_file.path
end

def get_top_fasta(peaks_filename, top_peaks_config)
  num_peaks = top_peaks_config[:num_peaks]
  return peaks_filename  if num_peaks == 'all'
  raise "FASTA files cannot be sorted by custom field"  if top_peaks_config[:order]

  tmp_file = register_new_tempfile('top_peaks.fa')
  File.readlines(peaks_filename).slice_before{|l|
    l.start_with?(">")
  }.first(num_peaks).each{|lns|
    lns.each{|l|
      tmp_file.puts(l)
    }
  }

  tmp_file.close
  tmp_file.path
end

def obtain_and_preprocess_peak_sequences!(opts, assembly_infos)
  if opts[:peaks_fn]
    if !opts[:peaks_url]
      peaks_filename = opts[:peaks_fn]
    else
      raise "Specify only one of peaks URL and peaks filename"
    end
  else # !opts[:peaks_fn]
    if opts[:peaks_url]
      peaks_filename = download_file(opts[:peaks_url])
    else
      raise "Specify one of peaks URL and peaks filename"
    end
  end

  peaks_format_info = refine_peaks_format_guess(guess_peaks_format(peaks_filename), opts)
  peaks_format = peaks_format_info[:peaks_format] # bed/narrowPeak/custom/fasta

  peaks_format_config = infer_peaks_format_config(peaks_format, opts) # center/summit/entire or NOP config
  mode = peaks_format_config[:mode]

  peaks_filename = decompress_file(peaks_filename, peaks_format_info[:peaks_compression])

  if mode == :nop
    # do nothing (e.g. for FASTA peak format)
    peaks_filename = get_top_fasta(peaks_filename, opts[:top_peaks])
  else
    raise "Error! Specify assembly name or mount assembly files (preferably via /assembly folder)."  if !assembly_infos
    peaks_filename = get_top_peaks(peaks_filename, opts[:top_peaks])
    if mode == :center
      peaks_filename = peak_center(peaks_filename, peaks_format_config)
      peaks_filename = slop_peaks(peaks_filename, assembly_infos[:chromosome_sizes_fn], opts[:flank_size])
    elsif mode == :summit
      peaks_filename = peak_summit(peaks_filename, peaks_format_config)
      peaks_filename = slop_peaks(peaks_filename, assembly_infos[:chromosome_sizes_fn], opts[:flank_size])
    elsif mode == :entire
      peaks_filename = entire_peak(peaks_filename, peaks_format_config)
    else
      raise "Unknown mode `#{mode}`"
    end

    peaks_filename = fasta_by_bed_peaks(peaks_filename, assembly_infos[:fasta_fn])
  end

  positive_peaks_fasta_filename = fasta_with_lengths(peaks_filename)
  positive_peaks_fasta_filename
end
