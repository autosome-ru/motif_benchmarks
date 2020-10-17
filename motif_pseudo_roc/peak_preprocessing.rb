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
  {peaks_compression: peaks_compression, peaks_format: peaks_format}.tap{|result|
    $stderr.puts "For #{filename} guessed peaks format: #{result}"
  }
end

# override sequences format/compression based on command-line options
def refine_peaks_format_guess(guessed_format, opts)
  {
    peaks_compression: opts.fetch(:peaks_compression, guessed_format[:peaks_compression]),
    peaks_format: opts.fetch(:peaks_format, guessed_format[:peaks_format]),
  }.tap{|result|
    $stderr.puts "Refined peaks format: #{result}"
  }
end

def narrowPeak_summit(peaks_filename)
  $stderr.print "narrowPeak summit for #{peaks_filename}"
  tmp_file = Tempfile.new('narrowPeak_summit.bed')
  File.open(peaks_filename){|f|
    f.each_line.each{|line|
      next  if line.start_with?('#')
      row = line.chomp.split
      chr = row[0]
      start = Integer(row[1])
      peak_offset = Integer(row[9])
      peak_summit = start + peak_offset
      value = row[6]
      tmp_file.puts [chr, peak_summit, peak_summit + 1, '.', value].join("\t")
    }
  }
  tmp_file.close
  tmp_file.path.tap{|result|
    $stderr.puts " --> #{result}"
  }
end

def peak_center(peaks_filename)
  $stderr.print "peak centers for #{peaks_filename}"
  tmp_file = Tempfile.new('peak_center.bed')
  File.open(peaks_filename){|f|
    f.each_line.each{|line|
      next  if line.start_with?('#')
      row = line.chomp.split
      chr = row[0]
      start = Integer(row[1])
      stop = Integer(row[2])
      center = (start + stop) / 2
      value = row[4]
      tmp_file.puts [chr, center, center + 1, '.', value].join("\t")
    }
  }
  tmp_file.close
  tmp_file.path.tap{|result|
    $stderr.puts " --> #{result}"
  }
end

def get_single_points_bed(peaks_filename, peaks_format)
  $stderr.puts "extract peak centers #{peaks_filename} (format: #{peaks_format})"
  if peaks_format == :bed
    peak_center(peaks_filename)
  elsif peaks_format == :narrowPeak
    narrowPeak_summit(peaks_filename)
  else
    raise 'Incorrect peak format.'
  end
end

def find_mounted_peaks_file
  fasta_files = ['/peaks.fa', '/peaks.fa.gz', '/peaks.mfa', '/peaks.mfa.gz']
  bed_files = ['/peaks.bed', '/peaks.bed.gz']
  narrowPeak_files = ['/peaks.narrowPeak', '/peaks.narrowPeak.gz']
  no_format_files = ['/peaks', '/peaks.gz']
  acceptable_peak_files = [fasta_files, no_format_files, bed_files, narrowPeak_files].flatten
  existing_peak_files = acceptable_peak_files.select{|fn| File.exist?(fn) }

  if existing_peak_files.empty?
    raise 'Provide a file with peaks. Either mount to /peak or its counterparts, or pass it via URL.'
  elsif existing_peak_files.size > 1
    raise 'Provide the only file with peaks.'
  end
  existing_peak_files.first.tap{|result|
    $stderr.puts "Found peaks file #{result}"
  }
end

def slop_peaks(peaks_filename, assembly_sizes_fn, flank_size)
  $stderr.print "slop peaks in #{peaks_filename} to a certain length"
  tmp_file = Tempfile.new('slop_peaks.bed').tap(&:close)
  system("/app/bedtools slop -i #{peaks_filename.shellescape} -g #{assembly_sizes_fn.shellescape} -l #{flank_size - 1} -r #{flank_size} > #{tmp_file.path.shellescape}")
  tmp_file.path.tap{|result|
    $stderr.puts " --> #{result}"
  }
end

def fasta_by_bed_peaks(peaks_bed_filename, assembly_fn)
  $stderr.print "extract fasta from bed file #{peaks_bed_filename}"
  tmp_file = Tempfile.new('peaks.fa').tap(&:close)
  system("/app/bedtools getfasta -bed #{peaks_bed_filename.shellescape} -fi #{assembly_fn.shellescape} > #{tmp_file.path.shellescape}")
  tmp_file.path.tap{|result|
    $stderr.puts " --> #{result}"
  }
end

# control should be formatted FASTA (i.e. with seq length specified in header line)
def fasta_with_lengths(fasta_fn)
  $stderr.print "add sequence lengths to FASTA headers in #{fasta_fn}"
  tmp_file = Tempfile.new('peaks_formatted.fa')
  File.open(fasta_fn){|f|
    f.each_line.slice_before(/^>/).each{|hdr, *lines|
      seq = lines.map(&:strip).join
      tmp_file.puts "#{hdr.strip}:#{seq.size}"
      tmp_file.puts seq
    }
  }
  tmp_file.close
  tmp_file.path.tap{|result|
    $stderr.puts " --> #{result}"
  }
end

def obtain_and_preprocess_peak_sequences!(opts, assembly_infos)
  peaks_filename = opts[:peaks_url] ? download_file(opts[:peaks_url]) : find_mounted_peaks_file
  peaks_format_info = refine_peaks_format_guess(guess_peaks_format(peaks_filename), opts)

  case peaks_format_info[:peaks_format]
  when :fasta
    # do nothing
  when :bed, :narrowPeak
    $stderr.puts 'preprocess peaks'
    peaks_filename = decompress_file(peaks_filename, peaks_format_info[:peaks_compression])
    peaks_filename = get_single_points_bed(peaks_filename, peaks_format_info[:peaks_format])
    peaks_filename = slop_peaks(peaks_filename, assembly_infos[:chromosome_sizes_fn], opts[:flank_size])
    peaks_filename = fasta_by_bed_peaks(peaks_filename, assembly_infos[:fasta_fn])
  end

  peaks_filename = fasta_with_lengths(peaks_filename)
  FileUtils.cp(peaks_filename, "/workdir/positive.fa")
end
