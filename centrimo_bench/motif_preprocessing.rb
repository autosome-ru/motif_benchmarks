require_relative 'utils'

def read_matrix(fn, num_columns: 4)
  lines = File.readlines(fn).map(&:chomp)
  rows = lines.map{|l| l.split }
  name = File.basename(fn, File.extname(fn))
  unless (rows[0].size == num_columns) && rows[0].all?{|x| Float(x, exception: false) }
    hdr = lines.first
    rows.shift
    name = hdr.start_with?('>') ? hdr[1..-1].strip : hdr.strip
  end
  matrix = rows.map{|row|
    row.map{|x| Float(x) }
  }
  raise  if matrix.empty?
  raise  unless matrix.all?{|row| row.size == num_columns }
  {name: name, matrix: matrix}
end

def obtain_and_preprocess_motif!(opts, necessary_motif_type:, default_motif_type:)
  motif_filename = opts[:motif_url] ? download_file(opts[:motif_url]) : find_mounted_motif_file

  motif_format = refine_motif_format_guess(guess_motif_format(motif_filename, default_motif_type: default_motif_type), opts)
  case necessary_motif_type
  when :pcm
    motif_filename = get_pcm(motif_filename, motif_format, opts)
  when :pfm
    motif_filename = get_pfm(motif_filename, motif_format, opts)
  when :pwm
    motif_filename = get_pwm(motif_filename, motif_format, opts)
  else
    raise "Unknown motif format #{necessary_motif_type}"
  end
  FileUtils.cp(motif_filename, "/workdir/motif.#{necessary_motif_type}")
end

def find_mounted_motif_file
  pfm_files = ['/motif.pfm', '/motif.ppm', '/matrix.pfm', '/matrix.ppm']
  pcm_files = ['/motif.pcm', '/matrix.pcm']
  no_format_files = ['/motif', '/matrix']
  acceptable_motif_files = [no_format_files, pfm_files, pcm_files].flatten
  existing_motif_files = acceptable_motif_files.select{|fn| File.exist?(fn) }

  if existing_motif_files.size == 0
    raise 'Provide a file with positional frequencies/counts matrix. Either mount to /motif or its counterparts, or pass it via URL.'
  elsif existing_motif_files.size > 1
    raise 'Provide the only file with positional frequencies/counts matrix'
  end
  existing_motif_files.first
end

def guess_motif_format(filename, default_motif_type:)
  case File.extname(filename)
  when '.pwm'
    :pwm
  when '.pcm'
    :pcm
  when '.pfm', '.ppm'
    :pfm
  else
    default_motif_type
  end
end

# override motif format based on command-line options
def refine_motif_format_guess(guessed_format, opts)
  opts.fetch(:motif_format, guessed_format)
end

#############################

def matrix_as_string(model)
  res = [">#{model[:name]}"]
  res += model[:matrix].map{|row| row.join("\t") }
  res.join("\n")
end

def matrix_as_meme_string(model)
  res = [">#{model[:name]}"]
  res += [
    "MEME version 4",
    "",
    "ALPHABET= ACGT",
    "",
    "MOTIF #{model[:name]}",
    "letter-probability matrix: alength= 4 w= #{model[:matrix].length} nsites= #{model.fetch(:word_count, 100)} E= 0"
  ]
  res += model[:matrix].map{|row| row.join("\t") }
  res.join("\n")
end
#############################

def calculate_pseudocount(count, pseudocount: :log)
  case pseudocount
  when :log
    Math.log([count, 2].max);
  when :sqrt
    Math.sqrt(count)
  else Numeric
    pseudocount
  end
end

#############################

def pcm2pfm(pcm)
  word_count = pcm[:word_count] || pcm[:matrix].map(&:sum).max
  pfm_matrix = pcm[:matrix].map{|row|
    norm = row.sum
    row.map{|x| x.to_f / norm }
  }
  {name: pcm[:name], matrix: pfm_matrix, word_count: word_count}
end

def pfm2pcm(pfm, word_count: )
  word_count = word_count || pfm[:word_count] || 100
  pcm_matrix = pfm[:matrix].map{|row|
    row.map{|el| el * word_count }
  }
  {name: pfm[:name], matrix: pcm_matrix, word_count: word_count}
end

def pcm2pwm(pcm, pseudocount: :log)
  pwm_matrix = pcm[:matrix].map{|row|
    count = row.sum
    row.map{|el|
      pseudocount_value = calculate_pseudocount(count, pseudocount: pseudocount)
      numerator = el + 0.25 * pseudocount_value
      denominator = 0.25 * (count + pseudocount_value)
      Math.log(numerator / denominator)
    }
  }
  {name: pcm[:name], matrix: pwm_matrix}
end

#############################

def get_pfm(filename, format, opts)
  case format
  when :pcm
    pcm = read_matrix(filename, num_columns: 4)
    pfm = pcm2pfm(pcm)
    tempname{|tempfile|
      File.write(tempfile, matrix_as_meme_string(pfm))
    }
  when :pfm
    pfm = read_matrix(filename, num_columns: 4)
    tempname{|tempfile|
      File.write(tempfile, matrix_as_meme_string(pfm))
    }
  when :pwm
    raise "Can't convert pwm --> pfm"
  else
    raise "Unknown motif format #{format}"
  end
end

def get_pcm(filename, format, opts)
  case format
  when :pcm
    pcm = read_matrix(filename, num_columns: 4)
    tempname{|tempfile|
      File.write(tempfile, matrix_as_meme_string(pcm))
    }
  when :pfm
    pfm = read_matrix(filename, num_columns: 4)
    pcm = pfm2pcm(pfm, word_count: opts[:word_count])
    tempname{|tempfile|
      File.write(tempfile, matrix_as_meme_string(pcm))
    }
  when :pwm
    raise "Can't convert pwm --> pcm"
  else
    raise "Unknown motif format #{format}"
  end
end

def get_pwm(filename, format, opts)
  case format
  when :pcm
    pcm = read_matrix(filename, num_columns: 4)
    pwm = pcm2pwm(pcm, pseudocount: opts[:pseudocount])
    tempname{|tempfile|
      File.write(tempfile, matrix_as_meme_string(pwm))
    }
  when :pfm
    pfm = read_matrix(filename, num_columns: 4)
    pcm = pfm2pcm(pfm, word_count: opts[:word_count])
    pwm = pcm2pwm(pcm, pseudocount: opts[:pseudocount])
    tempname{|tempfile|
      File.write(tempfile, matrix_as_meme_string(pwm))
    }
  when :pwm
    pwm = read_matrix(filename, num_columns: 4)
    tempname{|tempfile|
      File.write(tempfile, matrix_as_meme_string(pwm))
    }
  else
    raise "Unknown motif format #{format}"
  end
end
