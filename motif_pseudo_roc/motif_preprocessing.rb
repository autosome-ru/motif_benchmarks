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
  if opts[:motif_fn]
    if !opts[:motif_url]
      motif_filename = opts[:motif_fn]
    else
      raise "Specify only one of motif URL and motif filename"
    end
  else # !opts[:motif_fn]
    if opts[:motif_url]
      motif_filename = download_file(opts[:motif_url])
    else
      raise "Specify one of motif URL and motif filename"
    end
  end

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
  motif_filename
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
  pfm_matrix = pcm[:matrix].map{|row|
    norm = row.sum
    row.map{|x| x.to_f / norm }
  }
  {name: pcm[:name], matrix: pfm_matrix}
end

def pfm2pcm(pfm, word_count: 100)
  pcm_matrix = pfm[:matrix].map{|row|
    row.map{|el| el * word_count }
  }
  {name: pfm[:name], matrix: pcm_matrix}
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
      File.write(tempfile, matrix_as_string(pfm))
    }
  when :pfm
    filename
  when :pwm
    raise "Can't convert pwm --> pfm"
  else
    raise "Unknown motif format #{format}"
  end
end

def get_pcm(filename, format, opts)
  case format
  when :pcm
    filename
  when :pfm
    pfm = read_matrix(filename, num_columns: 4)
    pcm = pfm2pcm(pfm, word_count: opts[:word_count])
    tempname{|tempfile|
      File.write(tempfile, matrix_as_string(pcm))
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
      File.write(tempfile, matrix_as_string(pwm))
    }
  when :pfm
    pfm = read_matrix(filename, num_columns: 4)
    pcm = pfm2pcm(pfm, word_count: opts[:word_count])
    pwm = pcm2pwm(pcm, pseudocount: opts[:pseudocount])
    tempname{|tempfile|
      File.write(tempfile, matrix_as_string(pwm))
    }
  when :pwm
    filename
  else
    raise "Unknown motif format #{format}"
  end
end

def get_motif_thresholds(motif_fn, background_type:, background:)
  thresholds_file = register_new_tempfile("motif.thr").tap(&:close)
  thresholds_fn = thresholds_file.path
  if background_type == :mono
    ape_class = 'ru.autosome.ape.PrecalculateThresholds'
    system("java -cp /app/ape.jar #{ape_class} #{motif_fn} --single-motif --background #{background} > #{thresholds_fn}")
  elsif background_type == :di
    ape_class = 'ru.autosome.ape.di.PrecalculateThresholds'
    system("java -cp /app/ape.jar #{ape_class} #{motif_fn} --single-motif --background #{background} --from-mono > #{thresholds_fn}")
  else
    raise "Should not be here"
  end
  thresholds_fn
end
