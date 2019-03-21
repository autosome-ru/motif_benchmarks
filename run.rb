require 'json'
require 'shellwords'

def read_matrix(fn, num_columns: 4)
  lines = File.readlines(fn).map(&:chomp)
  rows = lines.map{|l| l.split("\t") }
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

config_fn = "/workdir/config.json"
result_fn = "/workdir/persistent/result.json"
config = JSON.parse(File.read(config_fn))

control_fn = '/benchmark_specific_data/control.formatted.mfa' # control should be formatted FASTA (i.e. with seq length specified in header line)
pwm_fn = 'motif.pwm'
File.write(pwm_fn, config['motif'])

motif_length = read_matrix(pwm_fn, num_columns: 4)[:matrix].length
ape_class = 'ru.autosome.ape.PrecalculateThresholds'
sarus_class = 'ru.autosome.SARUS'

thresholds_fn = 'motif.thr'

system("java -cp /app/ape.jar #{ape_class} #{pwm_fn.shellescape} --single-motif --background uniform > #{thresholds_fn.shellescape}")

system("java -cp /app/sarus.jar #{sarus_class} #{control_fn.shellescape} #{pwm_fn.shellescape} besthit " + 
    " --output-scoring-mode pvalue --pvalues-file #{thresholds_fn.shellescape} --add-flanks" +
    " | ruby /app/calculate_auc.rb #{motif_length} - " +
    " > #{result_fn.shellescape}")
