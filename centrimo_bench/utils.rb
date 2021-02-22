require 'tempfile'
require 'tmpdir'
require 'shellwords'
require 'fileutils'

def tempname(prefix: '')
  path = nil
  while !path || File.exist?(path)
    t = Time.now.strftime("%Y%m%d")
    path = File.join(Dir.tmpdir, "#{prefix}#{t}-#{$$}-#{rand(0x100000000).to_s(36)}")
  end
  yield path  if block_given?
  path
end

$tempfile_list = []
def register_new_tempfile(*args, **kwargs)
  raise "Shouldn't pass block to a tempfile creation"  if block_given?
  file = Tempfile.new(*args, **kwargs)
  $tempfile_list << file
  file
end

def download_file(url)
  # We want to find out original name of the downloaded file.
  # So we create a folder and download the only file to it
  # and can get the name of that file
  $stderr.puts "download #{url}"
  dirname = Dir.mktmpdir('downloads')
  system("wget -P #{dirname.shellescape} #{url.shellescape}")
  Dir.glob("#{dirname}/*").first
end

def decompress_file(filename, compression, output_filename: nil)
  case compression
  when false
    filename
  when :gz
    tmp_dir = Dir.mktmpdir('decompress_file')
    basename = File.basename(filename, File.extname(filename))
    output_filename = File.join(tmp_dir, basename)  unless output_filename
    $stderr.puts "decompress #{filename} (compression: #{compression}) into #{output_filename}"
    system("gzip -cd #{filename.shellescape} > #{output_filename.shellescape}")
    output_filename
  else
    raise 'Unknown compression format'
  end
end
