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

def decompress_file(filename, compression)
  case compression
  when false
    filename
  when :gz
    $stderr.puts "decompress #{filename} (compression: #{compression})"
    tmp_file = Tempfile.new.tap(&:close)
    system("gzip -cd #{filename.shellescape} > #{tmp_file.path.shellescape}")
    tmp_file.path
  else
    raise 'Unknown compression format'
  end
end
