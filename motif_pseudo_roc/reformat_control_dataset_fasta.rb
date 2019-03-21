# Takes sequences from fasta files and rename them so that each sequence is named in the form 
# {dataset name}:{FASTA header}:{sequence length}
# Sequence length plays an important role in correcting p-values
dataset_fn = ARGV[0]
dataset_name = File.basename(dataset_fn, File.basename(dataset_fn))

File.readlines(dataset_fn).slice_before{|l|
  l.start_with?('>')
}.each{|hdr, *seqs|
  seq = seqs.map(&:strip).join
  hdr = hdr.chomp[1..-1]
  puts ">#{dataset_name}:#{hdr}:#{seq.length}"
  puts seq
}
