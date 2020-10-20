require_relative 'utils'

def obtain_and_preprocess_assembly!(opts, necessary: false)
  if opts[:assembly_name]
    if Dir.exist?("/assembly")
      assembly_fasta_fn = "/assembly/#{opts[:assembly_name]}.fa"
      assembly_sizes_fn = "/assembly/#{opts[:assembly_name]}.chrom.sizes"
    else
      assembly_fasta_fn = "/assembly.fa"
      assembly_sizes_fn = "/assembly.chrom.sizes"
    end
  else
    if Dir.exist?("/assembly") && File.exist?("/assembly/assembly.fa")
      assembly_fasta_fn = "/assembly/assembly.fa"
      assembly_sizes_fn = "/assembly/assembly.chrom.sizes"
    elsif File.exist?("/assembly.fa")
      assembly_fasta_fn = "/assembly.fa"
      assembly_sizes_fn = "/assembly.chrom.sizes"
    else
      if necessary
        raise "Error! Specify assembly name or mount assembly files (preferably via /assembly folder)."
      else
        return nil
      end
    end
  end

  if !File.exist?(assembly_fasta_fn) && File.exist?("#{assembly_fasta_fn}.gz")
    decompress_file("#{assembly_fasta_fn}.gz", :gz)
  end

  system("/app/download_assembly.sh #{opts[:assembly_name]} #{assembly_fasta_fn.shellescape}")  if !File.exist?(assembly_fasta_fn)

  if !File.exist?(assembly_sizes_fn)
    $stderr.puts "Chromosome sizes file #{assembly_sizes_fn} not found, generating..."
    system("/app/chrom_sizes #{assembly_fasta_fn.shellescape} > #{assembly_sizes_fn.shellescape}")
  end
  {fasta_fn: assembly_fasta_fn, chromosome_sizes_fn: assembly_sizes_fn}
end
