obtain_and_preprocess_assembly <- function(opts) {
  if (is.na(opts$assembly_name)) {
    if (dir.exists("/assembly") && file.exists("/assembly/assembly.fa")) {
      assembly_fasta_fn = "/assembly/assembly.fa"
      assembly_sizes_fn = "/assembly/assembly.chrom.sizes"
    } else if (file.exists("/assembly.fa")) {
      assembly_fasta_fn = "/assembly.fa"
      assembly_sizes_fn = "/assembly.chrom.sizes"
    } else {
      stop("Error! Specify assembly name or mount assembly files (preferably via /assembly folder)")
    }
  } else {
    if (dir.exists("/assembly")) {
      assembly_fasta_fn = file.path("/assembly", paste0(opts$assembly_name, ".fa"))
      assembly_sizes_fn = file.path("/assembly", paste0(opts$assembly_name, ".chrom.sizes"))
    } else {
      assembly_fasta_fn = "/assembly.fa"
      assembly_sizes_fn = "/assembly.chrom.sizes"
    }
  }

  if (!file.exists(assembly_fasta_fn)) {
    system(paste("/app/download_assembly.sh", opts$assembly_name, shQuote(assembly_fasta_fn)))
  }

  if (!file.exists(assembly_sizes_fn)) {
    writeLines(paste("Chromosome sizes file", assembly_sizes_fn, "not found, generating..."), con=stderr())
    system(paste("/app/chrom_sizes", shQuote(assembly_fasta_fn), " > ", shQuote(assembly_sizes_fn)))
  }
  return(list(fasta_fn=assembly_fasta_fn, sizes_fn=assembly_sizes_fn))
}
