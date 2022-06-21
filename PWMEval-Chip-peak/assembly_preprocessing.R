obtain_and_preprocess_assembly <- function(opts) {
  if (!is.na(opts$assembly_name)) {
    assembly_fasta_fn = file.path("/assembly", paste0(opts$assembly_name, ".fa"))
    assembly_sizes_fn = file.path("/assembly", paste0(opts$assembly_name, ".chrom.sizes"))
    if (!file.exists(assembly_fasta_fn)) {
      dir.create(dirname(assembly_fasta_fn), recursive=TRUE, showWarnings=FALSE)
      system(paste("/app/download_assembly.sh", opts$assembly_name, shQuote(assembly_fasta_fn)))
    }
  } else {
    if (is.na(opts$assembly_fasta_fn)) {
      stop("Error! Specify assembly name or mount assembly files (preferably via /assembly folder)")
    }
    assembly_fasta_fn = opts$assembly_fasta_fn
    if (!is.na(opts$assembly_sizes_fn)) {
      assembly_sizes_fn = opts$assembly_sizes_fn
    } else {
      assembly_sizes_fn = sub(".fa$", ".chrom.sizes", assembly_fasta_fn)
    }
  }

  if (!file.exists(assembly_sizes_fn)) {
    dir.create(dirname(assembly_sizes_fn), recursive=TRUE, showWarnings=FALSE)
    writeLines(paste("Chromosome sizes file", assembly_sizes_fn, "not found, generating..."), con=stderr())
    system(paste("/app/chrom_sizes", shQuote(assembly_fasta_fn), " > ", shQuote(assembly_sizes_fn)))
  }
  return(list(fasta_fn=assembly_fasta_fn, sizes_fn=assembly_sizes_fn))
}
