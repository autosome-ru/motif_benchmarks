source("/app/utils.R")

guess_seq_format <- function(filename) {
  if (endsWith(filename, ".gz")) {
    compression = "gz"
  } else {
    compression = "no"
  }

  filename = sub("\\.gz$", "", filename)

  if (endsWith(filename, ".fa") || endsWith(filename, ".fasta")) {
    seq_format = "fasta"
  } else if (endsWith(filename, ".fq") || endsWith(filename, ".fastq")) {
    seq_format = "fastq"
  } else { # default
    seq_format = "fasta"
  }
  return(list(compression=compression, seq_format=seq_format))
}

# override sequences format/compression based on command-line options
refine_seq_format_guess <- function(guessed_format, opts) {
  seq_format = guessed_format$seq_format
  compression = guessed_format$compression

  if (opts$seq_format_fasta) {
    seq_format = 'fasta'
  }
  if (opts$seq_format_fastq) {
    seq_format = 'fastq'
  }
  if (opts$compression_no) {
    compression = 'no'
  }
  if (opts$compression_gz) {
    compression = 'gz'
  }
  return(list(compression=compression, seq_format=seq_format))
}

fastq2fasta <- function(seq_filename) {
  tmp_fn = tempfile()
  system(paste("/app/fastq2fasta.sh", shQuote(seq_filename), " > ", shQuote(tmp_fn)))
  return(tmp_fn)
}

convert2fasta <- function(seq_filename, seq_format) {
  if (seq_format == 'fasta') {
    return(seq_filename)
  } else if (seq_format == 'fastq') {
    return(fastq2fasta(seq_filename))
  } else {
    simpleError("Incorrect format")
  }
}

# in addition to filtering it also joins FASTA spreaded on multiple lines into single-string format
filter_fasta <- function(seq_filename, opts) {
  only_acgt = 'yes'
  if (opts$allow_iupac) {
    only_acgt = 'no'
  }

  seq_length = 'no'
  if (!is.na(opts$seq_length)) {
    seq_length = opts$seq_length
  }

  tmp_fn = tempfile()
  system(paste("/app/filter_fasta", shQuote(seq_filename), seq_length,  only_acgt, " > ", shQuote(tmp_fn)))
  return(tmp_fn)
}

filter_redundant <- function(seq_filename, opts) {
  if (opts$non_redundant) {
    tmp_fn = tempfile()
    system(paste("/app/filter_redundant.sh", shQuote(seq_filename), " > ", shQuote(tmp_fn)))
    return(tmp_fn)
  } else {
    return(seq_filename)
  }
}

find_mounted_sequences_file <- function() {
  fasta_files = c('/seq.fasta', '/seq.fa', '/seq.fasta.gz', '/seq.fa.gz')
  fastq_files = c('/seq.fastq', '/seq.fq', '/seq.fastq.gz', '/seq.fq.gz')
  no_format_files = c('/seq', '/seq.gz')
  acceptable_seq_files = c(no_format_files, fasta_files, fastq_files)
  existing_seq_files = file.exists(acceptable_seq_files)

  if (sum(existing_seq_files) == 0) {
    simpleError("Provide a file with SELEX sequences. Either mount to /seq or its counterparts, or pass it via URL.")
  } else if (sum(existing_seq_files) > 1) {
    simpleError("Provide the only file with SELEX sequences.")
  }

  seq_filename = acceptable_seq_files[existing_seq_files][1]
  return(seq_filename)
}


obtain_and_preprocess_sequences <- function(opts) {
  if (!is.na(opts$seq_url)) {
    seq_filename = download_file(opts$seq_url)
  } else {
    seq_filename = find_mounted_sequences_file()
  }

  seq_format_info = refine_seq_format_guess(guess_seq_format(seq_filename), opts)

  # process sequences file into uncompressed FASTA file
  seq_filename = decompress_file(seq_filename, seq_format_info$compression)
  seq_filename = convert2fasta(seq_filename, seq_format_info$seq_format)
  seq_filename = filter_fasta(seq_filename, opts)
  seq_filename = filter_redundant(seq_filename, opts)
  file.copy(seq_filename, "/workdir/positive.fa")
}
