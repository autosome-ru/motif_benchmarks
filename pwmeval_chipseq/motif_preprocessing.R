source("/app/utils.R")

obtain_and_preprocess_motif <- function(opts) {
  if (!is.na(opts$motif_url)) {
    motif_filename = download_file(opts$motif_url)
  } else {
    motif_filename = find_mounted_motif_file()
  }

  motif_format = refine_motif_format_guess(guess_motif_format(motif_filename), opts)
  motif_filename = get_ppm(motif_filename, motif_format)
  file.copy(motif_filename, "/workdir/motif.ppm")
}

find_mounted_motif_file <- function() {
  ppm_files = c('/motif.ppm', '/motif.pfm', '/matrix.ppm', '/matrix.pfm')
  pcm_files = c('/motif.pcm', '/matrix.pcm')
  no_format_files = c('/motif', '/matrix')
  acceptable_motif_files = c(no_format_files, ppm_files, pcm_files)
  existing_motif_files = file.exists(acceptable_motif_files)

  if (sum(existing_motif_files) == 0) {
    simpleError("Provide a file with positional frequencies/counts matrix. Either mount to /motif or its counterparts, or pass it via URL.")
  } else if (sum(existing_motif_files) > 1) {
    simpleError("Provide the only file with positional frequencies/counts matrix")
  }

  motif_filename = acceptable_motif_files[existing_motif_files][1]
  return(motif_filename)
}

guess_motif_format <- function(filename) {
  if (endsWith(filename, ".pcm")) {
    seq_format = "pcm"
  } else if (endsWith(filename, ".pfm") || endsWith(filename, ".ppm")) {
    seq_format = "ppm"
  } else { # default
    seq_format = "pcm"
  }
  return(seq_format)
}

# override motif format based on command-line options
refine_motif_format_guess <- function(guessed_format, opts) {
  format = guessed_format
  if (opts$ppm) {
    format = 'ppm'
  } else if (opts$pcm) {
    format = 'pcm'
  }
  return(format)
}

get_ppm <- function(filename, format) {
  if (format == 'pcm') {
    tmp_fn = tempfile()
    system(paste("/app/pcm2ppm.R", shQuote(filename), " > ", shQuote(tmp_fn)))
    return(tmp_fn)
  } else if (format == 'ppm') {
    return(filename)
  } else {
    simpleError("Unknown motif format")
  }
}
