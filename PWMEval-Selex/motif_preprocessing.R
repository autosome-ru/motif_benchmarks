source("/app/utils.R")
source("/app/pcm2pfm_utils.R")

obtain_and_preprocess_motif <- function(opts) {
  if (is.na(opts$motif_url) && is.na(opts$motif_fn)) {
    stop("Specify motif file or motif URL.")
  }
  else if (!is.na(opts$motif_url) && !is.na(opts$motif_fn)) {
    stop("You should specify either motif file or motif URL, but not both.")
  } else if (!is.na(opts$motif_url)) {
    motif_filename = download_file(opts$motif_url)
  } else {
    motif_filename = opts$motif_fn
  }

  motif_format = refine_motif_format_guess(guess_motif_format(motif_filename), opts)
  pfm_motif_filename = get_pfm(motif_filename, motif_format)
  return(pfm_motif_filename)
}

guess_motif_format <- function(filename) {
  if (endsWith(filename, ".pcm")) {
    seq_format = "pcm"
  } else if (endsWith(filename, ".pfm") || endsWith(filename, ".ppm")) {
    seq_format = "pfm"
  } else { # default
    seq_format = "pcm"
  }
  return(seq_format)
}

# override motif format based on command-line options
refine_motif_format_guess <- function(guessed_format, opts) {
  format = guessed_format
  if (opts$pfm) {
    format = 'pfm'
  } else if (opts$pcm) {
    format = 'pcm'
  }
  return(format)
}

get_pfm <- function(filename, format) {
  if (format == 'pcm') {
    tmp_fn = tempfile()
    pcm2pfm_files(filename, tmp_fn)
    return(tmp_fn)
  } else if (format == 'pfm') {
    return(filename)
  } else {
    stop("Unknown motif format")
  }
}
