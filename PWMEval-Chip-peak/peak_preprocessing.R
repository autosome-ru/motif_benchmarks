source("/app/utils.R")

guess_peak_format <- function(filename) {
  if (endsWith(filename, ".gz")) {
    compression = "gz"
  } else {
    compression = "no"
  }

  filename = sub("\\.gz$", "", filename)

  if (endsWith(filename, ".bed")) {
    peak_format = "bed"
  } else if (endsWith(filename, ".narrowPeak")) {
    peak_format = "narrowPeak"
  } else { # default
    peak_format = "bed"
  }
  return(list(compression=compression, peak_format=peak_format))
}

# override sequences format/compression based on command-line options
refine_peaks_format_guess <- function(guessed_format, opts) {
  peak_format = guessed_format$peak_format
  compression = guessed_format$compression

  if (opts$peak_format_bed) {
    peak_format = 'bed'
  }
  if (opts$peak_format_narrowPeak) {
    peak_format = 'narrowPeak'
  }
  if (opts$compression_no_peaks) {
    compression = 'no'
  }
  if (opts$compression_gz_peaks) {
    compression = 'gz'
  }
  return(list(compression=compression, peak_format=peak_format))
}

narrowPeak_summit <- function(peak_filename) {
  tmp_fn = tempfile()
  system(paste("/app/narrowpeak_summit.sh", shQuote(peak_filename), " > ", shQuote(tmp_fn)))
  return(tmp_fn)
}

peak_center <- function(peak_filename) {
  tmp_fn = tempfile()
  system(paste("/app/peak_center.sh", shQuote(peak_filename), " > ", shQuote(tmp_fn)))
  return(tmp_fn)
}

get_single_points_bed <- function(peak_filename, peak_format) {
  if (peak_format == 'bed') {
    return(peak_center(peak_filename))
  } else if (peak_format == 'narrowPeak') {
    return(narrowPeak_summit(peak_filename))
  } else {
    stop("Incorrect peak format")
  }
}

obtain_and_preprocess_peak_centers <- function(opts) {
  if (is.na(opts$peaks_url) && is.na(opts$peaks_fn)) {
    stop("Specify peaks file or peaks URL.")
  } else if (!is.na(opts$peaks_url) && !is.na(opts$peaks_fn)) {
    stop("You should specify either peaks file or peaks URL, but not both.")
  } else if (!is.na(opts$peaks_url)) {
    peak_filename = download_file(opts$peaks_url)
  } else {
    peak_filename = opts$peaks_fn
  }

  peaks_format_info = refine_peaks_format_guess(guess_peak_format(peak_filename), opts)

  peak_filename = decompress_file(peak_filename, peaks_format_info$compression)
  peak_centers_filename = get_single_points_bed(peak_filename, peaks_format_info$peak_format)
  return(peak_centers_filename)
}

get_top_peaks <- function(peaks_filename, num_top_peaks) {
  tmp_fn = tempfile()
  system(paste("sort -k5,5nr ", peaks_filename, " | head -n", num_top_peaks, " > ", shQuote(tmp_fn)))
  return(tmp_fn)
}

extend_positive_peaks <- function(peak_centers_filename, assembly) {
  tmp_fn = tempfile()
  system(paste("/app/bedtools slop -i ", peak_centers_filename, " -g ", shQuote(assembly$sizes_fn), " -l 124 -r 125  > ", shQuote(tmp_fn)))
  return(tmp_fn)
}

extend_negative_peaks <- function(peak_centers_filename, assembly) {
  tmp_fn = tempfile()
  system(paste("/app/bedtools slop -i ", peak_centers_filename, " -g ", shQuote(assembly$sizes_fn), " -l -301 -r 550  > ", shQuote(tmp_fn)))
  system(paste("/app/bedtools slop -i ", peak_centers_filename, " -g ", shQuote(assembly$sizes_fn), " -l 549 -r -300  >> ", shQuote(tmp_fn)))
  return(tmp_fn)
}

get_peaks_fasta <- function(peaks_filename, assembly) {
  tmp_fn = tempfile()
  system(paste("/app/bedtools getfasta -bed ", peaks_filename, " -fi ", shQuote(assembly$fasta_fn), "  > ", shQuote(tmp_fn)))
  return(tmp_fn)
}
