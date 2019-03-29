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
    simpleError("Incorrect peak format")
  }
}

find_mounted_peaks_file <- function() {
  bed_files = c('/peaks.bed', '/peaks.bed.gz')
  narrowPeak_files = c('/peaks.narrowPeak', '/peaks.narrowPeak.gz')
  no_format_files = c('/peaks', '/peaks.gz')
  acceptable_peak_files = c(no_format_files, bed_files, narrowPeak_files)
  existing_peak_files = file.exists(acceptable_peak_files)

  if (sum(existing_peak_files) == 0) {
    simpleError("Provide a file with peaks. Either mount to /peak or its counterparts, or pass it via URL.")
  } else if (sum(existing_peak_files) > 1) {
    simpleError("Provide the only file with peaks.")
  }

  peak_filename = acceptable_peak_files[existing_peak_files][1]
  return(peak_filename)
}

obtain_and_preprocess_peaks <- function(opts) {
  if (!is.na(opts$peaks_url)) {
    peak_filename = download_file(opts$peaks_url)
  } else {
    peak_filename = find_mounted_peaks_file()
  }

  peaks_format_info = refine_peaks_format_guess(guess_peak_format(peak_filename), opts)

  peak_filename = decompress_file(peak_filename, peaks_format_info$compression)
  peak_filename = get_single_points_bed(peak_filename, peaks_format_info$peak_format)
  file.copy(peak_filename, "/workdir/peak_centers_scored.bed")
}
