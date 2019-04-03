#!/usr/bin/env Rscript
library('rjson')
read_json <- function(filename) {
  return(fromJSON(readChar(filename, file.info(filename)$size)))
}

print_matrix <- function(matrix) {
  for(row in matrix) {
    writeLines(paste(row, collapse='\t'))
  }
}

args = commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  filename = args[1]
  motif_type = args[2]
  json <- read_json(filename)
  if (motif_type == 'pcm') {
    if (!is.null(json$pcm)) {
      print_matrix(json$pcm)
    } else {
      stop("JSON has no `pcm` key")
    }
  } else if (motif_type == 'pwm') {
    if (!is.null(json$pwm)) {
      print_matrix(json$pwm)
    } else {
      stop("JSON has no `pwm` key")
    }
  } else if (motif_type == 'pfm' || motif_type == 'ppm') {
    if (!is.null(json$pfm)) {
      print_matrix(json$pfm)
    } else if (!is.null(json$ppm)) {
      print_matrix(json$ppm)
    } else {
      stop("JSON has no `pfm` or `ppm` key")
    }
  } else {
    stop("Unknown motif format")
  }
} else {
  stop("Specify json-file")
}
