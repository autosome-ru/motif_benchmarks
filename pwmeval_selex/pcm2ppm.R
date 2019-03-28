#!/usr/bin/env Rscript

str_is_numeric <- function(str) {
  return(!is.na(suppressWarnings(as.numeric(str))))
}

is_matrix_row <- function(str) {
  return(all(str_is_numeric(strsplit(str, '\\s+'))))
}

parse_matrix <- function(lines) {
  return(lapply(strsplit(lines, '\\s+'), as.numeric))
}

print_matrix <- function(matrix) {
  for(row in matrix) {
    writeLines(paste(row, collapse='\t'))
  }
}

pcm2ppm <- function(pcm) {
  return(lapply(pcm, function(row){row / sum(row)}))
}

args = commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  filepath = args[1]
  con = file(filepath, "r")
} else {
  con = file("stdin")
  open(con)
}
lines = suppressWarnings(readLines(con))
close(con)

has_header = !is_matrix_row(lines[1])
if (has_header) {
  writeLines(lines[1])
  pcm <- parse_matrix(lines[2:length(lines)])
} else {
  pcm <- parse_matrix(lines)
}

ppm <- pcm2ppm(pcm)
print_matrix(ppm)
