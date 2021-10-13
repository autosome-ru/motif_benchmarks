str_is_numeric <- function(str) {
  return(!is.na(suppressWarnings(as.numeric(str))))
}

is_matrix_row <- function(str) {
  terms <- unlist(strsplit(str, '\\s+'))
  return(all(str_is_numeric(terms)))
}

parse_matrix <- function(lines) {
  return(lapply(strsplit(lines, '\\s+'), as.numeric))
}

print_matrix <- function(matrix, con=stdout()) {
  for(row in matrix) {
    writeLines(paste(row, collapse='\t'), con=con)
  }
}

pcm2pfm <- function(pcm) {
  return(lapply(pcm, function(row){row / sum(row)}))
}

read_lines_file_or_stdin <- function(input_filename) {
  is_stdin = is.na(input_filename) || (input_filename == "stdin") || (input_filename == "-")
  if (is_stdin) {
    con = file("stdin")
  } else {
    con = file(input_filename, "r")
  }
  lines = suppressWarnings(readLines(con))
  if (!is_stdin) {
    close(con)
  }
  return(lines)
}

pcm2pfm_files <- function(input_filename, output_filename) {
  lines = read_lines_file_or_stdin(input_filename)
  has_header = !is_matrix_row(lines[1])
  is_stdout = is.na(output_filename) || (output_filename == "stdout") || (output_filename == "-")
  if (is_stdout) {
    output_con = stdout()
  } else {
    output_con = file(output_filename, "w")
  }
  if (has_header) {
    writeLines(lines[1], con=output_con)
    pcm <- parse_matrix(lines[2:length(lines)])
  } else {
    pcm <- parse_matrix(lines)
  }

  pfm <- pcm2pfm(pcm)
  print_matrix(pfm, con=output_con)
  if (!is_stdout) {
    close(output_con)
  }
}
