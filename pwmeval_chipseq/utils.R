download_file <- function(url) {
  # We want to find out original name of the downloaded file.
  # So we create a folder and download the only file to it
  # and can get the name of that file
  dirname = tempfile()
  dir.create(dirname)
  system(paste("wget -P", shQuote(dirname), shQuote(url)))
  original_fn = list.files(dirname)[1]
  return(file.path(dirname, original_fn))
}
