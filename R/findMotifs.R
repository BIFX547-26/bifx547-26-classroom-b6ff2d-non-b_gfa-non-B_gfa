#' Find non-B DNA forming motifs
#'
#' @param sequence A DNA sequence as a character string.
#' @param ... Optional parameters to override defaults.
#' @return A list of DataFrames containing found motifs.
#' @export
findMotifs <- function(sequence, ...) {
  params <- list(...)
  find_motifs_cpp(sequence, params)
}
