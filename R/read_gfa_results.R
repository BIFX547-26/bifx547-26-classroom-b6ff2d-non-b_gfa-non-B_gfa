#' Read GFA Analysis Results
#'
#' Reads and parses TSV output files produced by the GFA program for a given analysis run.
#'
#' @param out_prefix Character string specifying the prefix used when running [run_gfa()].
#'   Default is "gfa_output". The function will look for TSV files matching the pattern
#'   `{out_prefix}_{motif_type}.tsv`.
#'
#' @return A named list of data frames, where each name corresponds to a motif type
#'   and the data frame contains the results for that motif type. Possible names include:
#'   \itemize{
#'     \item `APR` - A-Phased Repeats: repeating pattern of A and T bases
#'     \item `DR` - Direct Repeats: identical sequences in direct orientation
#'     \item `GQ` - G-Quadruplexes: G4 structures with four guanine planes
#'     \item `IR` - Inverted Repeats: sequences followed by reverse complement
#'     \item `MR` - Mirror Repeats: palindromic repeat patterns
#'     \item `STR` - Short Tandem Repeats: microsatellites
#'     \item `Z` - Z-DNA: left-handed DNA conformation
#'   }
#'   Only data frames for motif types with existing output files are included.
#'   Each data frame contains tab-separated values with a header row.
#'
#' @details
#'   The function reads TSV files produced by [run_gfa()] and returns them as a named list.
#'   If the GFA program did not find any instances of a particular motif type,
#'   the corresponding TSV file may not exist, and it will not be included in the result.
#'
#'   The function requires at least one TSV file to exist for the given `out_prefix`.
#'   If no TSV files are found, an error is raised.
#'
#' @examples
#' \dontrun{
#'   # Run GFA analysis
#'   run_gfa("sequences.fasta", "my_run")
#'
#'   # Read the results
#'   results <- read_gfa_results("my_run")
#'
#'   # Examine the structure
#'   str(results)
#'   names(results)
#'
#'   # Access individual motif types
#'   apr_data <- results$APR
#'   gq_data <- results$GQ
#'
#'   # Check which motifs were found
#'   if (length(results) > 0) {
#'     cat("Found", length(results), "motif types\n")
#'   }
#' }
#'
#' @seealso [run_gfa()] for running the GFA analysis program.
#'
#' @keywords genomics dna-motifs io
#' @export
read_gfa_results <- function(out_prefix = "gfa_output") {
  suffixes <- c("APR", "DR", "GQ", "IR", "MR", "STR", "Z")
  files <- paste0(out_prefix, "_", suffixes, ".tsv")
  names(files) <- suffixes

  existing <- files[file.exists(files)]

  if (length(existing) == 0) {
    stop("No TSV files found for this output prefix.", call. = FALSE)
  }

  # Read files, filtering out empty ones
  results <- list()
  for (name in names(existing)) {
    file_path <- existing[[name]]
    file_size <- file.size(file_path)
    
    # Only read files that have content beyond header
    if (file_size > 0) {
      tryCatch({
        df <- read.delim(file_path, sep = "\t", header = TRUE)
        # Only include if there are actual data rows
        if (nrow(df) > 0) {
          results[[name]] <- df
        }
      }, error = function(e) {
        # Skip files that can't be read
        NULL
      })
    }
  }
  
  results
}