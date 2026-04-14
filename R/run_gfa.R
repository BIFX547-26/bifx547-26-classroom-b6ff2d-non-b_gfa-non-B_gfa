#' Run the GFA (Genomic Feature Analyzer) program
#'
#' Executes the compiled GFA tool on a FASTA input file to identify non-B DNA motifs
#' and other genomic features.
#'
#' @param seq_file Character string specifying the path to a FASTA-formatted input file.
#'   The file must exist and be readable. Sequences must be in standard FASTA format
#'   (header lines starting with ">", followed by sequence lines).
#'
#' @param out_prefix Character string specifying the prefix for output files.
#'   Default is "gfa_output". Output files will be named as:
#'   \itemize{
#'     \item `{out_prefix}_APR.tsv` - A-Phased Repeats
#'     \item `{out_prefix}_DR.tsv` - Direct Repeats
#'     \item `{out_prefix}_GQ.tsv` - G-Quadruplexes
#'     \item `{out_prefix}_IR.tsv` - Inverted Repeats
#'     \item `{out_prefix}_MR.tsv` - Mirror Repeats
#'     \item `{out_prefix}_STR.tsv` - Short Tandem Repeats
#'     \item `{out_prefix}_Z.tsv` - Z-DNA
#'   }
#'
#' @return Invisibly returns the output from the system call. The function is called
#'   primarily for its side effect of creating output files in the current working directory.
#'
#' @details
#'   The function verifies that:
#'   \enumerate{
#'     \item The input FASTA file exists and is readable
#'     \item The compiled GFA executable is available at `src/gfa`
#'   }
#'
#'   The GFA program is called with the arguments:
#'   `-skipWGET -seq <seq_file> -out <out_prefix>`
#'
#'   The command is executed in the current working directory, and output files are
#'   created in the same location.
#'
#' @examples
#' \dontrun{
#'   # Run GFA on a sample FASTA file
#'   run_gfa("sample.fasta", "my_analysis")
#'
#'   # Read the resulting output files
#'   results <- read_gfa_results("my_analysis")
#'   str(results)
#' }
#'
#' @seealso [read_gfa_results()] for reading and parsing the output TSV files.
#'
#' @keywords genomics dna-motifs
#' @export
run_gfa <- function(seq_file, out_prefix = "gfa_output") {
  if (!file.exists(seq_file)) {
    stop("Input FASTA file does not exist.", call. = FALSE)
  }

  if (!file.exists("src/gfa")) {
    stop("Compiled gfa executable not found in src. Build it first.", call. = FALSE)
  }

  result <- system2(
    "src/gfa",
    args = c("-skipWGET", "-seq", seq_file, "-out", out_prefix),
    stdout = TRUE,
    stderr = TRUE
  )

  invisible(result)
}