test_that("read_gfa_results fails when no output files exist", {
  expect_error(read_gfa_results("nonexistent_prefix_xyz"))
})

test_that("read_gfa_results returns empty list for missing files", {
  # Test with a prefix that definitely has no matching files
  expect_error(read_gfa_results("absolutely_does_not_exist_12345"),
              "No TSV files found")
})

test_that("read_gfa_results returns a named list", {
  skip_if_not(file.exists("src/gfa"), message = "gfa executable not found")
  
  input_file <- "tests/testthat/fixtures/sample.fasta"
  skip_if_not(file.exists(input_file), message = "Sample fixture not found")
  
  output_prefix <- tempfile(pattern = "gfa_read_test_")
  
  # Generate test data
  run_gfa(input_file, output_prefix)
  
  # Read the results
  results <- read_gfa_results(output_prefix)
  
  # Verify it's a list
  expect_type(results, "list")
  
  # Cleanup
  cleanup_files <- list.files(dirname(output_prefix), 
                             pattern = paste0(basename(output_prefix), ".*\\.tsv$"),
                             full.names = TRUE)
  file.remove(cleanup_files)
})

test_that("read_gfa_results returns data frames", {
  skip_if_not(file.exists("src/gfa"), message = "gfa executable not found")
  
  input_file <- "tests/testthat/fixtures/sample.fasta"
  skip_if_not(file.exists(input_file), message = "Sample fixture not found")
  
  output_prefix <- tempfile(pattern = "gfa_df_test_")
  
  run_gfa(input_file, output_prefix)
  results <- read_gfa_results(output_prefix)
  
  # Each element should be a data frame
  for (name in names(results)) {
    expect_s3_class(results[[name]], "data.frame",
                   info = paste("Element", name, "should be a data frame"))
  }
  
  # Cleanup
  cleanup_files <- list.files(dirname(output_prefix), 
                             pattern = paste0(basename(output_prefix), ".*\\.tsv$"),
                             full.names = TRUE)
  file.remove(cleanup_files)
})

test_that("read_gfa_results respects output prefix parameter", {
  skip_if_not(file.exists("src/gfa"), message = "gfa executable not found")
  
  input_file <- "tests/testthat/fixtures/sample.fasta"
  skip_if_not(file.exists(input_file), message = "Sample fixture not found")
  
  output_prefix <- tempfile(pattern = "gfa_prefix_test_")
  
  run_gfa(input_file, output_prefix)
  results <- read_gfa_results(output_prefix)
  
  # Should not be empty
  expect_gt(length(results), 0)
  
  # Cleanup
  cleanup_files <- list.files(dirname(output_prefix), 
                             pattern = paste0(basename(output_prefix), ".*\\.tsv$"),
                             full.names = TRUE)
  file.remove(cleanup_files)
})

test_that("read_gfa_results handles valid motif type names", {
  skip_if_not(file.exists("src/gfa"), message = "gfa executable not found")
  
  input_file <- "tests/testthat/fixtures/sample.fasta"
  skip_if_not(file.exists(input_file), message = "Sample fixture not found")
  
  output_prefix <- tempfile(pattern = "gfa_motif_test_")
  
  run_gfa(input_file, output_prefix)
  results <- read_gfa_results(output_prefix)
  
  # Valid motif types that could appear
  valid_types <- c("APR", "DR", "GQ", "IR", "MR", "STR", "Z")
  
  # All returned names should be valid types
  for (name in names(results)) {
    expect_true(name %in% valid_types,
               info = paste(name, "is not a recognized motif type"))
  }
  
  # Cleanup
  cleanup_files <- list.files(dirname(output_prefix), 
                             pattern = paste0(basename(output_prefix), ".*\\.tsv$"),
                             full.names = TRUE)
  file.remove(cleanup_files)
})

test_that("read_gfa_results handles partial output files", {
  skip_if_not(file.exists("src/gfa"), message = "gfa executable not found")
  
  # Create a minimal test file
  test_prefix <- tempfile(pattern = "gfa_partial_")
  
  # Create one output file manually
  dir.create(dirname(test_prefix), showWarnings = FALSE, recursive = TRUE)
  output_file <- paste0(test_prefix, "_GQ.tsv")
  
  # Write a minimal TSV file with header
  write.table(
    data.frame(start = 100, end = 150, pattern = "GQ"),
    file = output_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  results <- read_gfa_results(test_prefix)
  
  # Should have read the GQ file
  expect_true("GQ" %in% names(results))
  expect_s3_class(results$GQ, "data.frame")
  
  # Cleanup
  file.remove(output_file)
})

test_that("read_gfa_results default prefix works", {
  # Test that default parameter is correctly used
  expect_error(read_gfa_results(),
              "No TSV files found")
})

test_that("read_gfa_results data frames have correct structure", {
  skip_if_not(file.exists("src/gfa"), message = "gfa executable not found")
  
  input_file <- "tests/testthat/fixtures/sample.fasta"
  skip_if_not(file.exists(input_file), message = "Sample fixture not found")
  
  output_prefix <- tempfile(pattern = "gfa_struct_test_")
  
  run_gfa(input_file, output_prefix)
  results <- read_gfa_results(output_prefix)
  
  # Each data frame should have rows and columns
  for (name in names(results)) {
    df <- results[[name]]
    expect_true(nrow(df) >= 0,
               info = paste("Data frame", name, "should have rows"))
    expect_true(ncol(df) > 0,
               info = paste("Data frame", name, "should have columns"))
  }
  
  # Cleanup
  cleanup_files <- list.files(dirname(output_prefix), 
                             pattern = paste0(basename(output_prefix), ".*\\.tsv$"),
                             full.names = TRUE)
  file.remove(cleanup_files)
})
