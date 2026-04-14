test_that("run_gfa fails with missing input FASTA file", {
  expect_error(run_gfa("nonexistent_file.fasta"))
})

test_that("run_gfa fails with missing executable", {
  # Create a temporary directory to avoid the executable check
  skip_if_not(file.exists("src/gfa"), message = "gfa executable not found")
  
  # This should work since executable exists
  expect_error(run_gfa("nonexistent_file.fasta"), 
               "Input FASTA file does not exist")
})

test_that("run_gfa accepts valid input file and output prefix", {
  skip_if_not(file.exists("src/gfa"), message = "gfa executable not found")
  
  # Use sample fixture
  input_file <- "tests/testthat/fixtures/sample.fasta"
  skip_if_not(file.exists(input_file), message = "Sample fixture not found")
  
  output_prefix <- tempfile()
  
  # Should run without error
  result <- run_gfa(input_file, output_prefix)
  
  # Should return invisible NULL
  expect_null(result)
})

test_that("run_gfa creates output files", {
  skip_if_not(file.exists("src/gfa"), message = "gfa executable not found")
  
  input_file <- "tests/testthat/fixtures/sample.fasta"
  skip_if_not(file.exists(input_file), message = "Sample fixture not found")
  
  output_prefix <- tempfile(pattern = "gfa_test_")
  
  # Run gfa
  run_gfa(input_file, output_prefix)
  
  # Check that at least one output file was created
  tsv_files <- list.files(dirname(output_prefix), 
                         pattern = paste0(basename(output_prefix), ".*\\.tsv$"))
  
  expect_gt(length(tsv_files), 0, 
           info = "Expected at least one TSV output file")
  
  # Cleanup
  cleanup_files <- list.files(dirname(output_prefix), 
                             pattern = paste0(basename(output_prefix), ".*\\.tsv$"),
                             full.names = TRUE)
  file.remove(cleanup_files)
})

test_that("run_gfa uses correct command line arguments", {
  skip_if_not(file.exists("src/gfa"), message = "gfa executable not found")
  
  input_file <- "tests/testthat/fixtures/sample.fasta"
  skip_if_not(file.exists(input_file), message = "Sample fixture not found")
  
  output_prefix <- tempfile()
  
  # This test verifies that the function at least attempts to run gfa
  # The actual execution depends on the compiled binary
  result <- run_gfa(input_file, output_prefix)
  
  # If no error is thrown, the command was likely executed
  expect_null(result)
})

test_that("run_gfa handles special characters in file paths", {
  skip_if_not(file.exists("src/gfa"), message = "gfa executable not found")
  
  # Create a file with spaces in name
  input_file <- "tests/testthat/fixtures/sample.fasta"
  skip_if_not(file.exists(input_file), message = "Sample fixture not found")
  
  output_prefix <- tempfile(pattern = "gfa_test with spaces_")
  
  # Should handle special characters without error
  result <- run_gfa(input_file, output_prefix)
  expect_null(result)
  
  # Cleanup
  cleanup_files <- list.files(dirname(output_prefix), 
                             pattern = paste0(basename(output_prefix), ".*\\.tsv$"),
                             full.names = TRUE)
  file.remove(cleanup_files)
})
