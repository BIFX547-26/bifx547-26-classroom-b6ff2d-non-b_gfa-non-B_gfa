test_that("run_gfa fails for missing FASTA file", {
  expect_error(run_gfa("fake_file.fasta"))
})

test_that("read_gfa_results fails when no output files exist", {
  expect_error(read_gfa_results("not_real_prefix"))
})