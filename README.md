# nonBgfa

R package refactor of the legacy non-B_gfa codebase.

This package reorganizes the original project into an R package structure with:

- R functions in `R/`
- C source code in `src/`
- tests in `tests/`

## Functions

### `run_gfa()`

Runs the compiled `gfa` program on a FASTA input file.

Example:

```r
run_gfa("tests/testthat/fixtures/sample.fasta", "sample_run")
results <- read_gfa_results("sample_run")
head(results$STR)