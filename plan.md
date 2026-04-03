# Refactoring Plan: GFA Suite to R Package

## Objective
Refactor the original GFA (Genomic Feature Analyzer) C codebase into a modular library and create a functional R package (`nonBgfa`) while maintaining full functional parity with the original command-line tool.

## Completed Tasks

### 1. Research and Baseline
- Analyzed existing C source files and `Makefile`.
- Established a baseline using the provided `gfa_test.fasta` and benchmark output files.
- Verified the original tool's functionality.

### 2. Architectural Refactoring
- **Centralized Globals:** Moved large global arrays (DNA buffers and repeat results) into `src/globals.c` to manage state consistently across the library and CLI.
- **Modular API:** Created `src/gfa_api.h` and `src/gfa_api.c` to encapsulate the motif-finding logic.
    - `run_gfa_core`: The core engine that processes pre-populated buffers.
    - `run_gfa_analysis`: A high-level entry point for string-based sequences.
- **Header Compatibility:** Updated `src/gfa.h` with `extern "C"` guards and resolved naming conflicts with R headers (specifically `TRUE`/`FALSE` definitions).

### 3. R Package Development (`nonBgfa`)
- **Structure:** Set up standard R package directories: `R/`, `src/`, `DESCRIPTION`, `NAMESPACE`.
- **Rcpp Integration:** Implemented `src/Rcpp_gfa.cpp` to bridge C++ and R, converting internal `REP` structures into R DataFrames.
- **R Wrapper:** Created `R/findMotifs.R` to provide a user-friendly interface with support for parameter overrides.

### 4. CLI Tool Update
- Refactored `src/gfa_cli.c` to use the new `GFA_Params` and `run_gfa_core` API.
- Updated the `Makefile` to support the new modular source structure.

### 5. Verification and Validation
- **R Package Test:** Created a script to compare `nonBgfa::findMotifs()` output against benchmark counts. Results matched perfectly (e.g., 14 IRs, 7 GQs).
- **CLI Parity Test:** Verified that the refactored `./gfa` binary produces TSV/GFF files identical to the original version.

## Revised Memory Strategy
To improve scalability and prevent excessive static memory usage in the R package:
- **Transition to Dynamic Allocation:** Replace fixed-size static arrays (e.g., `dna[300000001]`) with dynamic pointers.
- **Lazy Initialization:** Allocate memory only when a sequence is being analyzed, based on the actual sequence length.
- **Memory Safety:** Ensure all dynamically allocated buffers are properly freed after analysis to prevent leaks, especially when called repeatedly from R.
- **Error Handling:** Implement checks for allocation failures and provide meaningful errors to the R environment.

## Project Structure
- `src/`: Core C library, Rcpp wrappers, and CLI source.
- `R/`: R package function definitions.
- `DESCRIPTION` / `NAMESPACE`: R package metadata.
- `Makefile`: Build instructions for the CLI tool.
- `README.md`: Updated with R package usage and installation instructions.

## Final Status
The project is now a dual-purpose repository:
1. A high-performance command-line tool for batch processing FASTA files.
2. A modern R package for interactive genomic analysis.
