# nonBgfa: Genomic Feature Analyzer R Package

An R package for identifying and analyzing non-B DNA motifs and genomic features in nucleic acid sequences. This is a refactored R wrapper around the legacy non-B_gfa codebase, organized into a modern R package structure.

## Overview

The **nonBgfa** package detects and characterizes various non-B DNA structures and sequence features, including:

- **APR** - A-Phased Repeats (alternating A/T base patterns)
- **DR** - Direct Repeats (identical sequences in direct orientation)
- **GQ** - G-Quadruplexes (four-stranded G4 structures)
- **IR** - Inverted Repeats (sequences followed by reverse complement)
- **MR** - Mirror Repeats (palindromic repeat patterns)
- **STR** - Short Tandem Repeats (microsatellites)
- **Z** - Z-DNA (left-handed DNA conformation)

## Installation

### Prerequisites

- R 4.0 or later
- A C compiler (for building from source)
- The FASTA input file containing your sequences

### Installing from source

```r
# Using devtools
devtools::install()

# Or using R CMD INSTALL
R CMD INSTALL .
```

## Package Structure

The package is organized as follows:

```
nonBgfa/
├── R/                    # R functions
│   ├── run_gfa.R
│   └── read_gfa_results.R
├── src/                  # C source code
│   ├── gfa.c            # Main entry point
│   ├── gfa.h            # Header with structure definitions
│   ├── findAPR.c        # A-Phased Repeat detection
│   ├── findDR.c         # Direct Repeat detection
│   ├── findGQ.c         # G-Quadruplex detection
│   ├── findIR.c         # Inverted Repeat detection
│   ├── findMR.c         # Mirror Repeat detection
│   ├── findSTR.c        # Short Tandem Repeat detection
│   ├── findZDNA.c       # Z-DNA detection
│   ├── read_fasta.c     # FASTA file reading
│   ├── print_tsv_file.c # TSV output formatting
│   └── ...              # Additional utility functions
├── tests/               # Test suite
│   └── testthat/
│       ├── test-run_gfa.R
│       ├── test-run_gfa_comprehensive.R
│       ├── test-read_gfa_results_comprehensive.R
│       └── fixtures/    # Test data files
└── man/                 # Generated documentation
```

## Usage

### Basic Workflow

1. **Prepare your sequences** in FASTA format
2. **Run GFA analysis** with `run_gfa()`
3. **Read and analyze results** with `read_gfa_results()`

### Example Analysis

```r
library(nonBgfa)

# Step 1: Run analysis on your sequences
run_gfa("my_sequences.fasta", "my_analysis")

# Step 2: Read the results
results <- read_gfa_results("my_analysis")

# Step 3: Explore the results
names(results)  # See which motif types were found

# Access specific motif results
if ("GQ" %in% names(results)) {
  gq_data <- results$GQ
  head(gq_data)
}

# Get summary statistics
lapply(results, nrow)  # Number of features per motif type
```

## Functions

### `run_gfa(seq_file, out_prefix = "gfa_output")`

Executes the GFA analysis on a FASTA input file.

**Parameters:**
- `seq_file` (character): Path to FASTA-formatted input file
- `out_prefix` (character, default: `"gfa_output"`): Prefix for output files

**Returns:** Invisibly returns the system call output (primarily called for side effects)

**Output files created:**
- `{out_prefix}_APR.tsv` - A-Phased Repeats
- `{out_prefix}_DR.tsv` - Direct Repeats  
- `{out_prefix}_GQ.tsv` - G-Quadruplexes
- `{out_prefix}_IR.tsv` - Inverted Repeats
- `{out_prefix}_MR.tsv` - Mirror Repeats
- `{out_prefix}_STR.tsv` - Short Tandem Repeats
- `{out_prefix}_Z.tsv` - Z-DNA

**Example:**

```r
# Run analysis with default output prefix
run_gfa("gfa_test.fasta")

# Run with custom output prefix
run_gfa("sequences.fasta", "my_analysis_results")
```

### `read_gfa_results(out_prefix = "gfa_output")`

Reads and parses TSV output files produced by the GFA analysis.

**Parameters:**
- `out_prefix` (character, default: `"gfa_output"`): Prefix used when running `run_gfa()`

**Returns:** A named list of data frames, one for each detected motif type

**Example:**

```r
# Read results from run with prefix "my_analysis"
results <- read_gfa_results("my_analysis")

# See which motif types were found
names(results)

# Access results by motif type
apr_repeats <- results$APR
str_repeats <- results$STR

# Get basic statistics
nrow(results$GQ)  # Number of G-Quadruplexes found
```

## Output File Format

The GFA tool generates tab-separated values (TSV) files for each motif type found. Each file contains a header row with column names followed by one row per detected motif instance.

### Common TSV Columns

While the exact columns depend on the motif type, typical columns include:

- **start** - Starting position in the sequence
- **end** - Ending position in the sequence
- **pattern** - The detected motif pattern
- **strand** - Which strand (+ or -)
- **score** - Statistical or structural score

### Example TSV File Structure

```
start	end	pattern	strand	score
100	150	ATATATATATAT	+	85
250	310	GCGCGCGCGC	+	92
450	520	AAAAAATTTT	-	78
```

## Troubleshooting

### "Input FASTA file does not exist"

**Problem:** The path to your FASTA file is incorrect or the file doesn't exist.

**Solution:** Verify the file path:
```r
file.exists("my_file.fasta")
# Check current directory
getwd()
# List files in current directory
list.files(pattern = "*.fasta")
```

### "Compiled gfa executable not found in src"

**Problem:** The C code hasn't been compiled yet.

**Solution:** Rebuild the package:
```r
devtools::load_all()  # Compiles C code
# Or install the package properly
devtools::install()
```

### "No TSV files found for this output prefix"

**Problem:** The GFA analysis didn't produce any output files.

**Possible causes:**
- The GFA program didn't find any motifs in your sequences
- The output prefix is incorrect
- The GFA program encountered an error

**Solution:** Check that the output files were created:
```r
# Check what files were created
list.files(pattern = "gfa_output.*")
# Try running with a different input file
```

### Empty Results

**Problem:** GFA ran successfully but no motifs were found.

**This is normal** if your sequences don't contain the specific motif types you're looking for. Not all sequences will contain all motif types.

## FASTA Input Format

Your input file must be in standard FASTA format:

```
>sequence_name_1
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>sequence_name_2
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
```

Requirements:
- Header lines start with `>`
- Sequence names can contain letters, numbers, hyphens, and underscores
- Sequences use standard IUPAC nucleotide codes (A, T, G, C, N)
- No length restrictions

## Testing

The package includes comprehensive test coverage for all functions:

```r
# Run the test suite
devtools::test()

# Run tests with detailed output
testthat::test_file("tests/testthat/test-run_gfa_comprehensive.R")
```

## License

GPL-3 (see LICENSE file)

## Citation

If you use this package in your research, please cite the original GFA work:

```
nonBgfa: An R package for genomic feature analysis
Version 0.0.0.9000
```

## Support & Contributing

For bug reports, feature requests, or contributions, please open an issue on the GitHub repository.

## References

- Non-B DNA structures and their biological significance
- IUPAC nucleotide codes
- FASTA format specification