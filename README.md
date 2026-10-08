# Quantum Dot Spectroscopy Analysis

**Julia-based scientific data processing, spectral analysis, and visualization for semiconductor quantum dots.**

## Overview

This project contains Julia functions for processing and visualizing numerical results from computational spectroscopy studies of semiconductor quantum dots. It reads structured simulation-output tables, reconstructs two-dimensional X-ray/infrared spectral probability distributions, applies numerical normalization and binning, extracts spectral features, and generates figures and summary tables.

**Physical chemistry application:** The analysis investigates the distribution of infrared emission frequencies associated with X-ray excitation in CdTe quantum dots. It supports comparisons between coherent and incoherent probability formulations, and between zero-order and electron–hole-corrected results, to examine how the calculated spectral response changes across computational treatments and dot sizes.

> **Scope:** This repository is a post-processing and visualization tool. It does not execute DFT or Hartree–Fock electronic-structure calculations; it analyzes numerical results produced by upstream calculations.

## Key Features

- **Structured data ingestion:** Reads whitespace-delimited numerical output with a header row and columns for X-ray energy, infrared energy, and probability values, then reconstructs a two-dimensional matrix.
- **Numerical processing:** Performs trapezoidal integration, probability normalization, and two-dimensional binning.
- **Spectral analysis:** Computes marginal and conditional distributions, identifies high-probability spectral features, and extracts peak locations.
- **Comparison of computational results:** Visualizes coherent and incoherent distributions, their differences, and zero-order versus electron–hole-corrected results.
- **Scientific visualization:** Produces 2D heatmaps, line plots, spectral slices, and 3D surfaces with CairoMakie.
- **Results export:** Saves figures (including PDF output) and writes numerical peak summaries and table-oriented text files.

## Tech Stack and Scientific Domain

| Category | Details |
| --- | --- |
| Language | Julia |
| Visualization | CairoMakie, LaTeXStrings |
| Other imported packages | Random, LinearAlgebra, SpecialFunctions, Dates, Distributions, Printf, Statistics |
| Scientific context | Computational spectroscopy; DFT and Hartree–Fock-derived results; semiconductor nanomaterials |
| System studied in the included workflows | CdTe quantum dots |
| Analysis methods | 2D probability distributions, numerical integration, binning, peak extraction, coherent/incoherent comparison |

## Input and Output

**Expected input:** A text file with one header line followed by whitespace-separated numeric data. The reader expects three columns in the order shown below, and it assumes the rows form a complete, consistently ordered two-dimensional grid.

```text
omega_x  omega_ir  probability
...      ...       ...
```

Energy axes in the supplied workflows are converted from atomic units to electronvolts for plotting using a conversion factor of `27.211`.

**Outputs:** Depending on the selected analysis function, the code saves probability heatmaps, spectral line plots, comparison figures, 3D surfaces, peak summaries, and text tables. Figure filenames and output directories are defined in the analysis routines.

## How to Use

1. Install Julia and the external packages imported by the script (`CairoMakie`, `LaTeXStrings`, `SpecialFunctions`, and `Distributions`). The remaining imports are Julia standard libraries.
2. Place compatible simulation-output files in the directories expected by the selected analysis routine, or edit the input paths to point to your data.
3. Save the Julia source as `spectroscopy_analysis.jl` (or use its existing filename).
4. Run the script from the repository directory:

   ```bash
   julia spectroscopy_analysis.jl
   ```

**Current execution behavior:** The supplied source ends by calling `test4A()`. This workflow reads both zero-order and electron–hole-corrected CdTe probability results and generates a comparison plot. Other analysis routines can be run by changing that final function call. The `test*` names are analysis workflows, not an automated unit-test suite.

**Reproducibility note:** The current script contains hard-coded relative input/output paths and selected index ranges. The original simulation outputs, a Julia environment lockfile, and runnable example data are not included with this README. Users must provide appropriately structured inputs and verify paths, grid sizes, and package versions before running it on a different dataset.

## Project Status

Research-oriented analysis script. The functions support multiple post-processing workflows, but the code has not yet been packaged as a general-purpose Julia library or validated as a turnkey application for arbitrary simulation data.
