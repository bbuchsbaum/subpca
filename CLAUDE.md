# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R package called `subpca` that provides methods for computing PCA within sub-blocks of data and then a second-level ("metapca") analysis over each cluster-wise PCA. The package extends the `multivarious` package functionality.

## Development Commands

### Testing
- Run all tests: `devtools::test()`
- Run specific test file: `testthat::test_file("tests/testthat/test_metapca.R")`
- Load package for development: `devtools::load_all()`

### Package Building
- Check package: `devtools::check()`
- Build documentation: `devtools::document()`
- Install package locally: `devtools::install()`

## Architecture

### Core Components

The package implements hierarchical PCA approaches through several main functions:

1. **clusterpca** (`R/clusterpca.R`) - Fits separate PCA models for each cluster/block of data, supporting both row-wise and column-wise clustering
2. **metapca** (`R/metapca.R`) - Performs PCA over multiple fitted PCA models, enabling hierarchical dimensionality reduction
3. **hcluspca** (`R/hcluspca.R`) - Hierarchical clustering-based PCA implementation
4. **subpca** (`R/subpca.R`) - Base implementation for sub-block PCA analysis
5. **musubpca** (`R/musubpca.R`) - Multi-unit sub-block PCA variant

### Key Dependencies

- **multivarious**: Core PCA functionality and bi_projector classes (imported from GitHub: bbuchsbaum/multivarious)
- **genpca**: Generalized PCA implementation
- **assertthat**: Input validation
- **furrr**: Parallel processing support via future_map
- **dclust**: Clustering functionality
- **Matrix**: Sparse matrix support

### S3 Methods

The package implements S3 methods for its custom classes (clusterpca, metapca) including:
- `coef()`, `ncomp()`, `scores()`, `sdev()`, `shape()`
- `project()`, `project_block()`, `partial_project()`
- `reconstruct()`, `residuals()`, `truncate()`

### Generic Functions

All generic functions from multivarious are re-exported in `R/all_generic.R` to maintain API consistency.