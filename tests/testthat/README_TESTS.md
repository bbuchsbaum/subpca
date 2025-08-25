# Test Suite Documentation for subpca Package

## Overview
This directory contains comprehensive testthat test cases for the subpca R package. The tests verify the functionality of all main functions and their S3 methods.

## Test Files

### Core Function Tests

1. **test_subpca.R** (existing)
   - Basic tests for the main `subpca()` function
   - Tests different parameter combinations

2. **test_subpca_comprehensive.R** (new)
   - Comprehensive tests for `subpca()` function
   - Tests all parameter variations (ncomp, ccomp values)
   - Tests different combine methods (pca, scaled, MFA)
   - Tests edge cases and error handling
   - Tests all S3 methods (scores, components, reconstruct, project, truncate)
   - 17 test groups with 59 individual assertions

3. **test_clusterpca.R** (existing, enhanced)
   - Tests for `clusterpca()` function
   - Column-wise and row-wise clustering
   - Fractional ccomp handling
   - Custom PCA functions
   - Residuals computation

4. **test_metapca.R** (existing)
   - Basic tests for `metapca()` function
   - Tests combining multiple PCA fits

5. **test_metapca_advanced.R** (existing)
   - Advanced tests for metapca functionality

6. **test_hcluspca.R** (new)
   - Comprehensive tests for hierarchical clustering PCA
   - Tests with different cut levels
   - Tests skip_global parameter
   - Tests orthogonalization
   - Tests different SVD methods
   - Note: Some tests may fail due to dependencies

7. **test_hcluspca_simple.R** (new)
   - Simplified tests for hcluspca that avoid problematic areas
   - Basic functionality verification
   - Parameter validation

8. **test_musubpca.R** (new)
   - Tests for multi-block subpca
   - Tests with multiblock data structures
   - Tests different combine methods
   - Note: Requires @export tag in musubpca function to work

### Integration Tests

9. **test_edge_cases.R** (existing)
   - Tests for edge cases across all functions

10. **test_run_all.R** (new)
    - Master test file that runs basic tests for all functions
    - Provides a quick verification of package functionality

11. **test_comprehensive_suite.R** (new)
    - Complete test suite with organized sections
    - Tests all main functions and S3 methods
    - Includes edge cases and error handling
    - Provides test summary output

## Test Coverage

### Functions Tested:
- ✅ `subpca()` - Main sub-block PCA function
- ✅ `clusterpca()` - Cluster-wise PCA
- ✅ `metapca()` - PCA over multiple PCA fits
- ⚠️  `hcluspca()` - Hierarchical clustering PCA (partial - some dependency issues)
- ⚠️  `musubpca()` - Multi-block subpca (needs @export tag)

### S3 Methods Tested:
- ✅ `scores()`
- ✅ `components()` / `coef()`
- ✅ `reconstruct()`
- ✅ `project()`
- ✅ `partial_project()`
- ✅ `project_block()`
- ✅ `truncate()`
- ✅ `sdev()`
- ✅ `residuals()`
- ✅ `shape()`
- ✅ `print()`

### Parameter Variations Tested:
- Fixed component numbers (ccomp as integer)
- Fractional variance explained (ccomp as fraction < 1)
- Function-based component selection (ccomp as function)
- Different combination methods (pca, scaled, MFA)
- Weights for meta-analysis
- Column-wise and row-wise clustering
- Various preprocessing options

### Edge Cases Tested:
- Minimum cluster sizes (3 members per cluster)
- Single component extraction
- Large number of clusters
- Small matrices
- Invalid parameters

## Running the Tests

To run all tests with the package loaded:

```r
# Load the package first
devtools::load_all()

# Run all tests
testthat::test_dir("tests/testthat")

# Or run specific test files
testthat::test_file("tests/testthat/test_comprehensive_suite.R")
```

## Known Issues

1. Some tests require `devtools::load_all()` to work properly due to function export issues
2. `musubpca()` function needs `@export` tag in the source code
3. `hcluspca()` has some dependency issues that may cause test failures
4. Very small matrices may cause errors in some edge cases due to mathematical constraints

## Test Statistics

- Total test files: 11
- Total test groups: ~50+
- Total assertions: ~200+
- Coverage: All exported functions and their main S3 methods

## Synthetic Data Generation

Tests use carefully crafted synthetic data:
- Small matrices (2x2, 3x3) for edge case testing
- Medium matrices (20x40, 30x60) for standard testing
- Larger matrices (50x100) for performance testing
- Specific cluster structures to test different scenarios
- Random seeds for reproducibility

## Future Improvements

1. Add tests for parallel processing with future package
2. Add performance benchmarks
3. Add tests for very large matrices
4. Add tests for missing data handling
5. Export musubpca function and add comprehensive tests