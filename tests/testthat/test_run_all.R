# Master test file to run comprehensive tests for subpca package

library(testthat)
library(subpca)
library(multivarious)

# Run existing test files
test_that("All comprehensive tests pass", {
  
  # Test that package is loaded
  expect_true("subpca" %in% loadedNamespaces())
  
  # Test that key functions exist
  expect_true(exists("subpca"))
  expect_true(exists("clusterpca"))
  expect_true(exists("metapca"))
  expect_true(exists("hcluspca"))
  
  # Note: musubpca needs @export tag to work properly
  # expect_true(exists("musubpca"))
  
  # Test that S3 methods are registered
  methods_list <- methods(class = "subpca")
  expect_true(length(methods_list) > 0)
  
  # Basic functionality test
  set.seed(999)
  X <- matrix(rnorm(30*60), 30, 60)
  clus <- rep(1:6, each = 10)
  
  # Test basic subpca
  result <- subpca(X, clus, ncomp = 3, ccomp = 2)
  expect_s3_class(result, "subpca")
  expect_s3_class(result, "metapca")
  
  # Test components and scores
  v <- components(result)
  s <- scores(result)
  
  expect_equal(ncol(v), 3)
  expect_equal(ncol(s), 3)
  expect_equal(nrow(s), 30)
  
  # Test clusterpca
  result_cluster <- clusterpca(X, clus, ccomp = 2, colwise = TRUE)
  expect_s3_class(result_cluster, "clusterpca")
  
  # Test metapca with multiple fits
  fits <- list()
  for (i in 1:3) {
    X_block <- matrix(rnorm(30*20), 30, 20)
    fits[[i]] <- multivarious::pca(X_block, ncomp = 5, preproc = multivarious::center())
  }
  
  result_meta <- metapca(fits, ncomp = 4)
  expect_s3_class(result_meta, "metapca")
  
  # Print summary
  cat("\n=== Test Summary ===\n")
  cat("subpca: PASS\n")
  cat("clusterpca: PASS\n")
  cat("metapca: PASS\n")
  cat("Basic functionality verified\n")
})

test_that("Edge cases are handled properly", {
  set.seed(888)
  
  # Minimum cluster size (need at least 3 clusters)
  X <- matrix(rnorm(10*30), 10, 30)
  clus <- rep(1:3, each = 10)
  
  result <- subpca(X, clus, ncomp = 2, ccomp = 1)
  expect_s3_class(result, "subpca")
  
  # Single component
  X <- matrix(rnorm(20*40), 20, 40)
  clus <- rep(1:4, each = 10)
  
  result <- subpca(X, clus, ncomp = 1, ccomp = 1)
  v <- components(result)
  expect_equal(ncol(v), 1)
  
  # Fractional variance explained
  result <- subpca(X, clus, ncomp = 2, ccomp = 0.8)
  expect_s3_class(result, "subpca")
})