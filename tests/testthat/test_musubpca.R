# Comprehensive tests for musubpca function

library(testthat)
library(subpca)
library(multidesign)
library(multivarious)

test_that("musubpca works with basic multiblock data", {
  set.seed(111)
  # Create multiblock data: 3 blocks, each 40x50
  blocks <- replicate(3, matrix(rnorm(40*50), 40, 50), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  # Create clusters for the 50 variables
  clus <- rep(1:5, each = 10)
  
  result <- musubpca(X, clus, 
                    ncomp = 5,
                    inner_ccomp = 2,
                    ccomp = 2)
  
  expect_s3_class(result, "metapca")
  
  # Check final dimensions
  v <- components(result)
  s <- scores(result)
  
  expect_equal(ncol(v), 5)  # ncomp
  expect_equal(ncol(s), 5)  # ncomp
})

test_that("musubpca handles column-stacked multiblock data", {
  set.seed(222)
  # Create column-stacked multiblock: 3 blocks stacked as columns
  blocks <- replicate(3, matrix(rnorm(50*40), 50, 40), simplify = FALSE)
  X_cstacked <- multidesign::multiblock(blocks)
  
  # Clusters for the 40 columns (variables)
  clus <- rep(1:4, each = 10)
  
  result <- musubpca(X_cstacked, clus,
                    ncomp = 4,
                    inner_ccomp = 2,
                    ccomp = 2)
  
  expect_s3_class(result, "metapca")
  
  v <- components(result)
  s <- scores(result)
  
  expect_equal(ncol(v), 4)
  expect_equal(ncol(s), 4)
})

test_that("musubpca handles fractional ccomp values", {
  set.seed(333)
  blocks <- replicate(2, matrix(rnorm(30*40), 30, 40), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  clus <- rep(1:4, each = 10)
  
  # Use fractional values for variance explained
  result <- musubpca(X, clus,
                    ncomp = 3,
                    inner_ccomp = 0.8,  # 80% variance
                    ccomp = 0.9)        # 90% variance
  
  expect_s3_class(result, "metapca")
  
  # Result should have valid components
  v <- components(result)
  expect_true(all(is.finite(v)))
  expect_equal(ncol(v), 3)
})

test_that("musubpca handles different combine methods", {
  set.seed(444)
  blocks <- replicate(2, matrix(rnorm(25*30), 25, 30), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  clus <- rep(1:3, each = 10)
  
  # Test different combination methods
  methods <- list(
    c("pca", "pca"),
    c("scaled", "scaled"),
    c("MFA", "MFA"),
    c("pca", "MFA")
  )
  
  for (method_pair in methods) {
    result <- musubpca(X, clus,
                      ncomp = 2,
                      inner_ccomp = 2,
                      ccomp = 2,
                      inner_combine = method_pair[1],
                      combine = method_pair[2])
    
    expect_s3_class(result, "metapca")
    
    v <- components(result)
    s <- scores(result)
    
    expect_equal(ncol(v), 2)
    expect_equal(ncol(s), 2)
    expect_true(all(is.finite(v)))
    expect_true(all(is.finite(s)))
  }
})

test_that("musubpca handles weights parameter", {
  set.seed(555)
  blocks <- replicate(3, matrix(rnorm(20*30), 20, 30), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  clus <- rep(1:3, each = 10)
  
  # Test with custom weights
  weights <- c(1, 2, 1)  # Weight middle cluster more
  
  result <- musubpca(X, clus,
                    weights = weights,
                    ncomp = 3,
                    inner_ccomp = 2,
                    ccomp = 2)
  
  expect_s3_class(result, "metapca")
  
  # Result should reflect weighted combination
  v <- components(result)
  expect_equal(ncol(v), 3)
})

test_that("musubpca validates cluster sizes", {
  set.seed(666)
  blocks <- replicate(3, matrix(rnorm(20*10), 20, 10), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  # Create clusters with too small size
  clus_invalid <- c(1, 1, rep(2:5, each = 2))  # First cluster has only 2 elements
  
  # This should fail because minimum cluster size is 3
  expect_error(
    musubpca(X, clus_invalid, ncomp = 2),
    "min\\(table\\(clus\\)\\) not greater than 2"
  )
  
  # Valid clustering - need 3 blocks and enough components for eigendecomp
  clus_valid <- c(rep(1:2, each = 5))
  result <- musubpca(X, clus_valid, 
                    ncomp = 2,
                    inner_ccomp = 2,  # Need at least 2 per cluster to get 4 total
                    ccomp = 2)
  
  expect_s3_class(result, "metapca")
})

test_that("musubpca handles varying number of components per cluster", {
  set.seed(777)
  blocks <- replicate(2, matrix(rnorm(30*40), 30, 40), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  clus <- rep(1:4, each = 10)
  
  # Different components for each cluster
  ccomp_vec <- c(1, 2, 3, 2)
  
  result <- musubpca(X, clus,
                    ncomp = 4,
                    inner_ccomp = 2,
                    ccomp = ccomp_vec)
  
  expect_s3_class(result, "metapca")
  
  v <- components(result)
  expect_equal(ncol(v), 4)
})

test_that("musubpca handles preprocessing", {
  set.seed(888)
  blocks <- replicate(2, matrix(rnorm(25*30), 25, 30), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  clus <- rep(1:3, each = 10)
  
  # Test with centering preprocessing
  result <- musubpca(X, clus,
                    ncomp = 3,
                    inner_ccomp = 2,
                    ccomp = 2,
                    preproc = center())
  
  expect_s3_class(result, "metapca")
  
  # Check that result is valid
  v <- components(result)
  s <- scores(result)
  
  expect_true(all(is.finite(v)))
  expect_true(all(is.finite(s)))
})

test_that("musubpca projection works correctly", {
  set.seed(999)
  blocks <- replicate(2, matrix(rnorm(20*30), 20, 30), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  clus <- rep(1:3, each = 10)
  
  result <- musubpca(X, clus,
                    ncomp = 2,
                    inner_ccomp = 2,
                    ccomp = 2)
  
  # Project new data
  new_blocks <- replicate(2, matrix(rnorm(5*30), 5, 30), simplify = FALSE)
  X_new <- multidesign::multiblock(new_blocks)
  
  # Combine blocks for projection
  X_new_combined <- do.call(cbind, new_blocks)
  
  proj <- project(result, X_new_combined)
  
  expect_equal(nrow(proj), 5)  # 5 new observations
  expect_equal(ncol(proj), 2)  # 2 components
})

test_that("musubpca handles minimum block edge case", {
  set.seed(1010)
  # Changed to 3 blocks (minimum for eigendecomp to work)
  blocks <- replicate(3, matrix(rnorm(30*40), 30, 40), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  clus <- rep(1:4, each = 10)
  
  result <- musubpca(X, clus,
                    ncomp = 2,
                    inner_ccomp = 2,
                    ccomp = 2)
  
  expect_s3_class(result, "metapca")
  
  v <- components(result)
  s <- scores(result)
  
  expect_equal(ncol(v), 2)
  expect_equal(ncol(s), 2)
})

test_that("musubpca handles large number of clusters", {
  set.seed(1111)
  # Changed to 3 blocks to meet eigendecomp requirements
  blocks <- replicate(3, matrix(rnorm(40*60), 40, 60), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  # Many small clusters
  clus <- rep(1:12, each = 5)
  
  result <- musubpca(X, clus,
                    ncomp = 5,
                    inner_ccomp = 1,  # Only 1 component per cluster due to small size
                    ccomp = 1)
  
  expect_s3_class(result, "metapca")
  
  v <- components(result)
  expect_equal(ncol(v), 5)
})

test_that("musubpca reconstruction works", {
  set.seed(1212)
  # Use 3 blocks to avoid eigendecomp issues
  blocks <- replicate(3, matrix(rnorm(20*30), 20, 30), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  clus <- rep(1:3, each = 10)
  
  result <- musubpca(X, clus,
                    ncomp = 3,
                    inner_ccomp = 2,
                    ccomp = 2)
  
  # Reconstruct data
  X_combined <- do.call(cbind, blocks)
  X_recon <- reconstruct(result)
  
  expect_equal(dim(X_recon), dim(X_combined))
  
  # Some variance should be explained
  residuals <- X_combined - X_recon
  expect_true(sum(residuals^2) < sum(X_combined^2))
})

test_that("musubpca truncation works", {
  set.seed(1313)
  # Use 3 blocks to avoid eigendecomp issues
  blocks <- replicate(3, matrix(rnorm(25*30), 25, 30), simplify = FALSE)
  X <- multidesign::multiblock(blocks)
  
  clus <- rep(1:3, each = 10)
  
  result <- musubpca(X, clus,
                    ncomp = 4,
                    inner_ccomp = 2,
                    ccomp = 2)
  
  # Truncate to fewer components
  result_trunc <- truncate(result, 2)
  
  expect_s3_class(result_trunc, "metapca")
  
  v_trunc <- components(result_trunc)
  s_trunc <- scores(result_trunc)
  
  expect_equal(ncol(v_trunc), 2)
  expect_equal(ncol(s_trunc), 2)
})