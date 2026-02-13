# Comprehensive tests for subpca function and S3 methods

library(testthat)
library(subpca)

test_that("subpca works with basic parameters", {
  set.seed(2001)
  X <- matrix(rnorm(50*100), 50, 100)
  clus <- rep(1:10, each = 10)
  
  result <- subpca(X, clus, ncomp = 5, ccomp = 2)
  
  expect_s3_class(result, "subpca")
  expect_s3_class(result, "metapca")
  expect_s3_class(result, "bi_projector")
  
  # Check dimensions
  v <- components(result)
  s <- scores(result)
  
  expect_equal(nrow(v), 100)  # number of variables
  expect_equal(ncol(v), 5)    # ncomp
  expect_equal(nrow(s), 50)   # number of observations
  expect_equal(ncol(s), 5)    # ncomp
})

test_that("subpca handles different combine methods", {
  set.seed(2002)
  X <- matrix(rnorm(40*80), 40, 80)
  clus <- rep(1:8, each = 10)
  
  methods <- c("pca", "scaled", "MFA")
  
  results <- list()
  for (method in methods) {
    result <- subpca(X, clus, 
                    ncomp = 4, 
                    ccomp = 2,
                    combine = method)
    
    expect_s3_class(result, "subpca")
    
    v <- components(result)
    s <- scores(result)
    
    expect_equal(ncol(v), 4)
    expect_equal(ncol(s), 4)
    expect_true(all(is.finite(v)))
    expect_true(all(is.finite(s)))
    
    results[[method]] <- result
  }
  
  # Different methods should give different results
  v_pca <- components(results[["pca"]])
  v_scaled <- components(results[["scaled"]])
  v_mfa <- components(results[["MFA"]])
  
  expect_false(isTRUE(all.equal(v_pca, v_scaled)))
  expect_false(isTRUE(all.equal(v_pca, v_mfa)))
})

test_that("subpca handles weights parameter", {
  set.seed(2003)
  X <- matrix(rnorm(30*60), 30, 60)
  clus <- rep(1:6, each = 10)
  
  # Equal weights
  weights_equal <- rep(1, 6)
  result_equal <- subpca(X, clus, 
                        weights = weights_equal,
                        ncomp = 3, 
                        ccomp = 2)
  
  # Different weights
  weights_diff <- c(2, 1, 1, 1, 1, 3)
  result_diff <- subpca(X, clus, 
                       weights = weights_diff,
                       ncomp = 3, 
                       ccomp = 2)
  
  expect_s3_class(result_equal, "subpca")
  expect_s3_class(result_diff, "subpca")
  
  # Results should be different with different weights
  v_equal <- components(result_equal)
  v_diff <- components(result_diff)
  
  expect_false(isTRUE(all.equal(v_equal, v_diff)))
})

test_that("subpca handles function-based ccomp", {
  set.seed(2004)
  X <- matrix(rnorm(35*70), 35, 70)
  clus <- rep(1:7, each = 10)
  
  # Function that determines components based on variance explained
  f <- function(p) {
    ev <- sdev(p)^2
    v <- cumsum(ev) / sum(ev)
    min(which(v > 0.8))  # 80% variance
  }
  
  result <- subpca(X, clus, ncomp = 4, ccomp = f)
  
  expect_s3_class(result, "subpca")
  
  v <- components(result)
  expect_equal(ncol(v), 4)
  expect_true(all(is.finite(v)))
})

test_that("subpca handles fractional ccomp", {
  set.seed(2005)
  X <- matrix(rnorm(40*80), 40, 80)
  clus <- rep(1:8, each = 10)
  
  # Use 50% of variance
  result <- subpca(X, clus, ncomp = 5, ccomp = 0.5)
  
  expect_s3_class(result, "subpca")
  
  v <- components(result)
  s <- scores(result)
  
  expect_equal(ncol(v), 5)
  expect_equal(ncol(s), 5)
})

test_that("subpca validates cluster sizes", {
  set.seed(2006)
  X <- matrix(rnorm(20*30), 20, 30)
  
  # Invalid: cluster with only 2 members
  clus_invalid <- c(1, 1, rep(2:10, each = 3), 10)
  
  expect_error(
    subpca(X, clus_invalid, ncomp = 2),
    "Each cluster must have at least 3 observations"
  )
  
  # Valid: all clusters have at least 3 members
  clus_valid <- rep(1:10, each = 3)
  result <- subpca(X, clus_valid, ncomp = 2, ccomp = 1)
  
  expect_s3_class(result, "subpca")
})

test_that("subpca print method works", {
  set.seed(2007)
  X <- matrix(rnorm(25*50), 25, 50)
  clus <- rep(1:5, each = 10)
  
  result <- subpca(X, clus, ncomp = 3, ccomp = 2)
  
  # Capture print output
  output <- capture.output(print(result))
  
  expect_true(any(grepl("subpca", output)))
  expect_true(any(grepl("number of clusters:", output)))
  expect_true(any(grepl("input dim:", output)))
  expect_true(any(grepl("output dim:", output)))
})

test_that("subpca project method works", {
  set.seed(2008)
  X <- matrix(rnorm(30*60), 30, 60)
  clus <- rep(1:6, each = 10)
  
  result <- subpca(X, clus, ncomp = 4, ccomp = 2)
  
  # New data for projection
  X_new <- matrix(rnorm(10*60), 10, 60)
  
  proj <- project(result, X_new)
  
  expect_equal(nrow(proj), 10)  # 10 new observations
  expect_equal(ncol(proj), 4)   # 4 components
  expect_true(all(is.finite(proj)))
})

test_that("subpca partial_project method works", {
  set.seed(2009)
  X <- matrix(rnorm(25*50), 25, 50)
  clus <- rep(1:5, each = 10)
  
  result <- subpca(X, clus, ncomp = 3, ccomp = 2)
  
  # New data for partial projection (only first 20 columns)
  X_new_partial <- matrix(rnorm(8*20), 8, 20)
  colind <- 1:20
  
  proj_partial <- partial_project(result, X_new_partial, colind = colind)
  
  expect_equal(nrow(proj_partial), 8)
  expect_equal(ncol(proj_partial), 3)
  expect_true(all(is.finite(proj_partial)))
})

test_that("subpca project_block method works", {
  set.seed(2010)
  X <- matrix(rnorm(30*60), 30, 60)
  clus <- rep(1:6, each = 10)
  
  result <- subpca(X, clus, ncomp = 4, ccomp = 2)
  
  # New data for single block (block 2 has 10 variables)
  X_new_block <- matrix(rnorm(5*10), 5, 10)
  
  proj_block <- project_block(result, X_new_block, block = 2)
  
  expect_equal(nrow(proj_block), 5)
  expect_equal(ncol(proj_block), 4)
})

test_that("subpca reconstruct method works", {
  set.seed(2011)
  X <- matrix(rnorm(20*40), 20, 40)
  clus <- rep(1:4, each = 10)
  
  result <- subpca(X, clus, ncomp = 3, ccomp = 2)
  
  # Full reconstruction
  X_recon <- reconstruct(result)
  
  expect_equal(dim(X_recon), dim(X))
  
  # Partial reconstruction with subset of components
  X_recon_partial <- reconstruct(result, comp = 1:2)
  
  expect_equal(dim(X_recon_partial), dim(X))
  
  # Reconstruction with subset of rows
  X_recon_rows <- reconstruct(result, rowind = 1:10)
  
  expect_equal(nrow(X_recon_rows), 10)
  expect_equal(ncol(X_recon_rows), 40)
  
  # Reconstruction with subset of columns
  X_recon_cols <- reconstruct(result, colind = 1:20)
  
  expect_equal(nrow(X_recon_cols), 20)
  expect_equal(ncol(X_recon_cols), 20)
})

test_that("subpca truncate method works", {
  set.seed(2012)
  X <- matrix(rnorm(25*50), 25, 50)
  clus <- rep(1:5, each = 10)
  
  result <- subpca(X, clus, ncomp = 5, ccomp = 2)
  
  # Truncate to fewer components
  result_trunc <- truncate(result, 3)
  
  expect_s3_class(result_trunc, "metapca")
  expect_s3_class(result_trunc, "metapca")
  
  v_trunc <- components(result_trunc)
  s_trunc <- scores(result_trunc)
  
  expect_equal(ncol(v_trunc), 3)
  expect_equal(ncol(s_trunc), 3)
})

test_that("subpca sdev method works", {
  set.seed(2013)
  X <- matrix(rnorm(30*60), 30, 60)
  clus <- rep(1:6, each = 10)
  
  result <- subpca(X, clus, ncomp = 4, ccomp = 2)
  
  sd <- sdev(result)
  
  expect_equal(length(sd), 4)
  expect_true(all(sd > 0))
  # Standard deviations should be in decreasing order
  expect_true(all(diff(sd) <= 0))
})

test_that("subpca handles preprocessing", {
  set.seed(2014)
  X <- matrix(rnorm(25*50, mean = 10, sd = 5), 25, 50)
  clus <- rep(1:5, each = 10)
  
  # Test with centering
  result_center <- subpca(X, clus, 
                         ncomp = 3, 
                         ccomp = 2,
                         preproc = center())
  
  expect_s3_class(result_center, "subpca")
  
  # Components should be valid
  v <- components(result_center)
  expect_true(all(is.finite(v)))
})

test_that("subpca handles edge case with minimum clusters", {
  set.seed(2015)
  # Minimum case: 2 clusters
  X <- matrix(rnorm(10*20), 10, 20)
  clus <- rep(1:2, each = 10)
  
  result <- subpca(X, clus, ncomp = 2, ccomp = 2)
  
  expect_s3_class(result, "subpca")
  expect_equal(attr(result, "nclus"), 2)
})

test_that("subpca handles large number of clusters", {
  set.seed(2016)
  X <- matrix(rnorm(50*150), 50, 150)
  clus <- rep(1:50, each = 3)  # 50 clusters, 3 variables each
  
  result <- subpca(X, clus, ncomp = 10, ccomp = 1)
  
  expect_s3_class(result, "subpca")
  expect_equal(attr(result, "nclus"), 50)
  
  v <- components(result)
  expect_equal(ncol(v), 10)
})

test_that("subpca residuals method works", {
  set.seed(2017)
  X <- matrix(rnorm(20*40), 20, 40)
  clus <- rep(1:4, each = 10)
  
  result <- subpca(X, clus, ncomp = 3, ccomp = 2)
  
  # Compute residuals
  X_recon <- reconstruct(result)
  resid <- X - X_recon
  
  expect_equal(dim(resid), dim(X))
  
  # Residuals should have reduced variance
  expect_true(sum(resid^2) < sum(X^2))
  
  # Mean of residuals should be near zero
  expect_true(abs(mean(resid)) < 0.1)
})
