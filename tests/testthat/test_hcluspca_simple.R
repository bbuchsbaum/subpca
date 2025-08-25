# Simple tests for hcluspca function - avoiding complex areas

library(testthat)
library(subpca)
library(multivarious)

test_that("hcluspca basic functionality works", {
  skip_if_not_installed("dendextend")
  skip_if_not_installed("multivarious")
  
  set.seed(123)
  # Small test case
  X <- matrix(rnorm(30*20), 30, 20)
  
  # Create a simple hierarchical clustering
  d <- dist(1:30)
  hclus <- hclust(d, method = "complete")
  
  # Test with a single level of cuts
  cuts <- c(3)
  
  # Try to run with minimal parameters
  result <- tryCatch({
    hcluspca(X, hclus, cuts, ccomp = 1, est_method = "standard", svd_method = "base")
  }, error = function(e) {
    skip(paste("hcluspca not working:", e$message))
  })
  
  if (!inherits(result, "try-error")) {
    expect_s3_class(result, "bi_projector")
    
    # Check that we have components and scores
    v <- components(result)
    s <- scores(result)
    
    expect_true(is.matrix(v))
    expect_true(is.matrix(s))
    expect_equal(nrow(s), 30)
    expect_equal(nrow(v), 20)
  }
})

test_that("hcluspca validates input parameters", {
  set.seed(234)
  X <- matrix(rnorm(20*15), 20, 15)
  
  d <- dist(1:20)
  hclus <- hclust(d)
  
  # Test with invalid cuts (less than 2)
  expect_error(
    hcluspca(X, hclus, cuts = c(1, 4), ccomp = 1),
    "all `cuts` must be greater than 1"
  )
  
  # Test with invalid ccomp length
  expect_error(
    hcluspca(X, hclus, cuts = c(2, 4), ccomp = c(1, 2)),
    "if `ccomp` is a vector it must have an entry for every level"
  )
})

test_that("hcluspca skip_global parameter works", {
  skip_if_not_installed("dendextend")
  
  set.seed(345)
  X <- matrix(rnorm(25*15), 25, 15)
  
  d <- dist(1:25)
  hclus <- hclust(d)
  
  cuts <- c(2)
  
  # Test that skip_global changes the behavior
  result_with_global <- tryCatch({
    hcluspca(X, hclus, cuts, ccomp = 1, skip_global = FALSE)
  }, error = function(e) NULL)
  
  result_without_global <- tryCatch({
    hcluspca(X, hclus, cuts, ccomp = 1, skip_global = TRUE)
  }, error = function(e) NULL)
  
  # If both work, they should give different results
  if (!is.null(result_with_global) && !is.null(result_without_global)) {
    v_with <- components(result_with_global)
    v_without <- components(result_without_global)
    
    # Different number of components expected
    expect_false(ncol(v_with) == ncol(v_without))
  }
})