# Tests are run with the package already loaded by testthat
# No need to explicitly load libraries here

test_that("can run a subpca analysis", {
  X1 <- matrix(rnorm(20*10), 20, 10)
  clus <- rep(1:3, length.out=ncol(X1))
  ret <- subpca(X1, clus, ncomp=6)
  expect_true(!is.null(ret))
})

test_that("can run a subpca analysis with weights", {
  X1 <- matrix(rnorm(20*10), 20, 10)
  clus <- rep(1:3, length.out=ncol(X1))
  weights <- c(1,2,3)
  ret <- subpca(X1, clus, ncomp=6, weights=weights)
  expect_true(!is.null(ret))

})

test_that("can run a subpca analysis with ccomp", {
  X1 <- matrix(rnorm(20*10), 20, 10)
  clus <- rep(1:3, length.out=ncol(X1))
  weights <- c(1,2,3)
  ccomp <- c(1,1,1)
  ret <- subpca(X1, clus, ncomp=6, weights=weights, ccomp=ccomp)
  expect_true(!is.null(ret))
})

test_that("subpca produces valid output", {
  X1 <- matrix(rnorm(20*10), 20, 10)
  clus <- rep(1:3, length.out=ncol(X1))
  weights <- c(1,1,1)
  ccomp <- 3
  ret <- subpca(X1, clus, ncomp=5, weights=weights, ccomp=ccomp, preproc=center())
  
  # Check that scores have correct dimensions
  expect_equal(nrow(scores(ret)), nrow(X1))
  expect_lte(ncol(scores(ret)), 5)  # At most 5 components as requested
  
  # Check that components have correct dimensions  
  expect_equal(nrow(components(ret)), ncol(X1))
  expect_equal(ncol(components(ret)), ncol(scores(ret)))
})

