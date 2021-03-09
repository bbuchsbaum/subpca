library(multivarious)
library(assertthat)
library(subpca)

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

test_that("subpca is like regular pca", {
  X1 <- matrix(rnorm(20*10), 20, 10)
  clus <- rep(1:3, length.out=ncol(X1))
  weights <- c(1,2,3)
  ccomp <- 4
  ret <- subpca(X1, clus, ncomp=10, weights=weights, ccomp=ccomp, preproc=center())
  ret2 <- pca(X1, ncomp=10, preproc=center())
  expect_true(sum(abs(scores(ret)) - abs(scores(ret2))) < 1e-5)
})

