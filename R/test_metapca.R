library(multivarious)

test_that("can run a metapca analysis", {
  X1 <- matrix(rnorm(20*10), 20, 10)
  X2 <- matrix(rnorm(20*20), 20, 20)
  X3 <- matrix(rnorm(20*30), 20, 30)

  pc1 <- pca(X1, ncomp=10, preproc=center())
  pc2 <- pca(X2, ncomp=10, preproc=center())
  pc3 <- pca(X3, ncomp=10, preproc=center())

  fits <- list(pc1,pc2,pc3)
  pfit <- metapca(fits, ncomp=15)
  testthat::expect_true(!is.null(pfit))

})
