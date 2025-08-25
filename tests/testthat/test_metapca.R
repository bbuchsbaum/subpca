# Tests are run with the package already loaded by testthat
# No need to explicitly load libraries here

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

test_that("metapca is and pca match up", {
  X1 <- matrix(rnorm(20*10), 20, 10)
  X2 <- matrix(rnorm(20*20), 20, 20)
  X3 <- matrix(rnorm(20*30), 20, 30)

  pc1 <- pca(X1, ncomp=10, preproc=center())
  pc2 <- pca(X2, ncomp=20, preproc=center())
  pc3 <- pca(X3, ncomp=30, preproc=center())

  fits <- list(pc1,pc2,pc3)
  pfit <- metapca(fits, ncomp=15)
  pfit2 <- pca(cbind(X1,X2,X3), preproc=center(), ncomp=15)

  sc1 <- scores(pfit)
  sc2 <- scores(pfit2)

  p1 <- project(pfit, new_data=cbind(X1,X2,X3))
  p2 <- project(pfit2, new_data=cbind(X1,X2,X3))

  testthat::expect_true(sum(abs(sc1) - abs(sc2)) < 1e-5)
  testthat::expect_true(sum(abs(p1) - abs(p2)) < 1e-5)

})

test_that("metapca with different combining methods", {
  X1 <- matrix(rnorm(20*10), 20, 10)
  X2 <- matrix(rnorm(20*20), 20, 20)
  X3 <- matrix(rnorm(20*30), 20, 30)

  pc1 <- pca(X1, ncomp=10, preproc=center())
  pc2 <- pca(X2, ncomp=10, preproc=center())
  pc3 <- pca(X3, ncomp=10, preproc=center())

  fits <- list(pc1,pc2,pc3)
  pfit1 <- metapca(fits, ncomp=15, combine="pca")
  pfit2 <- metapca(fits, ncomp=15, combine="MFA")
  pfit3 <- metapca(fits, ncomp=15, combine="scaled")

  testthat::expect_true(!is.null(pfit1))
  testthat::expect_true(!is.null(pfit2))
  testthat::expect_true(!is.null(pfit3))

})


test_that("metapca with weights", {
  X1 <- matrix(rnorm(20*10), 20, 10)
  X2 <- matrix(rnorm(20*20), 20, 20)
  X3 <- matrix(rnorm(20*30), 20, 30)

  pc1 <- pca(X1, ncomp=10, preproc=center())
  pc2 <- pca(X2, ncomp=10, preproc=center())
  pc3 <- pca(X3, ncomp=10, preproc=center())

  fits <- list(pc1,pc2,pc3)
  pfit1 <- metapca(fits, ncomp=15,
                   weights=c(1,2,3), combine="pca")

  testthat::expect_true(!is.null(pfit1))

})


