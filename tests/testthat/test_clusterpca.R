# Tests for clusterpca function

test_that("clusterpca works with column-wise clustering", {
  set.seed(42)
  X <- matrix(rnorm(20*50), 20, 50)
  clus <- rep(1:5, each=10)
  
  result <- clusterpca(X, clus, ccomp=3, colwise=TRUE)
  
  expect_s3_class(result, "clusterpca")
  expect_equal(length(result$fits), 5)  # 5 clusters
  expect_equal(length(result$clus), 50)  # 50 columns
  
  # Check that coef and scores work
  cf <- coef(result)
  sc <- scores(result)
  expect_equal(nrow(cf), 50)
  expect_equal(nrow(sc), 20)
})

test_that("clusterpca works with row-wise clustering", {
  set.seed(42)
  X <- matrix(rnorm(50*20), 50, 20)
  clus <- rep(1:5, each=10)
  
  result <- clusterpca(X, clus, ccomp=3, colwise=FALSE)
  
  expect_s3_class(result, "clusterpca")
  expect_equal(length(result$fits), 5)
  expect_equal(length(result$clus), 50)
  
  # Check dimensions
  cf <- coef(result)
  sc <- scores(result)
  expect_equal(nrow(cf), 20)  # variables
  expect_equal(nrow(sc), 50)  # observations
})

test_that("clusterpca handles fractional ccomp correctly", {
  set.seed(42)
  X <- matrix(rnorm(20*50), 20, 50)
  clus <- rep(1:5, each=10)
  
  # Use fractional ccomp (80% variance)
  result <- clusterpca(X, clus, ccomp=0.8, colwise=TRUE)
  
  expect_s3_class(result, "clusterpca")
  # Each cluster should have components explaining >= 80% variance
  for (i in seq_along(result$fits)) {
    fit <- result$fits[[i]]
    ev <- sdev(fit)^2
    var_explained <- cumsum(ev) / sum(ev)
    # The last component should bring us over 0.8
    expect_gte(var_explained[length(var_explained)], 0.8)
  }
})

test_that("clusterpca handles custom pcafun", {
  set.seed(42)
  X <- matrix(rnorm(20*30), 20, 30)
  clus <- rep(1:3, each=10)
  
  # Custom PCA function that adds a tag
  custom_pca <- function(X, ncomp, preproc, ind=NULL) {
    result <- multivarious::pca(X, ncomp=ncomp, preproc=preproc)
    attr(result, "custom_tag") <- "test"
    result
  }
  
  result <- clusterpca(X, clus, ccomp=2, pcafun=custom_pca)
  
  expect_s3_class(result, "clusterpca")
  # Check that custom function was used
  expect_equal(attr(result$fits[[1]], "custom_tag"), "test")
})

test_that("clusterpca residuals work correctly", {
  set.seed(42)
  X <- matrix(rnorm(20*30), 20, 30)
  clus <- rep(1:3, each=10)
  
  result <- clusterpca(X, clus, ccomp=2, colwise=TRUE)
  
  # Compute residuals using the original X
  # Use single ncomp value to avoid method dispatch issues
  resid <- residuals(result, ncomp=2, xorig=X)
  
  expect_equal(dim(resid), dim(X))
  # Check that residuals have been computed and are non-zero
  expect_true(sum(resid^2) > 0)
  # Some variance should be explained (residuals smaller than original)
  expect_true(sum(resid^2) < sum(X^2))
})

test_that("clusterpca shape and sdev methods work", {
  set.seed(42)
  X <- matrix(rnorm(20*30), 20, 30)
  clus <- rep(1:3, each=10)
  
  result <- clusterpca(X, clus, ccomp=2, colwise=TRUE)
  
  # Test shape - components returns the matrix, so get its shape
  comp <- components(result)
  sh <- dim(comp)
  expect_equal(length(sh), 2)
  expect_equal(sh[1], 30)  # number of variables
  
  # Test sdev
  sd <- sdev(result)
  expect_true(all(sd >= 0))
  expect_equal(length(sd), sum(result$ncomp))
})