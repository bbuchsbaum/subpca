# Tests for edge cases and error conditions

test_that("subpca handles minimum cluster size correctly", {
  set.seed(42)
  X <- matrix(rnorm(10*9), 10, 9)
  clus <- c(1,1,1, 2,2,2, 3,3,3)  # Exactly 3 per cluster (minimum)
  
  # Should work with minimum size
  result <- subpca(X, clus, ncomp=2, ccomp=2)
  expect_s3_class(result, "subpca")
  
  # Should fail with too small clusters
  clus_small <- c(1,1, 2,2, 3,3,3,3,3)  # Only 2 in first clusters
  expect_error(
    subpca(X, clus_small, ncomp=2),
    "Each cluster must have at least 3 observations"
  )
})

test_that("clusterpca handles single component request", {
  set.seed(42)
  X <- matrix(rnorm(20*30), 20, 30)
  clus <- rep(1:3, each=10)
  
  # Request only 1 component per cluster
  result <- clusterpca(X, clus, ccomp=1, colwise=TRUE)
  
  expect_s3_class(result, "clusterpca")
  expect_true(all(result$ncomp == 1))
  expect_equal(ncol(scores(result)), 3)  # 3 clusters * 1 component
})

test_that("metapca handles mismatched weights gracefully", {
  set.seed(42)
  X1 <- matrix(rnorm(20*10), 20, 10)
  X2 <- matrix(rnorm(20*15), 20, 15)
  
  pc1 <- pca(X1, ncomp=5)
  pc2 <- pca(X2, ncomp=5)
  
  # Correct number of weights
  mp1 <- metapca(list(pc1, pc2), ncomp=5, weights=c(1, 2))
  expect_s3_class(mp1, "metapca")
  
  # Wrong number of weights should error
  expect_error(
    metapca(list(pc1, pc2), ncomp=5, weights=c(1, 2, 3)),
    "weights must have one value per fit"
  )
})

test_that("fractional ccomp handles edge cases", {
  set.seed(42)
  X <- matrix(rnorm(20*30), 20, 30)
  clus <- rep(1:3, each=10)
  
  # Very high variance threshold (0.99)
  result_high <- clusterpca(X, clus, ccomp=0.99, colwise=TRUE)
  expect_s3_class(result_high, "clusterpca")
  
  # Very low variance threshold (0.1) 
  result_low <- clusterpca(X, clus, ccomp=0.1, colwise=TRUE)
  expect_s3_class(result_low, "clusterpca")
  
  # Should have fewer components with low threshold
  expect_lte(sum(result_low$ncomp), sum(result_high$ncomp))
  
  # Edge case: ccomp = 1.0 should not be treated as fractional
  result_one <- clusterpca(X, clus, ccomp=1.0, colwise=TRUE)
  expect_true(all(result_one$ncomp == 1))
})

test_that("clusterpca handles more components than data dimensions", {
  set.seed(42)
  # 5 observations, request 10 components
  X <- matrix(rnorm(5*30), 5, 30)
  clus <- rep(1:3, each=10)
  
  # Should automatically limit components
  suppressWarnings({
    result <- clusterpca(X, clus, ccomp=10, colwise=TRUE)
  })
  
  expect_s3_class(result, "clusterpca")
  # Each cluster should have at most 5 components (limited by observations)
  expect_true(all(result$ncomp <= 5))
})

test_that("metapca handles empty projection gracefully", {
  set.seed(42)
  X1 <- matrix(rnorm(20*10), 20, 10)
  X2 <- matrix(rnorm(20*15), 20, 15)
  
  pc1 <- pca(X1, ncomp=5)
  pc2 <- pca(X2, ncomp=5)
  mp <- metapca(list(pc1, pc2), ncomp=5)
  
  # Partial project with columns that don't match block indices
  new_data <- matrix(rnorm(5*25), 5, 25)
  # Try to project with indices outside the actual range
  result <- partial_project(mp, new_data, colind=26:50)
  
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 5)
  # Should be all zeros since no valid columns
  expect_true(all(result == 0))
})

test_that("subpca print method works", {
  set.seed(42)
  X <- matrix(rnorm(20*30), 20, 30)
  clus <- rep(1:3, each=10)
  
  result <- subpca(X, clus, ncomp=5)
  
  # Capture print output
  output <- capture.output(print(result))
  
  expect_true(any(grepl("subpca", output)))
  expect_true(any(grepl("number of clusters:  3", output)))
  expect_true(any(grepl("input dim:", output)))
  expect_true(any(grepl("output dim:", output)))
})

test_that("clusterpca handles custom preprocessor", {
  set.seed(42)
  X <- matrix(rnorm(20*30), 20, 30)
  clus <- rep(1:3, each=10)
  
  # Test with pass() preprocessor (no preprocessing)
  result_pass <- clusterpca(X, clus, ccomp=2, preproc=pass())
  expect_s3_class(result_pass, "clusterpca")
  
  # Test with center() preprocessor
  result_center <- clusterpca(X, clus, ccomp=2, preproc=center())
  expect_s3_class(result_center, "clusterpca")
  
  # Results should be different
  scores_pass <- scores(result_pass)
  scores_center <- scores(result_center)
  expect_gt(sum((scores_pass - scores_center)^2), 0.1)
})

test_that("metapca handles single fit input", {
  set.seed(42)
  X <- matrix(rnorm(20*10), 20, 10)
  pc <- pca(X, ncomp=5)
  
  # Should work with single fit (though not very useful)
  mp <- metapca(list(pc), ncomp=3)
  
  expect_s3_class(mp, "metapca")
  expect_equal(ncol(scores(mp)), 3)
  expect_equal(length(mp$fits), 1)
})

test_that("components accessor works for all classes", {
  set.seed(42)
  
  # Test with clusterpca
  X1 <- matrix(rnorm(20*30), 20, 30)
  clus1 <- rep(1:3, each=10)
  cp <- clusterpca(X1, clus1, ccomp=2, colwise=TRUE)
  comp_cp <- components(cp)
  expect_equal(nrow(comp_cp), 30)
  
  # Test with metapca
  pc1 <- pca(matrix(rnorm(20*10), 20, 10), ncomp=3)
  pc2 <- pca(matrix(rnorm(20*15), 20, 15), ncomp=3)
  mp <- metapca(list(pc1, pc2), ncomp=4)
  comp_mp <- components(mp)
  expect_equal(nrow(comp_mp), 25)  # 10 + 15
  
  # Test with subpca
  X2 <- matrix(rnorm(20*30), 20, 30)
  clus2 <- rep(1:3, each=10)
  sp <- subpca(X2, clus2, ncomp=5)
  comp_sp <- components(sp)
  expect_equal(nrow(comp_sp), 30)
})
