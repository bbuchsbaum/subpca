# Advanced tests for metapca functions

test_that("partial_project.metapca works correctly", {
  set.seed(123)
  # Create test data with different sized blocks
  X1 <- matrix(rnorm(20*10), 20, 10)
  X2 <- matrix(rnorm(20*15), 20, 15)
  X3 <- matrix(rnorm(20*20), 20, 20)
  
  # Create PCA fits
  pc1 <- pca(X1, ncomp=5, preproc=center())
  pc2 <- pca(X2, ncomp=5, preproc=center())
  pc3 <- pca(X3, ncomp=5, preproc=center())
  
  # Create metapca
  fits <- list(pc1, pc2, pc3)
  mp <- metapca(fits, ncomp=10)
  
  # Test partial projection on subset of columns
  new_data <- matrix(rnorm(5*45), 5, 45)  # 45 = 10+15+20
  
  # Project using only first and third blocks (columns 1:10 and 26:45)
  colind <- c(1:10, 26:45)
  partial_proj <- partial_project(mp, new_data[, colind], colind=colind)
  
  expect_equal(nrow(partial_proj), 5)
  expect_equal(ncol(partial_proj), ncol(scores(mp)))
  
  # Project using only middle block
  colind2 <- 11:25
  partial_proj2 <- partial_project(mp, new_data[, colind2], colind=colind2)
  
  expect_equal(nrow(partial_proj2), 5)
  expect_equal(ncol(partial_proj2), ncol(scores(mp)))
})

test_that("project_block.metapca works correctly", {
  set.seed(123)
  X1 <- matrix(rnorm(20*10), 20, 10)
  X2 <- matrix(rnorm(20*15), 20, 15)
  
  pc1 <- pca(X1, ncomp=5, preproc=center())
  pc2 <- pca(X2, ncomp=5, preproc=center())
  
  mp <- metapca(list(pc1, pc2), ncomp=8)
  
  # Test projecting block 1
  new_block1 <- matrix(rnorm(3*10), 3, 10)
  proj1 <- project_block(mp, new_block1, block=1)
  
  expect_equal(nrow(proj1), 3)
  expect_equal(ncol(proj1), ncol(scores(mp)))
  
  # Test projecting block 2
  new_block2 <- matrix(rnorm(3*15), 3, 15)
  proj2 <- project_block(mp, new_block2, block=2)
  
  expect_equal(nrow(proj2), 3)
  expect_equal(ncol(proj2), ncol(scores(mp)))
})

test_that("reconstruct.metapca works correctly", {
  set.seed(123)
  X1 <- matrix(rnorm(20*10), 20, 10)
  X2 <- matrix(rnorm(20*15), 20, 15)
  
  pc1 <- pca(X1, ncomp=5, preproc=center())
  pc2 <- pca(X2, ncomp=5, preproc=center())
  
  mp <- metapca(list(pc1, pc2), ncomp=8)
  
  # Test full reconstruction
  recon_full <- reconstruct(mp)
  expect_equal(dim(recon_full), c(20, 25))  # 20 rows, 10+15 columns
  
  # Test partial reconstruction with subset of components
  recon_partial <- reconstruct(mp, comp=1:3)
  expect_equal(dim(recon_partial), c(20, 25))
  
  # Reconstruction with fewer components should have less or equal variance
  var_full <- sum(recon_full^2)
  var_partial <- sum(recon_partial^2)
  expect_lte(var_partial, var_full + 1e-10)  # Allow for numerical tolerance
  
  # Test reconstruction with subset of rows
  recon_subset <- reconstruct(mp, comp=1:3, rowind=1:5)
  expect_equal(nrow(recon_subset), 5)
  expect_equal(ncol(recon_subset), 25)
  
  # Test reconstruction with subset of columns
  recon_cols <- reconstruct(mp, comp=1:3, colind=1:10)
  expect_equal(dim(recon_cols), c(20, 10))
})

test_that("truncate.metapca preserves structure", {
  set.seed(123)
  X1 <- matrix(rnorm(20*10), 20, 10)
  X2 <- matrix(rnorm(20*15), 20, 15)
  
  pc1 <- pca(X1, ncomp=5, preproc=center())
  pc2 <- pca(X2, ncomp=5, preproc=center())
  
  mp <- metapca(list(pc1, pc2), ncomp=8)
  
  # Truncate to fewer components
  mp_trunc <- truncate(mp, 3)
  
  expect_s3_class(mp_trunc, "metapca")
  expect_equal(ncol(scores(mp_trunc)), 3)
  expect_equal(ncol(components(mp_trunc)), 3)
  expect_equal(length(sdev(mp_trunc)), 3)
  
  # Check that block indices are preserved
  expect_equal(mp_trunc$outer_block_indices, mp$outer_block_indices)
  expect_equal(mp_trunc$inner_block_indices, mp$inner_block_indices)
  
  # Check that fits are preserved
  expect_equal(length(mp_trunc$fits), length(mp$fits))
})

test_that("metapca handles scaled combine correctly", {
  set.seed(123)
  # Create blocks with very different variances
  X1 <- matrix(rnorm(20*10, sd=0.1), 20, 10)  # Low variance
  X2 <- matrix(rnorm(20*10, sd=10), 20, 10)   # High variance
  
  pc1 <- pca(X1, ncomp=5, preproc=center())
  pc2 <- pca(X2, ncomp=5, preproc=center())
  
  # Without scaling, high variance block dominates
  mp_unscaled <- metapca(list(pc1, pc2), ncomp=5, combine="pca")
  
  # With scaling, both blocks contribute more equally
  mp_scaled <- metapca(list(pc1, pc2), ncomp=5, combine="scaled")
  
  # The scaled version should have different scores
  scores_unscaled <- scores(mp_unscaled)
  scores_scaled <- scores(mp_scaled)
  
  # They should be different (not just sign flips)
  diff <- min(sum((scores_unscaled - scores_scaled)^2),
              sum((scores_unscaled + scores_scaled)^2))
  expect_gt(diff, 0.1)  # Substantial difference
})

test_that("metapca MFA combine works correctly", {
  set.seed(123)
  X1 <- matrix(rnorm(20*10), 20, 10)
  X2 <- matrix(rnorm(20*15), 20, 15)
  
  pc1 <- pca(X1, ncomp=5, preproc=center())
  pc2 <- pca(X2, ncomp=5, preproc=center())
  
  # Test MFA combining
  mp_mfa <- metapca(list(pc1, pc2), ncomp=5, combine="MFA")
  
  expect_s3_class(mp_mfa, "metapca")
  expect_equal(ncol(scores(mp_mfa)), 5)
  
  # MFA should weight by inverse of first singular value
  # Check that the metafit has the expected structure
  expect_true(!is.null(mp_mfa$metafit))
})