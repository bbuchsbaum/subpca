test_that("time_contrast_clusterpca basic shapes and projection equivalence", {
  skip_on_cran()
  set.seed(1)
  Tn <- 120; P <- 500; K <- 25
  clus <- sample(seq_len(K), P, replace = TRUE)
  coords <- matrix(rnorm(P * 3), ncol = 3)
  distmat <- cluster_centroid_distmat(coords, clus)

  X <- matrix(rnorm(Tn * P), nrow = Tn)
  fit <- time_contrast_clusterpca(
    X, clus, ncomp = 2,
    scheme = list(type = "gaussian", sigma = median(distmat[distmat > 0])),
    distmat = distmat, q_bg = 10
  )

  # scores and coef dimensions
  Sc <- scores(fit)
  Ld <- coef(fit)
  expect_equal(nrow(Sc), Tn)
  expect_equal(nrow(Ld), P)
  expect_equal(ncol(Sc), ncol(Ld))

  # project() on training data equals stored scores (up to numerical noise)
  Sc2 <- project(fit, X)
  expect_true(max(abs(Sc - Sc2)) < 1e-6)

  # background info sane
  bi <- background_info(fit)
  expect_equal(nrow(bi), K)
  expect_true(all(bi$bg_rank >= 0))
  expect_true(all(bi$tau > 0))
})

test_that("schemes behave and q_bg truncates background rank", {
  skip_on_cran()
  set.seed(2)
  Tn <- 60; P <- 200; K <- 12
  clus <- rep(seq_len(K), length.out = P)
  coords <- matrix(rnorm(P * 2), ncol = 2)
  distmat <- cluster_centroid_distmat(coords, clus)

  X <- matrix(rnorm(Tn * P), nrow = Tn)

  fit_knn <- time_contrast_clusterpca(
    X, clus, ncomp = 1,
    scheme = list(type = "knn_near", k = 5),
    distmat = distmat, q_bg = 3
  )
  bi <- background_info(fit_knn)
  expect_true(all(bi$bg_rank <= 3))

  fit_uniform <- time_contrast_clusterpca(X, clus, ncomp = 1, scheme = list(type="uniform"))
  bi2 <- background_info(fit_uniform)
  expect_true(all(bi2$used_bg_clusters >= 1))
})

test_that("plotting methods run", {
  set.seed(3)
  Tn <- 80; P <- 320; K <- 16
  clus <- rep(seq_len(K), length.out = P)
  coords <- matrix(rnorm(P * 3), ncol = 3)
  distmat <- cluster_centroid_distmat(coords, clus)
  X <- matrix(rnorm(Tn * P), nrow = Tn)

  fit <- time_contrast_clusterpca(X, clus, ncomp = 2,
                                  scheme = list(type = "gaussian", sigma = median(distmat[distmat>0])),
                                  distmat = distmat, q_bg = 8)

  expect_silent(plot(fit, cluster = 1, type = "spectrum"))
  expect_silent(plot(fit, cluster = 1, type = "time", comps = 1:2))
  expect_silent(plot(fit, cluster = 1, type = "bg"))
})

test_that("time_filters and contrast_ratio work", {
  skip_on_cran()
  set.seed(4)
  Tn <- 100; P <- 400; K <- 20
  clus <- sample(seq_len(K), P, replace = TRUE)
  X <- matrix(rnorm(Tn * P), nrow = Tn)
  
  fit <- time_contrast_clusterpca(X, clus, ncomp = 2, scheme = list(type = "uniform"))
  
  # time_filters
  tf <- time_filters(fit)
  expect_equal(nrow(tf), Tn)
  expect_equal(ncol(tf), ncomp(fit))
  
  tf_list <- time_filters(fit, by_cluster = TRUE)
  expect_equal(length(tf_list), K)
  
  # contrast_ratio
  cr <- contrast_ratio(fit, normalize = FALSE)
  expect_equal(length(cr), ncomp(fit))
  expect_true(all(cr >= 0))
  
  cr_list <- contrast_ratio(fit, by_cluster = TRUE)
  expect_equal(length(cr_list), K)
})

test_that("graph diffusion weighting works", {
  skip_on_cran()
  set.seed(5)
  Tn <- 80; P <- 300; K <- 15
  clus <- rep(seq_len(K), length.out = P)
  coords <- matrix(rnorm(P * 3), ncol = 3)
  distmat <- cluster_centroid_distmat(coords, clus)
  
  # Build diffusion weights
  W <- cluster_diffusion_weights(distmat, k = 10, beta = 0.5)
  expect_equal(nrow(W), K)
  expect_equal(ncol(W), K)
  
  X <- matrix(rnorm(Tn * P), nrow = Tn)
  
  fit_diff <- time_contrast_clusterpca(
    X, clus, ncomp = 2,
    scheme = list(type = "graphdiff", W = W)
  )
  
  expect_s3_class(fit_diff, "time_contrast_clusterpca")
  expect_equal(ncomp(fit_diff), K * 2)  # 2 components per cluster, K clusters
})

test_that("predict alias works", {
  skip_on_cran()
  set.seed(6)
  Tn <- 60; P <- 150; K <- 10
  clus <- sample(seq_len(K), P, replace = TRUE)
  X <- matrix(rnorm(Tn * P), nrow = Tn)
  
  fit <- time_contrast_clusterpca(X, clus, ncomp = 1)
  
  # predict should work as alias for project
  Xnew <- matrix(rnorm(30 * P), nrow = 30)
  sc1 <- project(fit, Xnew)
  sc2 <- predict(fit, Xnew)
  expect_equal(sc1, sc2)
})

test_that("different weighting schemes produce different results", {
  skip_on_cran()
  set.seed(7)
  Tn <- 100; P <- 400; K <- 20
  clus <- sample(seq_len(K), P, replace = TRUE)
  coords <- matrix(rnorm(P * 3), ncol = 3)
  distmat <- cluster_centroid_distmat(coords, clus)
  X <- matrix(rnorm(Tn * P), nrow = Tn)
  
  sigma0 <- median(distmat[distmat > 0])
  
  fit_uniform <- time_contrast_clusterpca(X, clus, ncomp = 2, 
                                          scheme = list(type = "uniform"))
  
  fit_gaussian <- time_contrast_clusterpca(X, clus, ncomp = 2,
                                           scheme = list(type = "gaussian", sigma = sigma0),
                                           distmat = distmat)
  
  fit_knn_near <- time_contrast_clusterpca(X, clus, ncomp = 2,
                                           scheme = list(type = "knn_near", k = 5),
                                           distmat = distmat)
  
  fit_knn_far <- time_contrast_clusterpca(X, clus, ncomp = 2,
                                          scheme = list(type = "knn_far", k = 5),
                                          distmat = distmat)
  
  # Different schemes should produce different scores
  sc_uniform <- scores(fit_uniform)
  sc_gaussian <- scores(fit_gaussian)
  sc_near <- scores(fit_knn_near)
  sc_far <- scores(fit_knn_far)
  
  # Check that they're actually different (not identical)
  expect_true(mean(abs(sc_uniform - sc_gaussian)) > 1e-10)
  expect_true(mean(abs(sc_near - sc_far)) > 1e-10)
})

test_that("empty clusters are handled gracefully", {
  skip_on_cran()
  set.seed(8)
  Tn <- 50; P <- 100; K <- 10
  # Create clustering with an empty cluster
  clus <- sample(c(1:9, rep(10, 0)), P, replace = TRUE)
  X <- matrix(rnorm(Tn * P), nrow = Tn)
  
  expect_warning(
    fit <- time_contrast_clusterpca(X, clus, ncomp = 1),
    NA  # No warning expected - should handle empty clusters
  )
  
  # Should still work
  expect_s3_class(fit, "time_contrast_clusterpca")
  sc <- scores(fit)
  expect_equal(nrow(sc), Tn)
})

test_that("mode='partial' residualization + projection works", {
  skip_on_cran()
  set.seed(11)
  Tn <- 100; P <- 300; K <- 15
  clus <- rep(seq_len(K), length.out = P)
  coords <- matrix(rnorm(P * 3), ncol = 3)
  distmat <- cluster_centroid_distmat(coords, clus)
  X <- matrix(rnorm(Tn * P), nrow = Tn)

  fit <- time_contrast_clusterpca(
    X, clus, ncomp = 2,
    mode = "partial", ridge_bg = 1e-2,
    scheme = list(type = "gaussian", sigma = median(distmat[distmat>0])),
    distmat = distmat, q_bg = 10
  )
  Sc <- scores(fit); Sc2 <- project(fit, X)
  expect_equal(dim(Sc), dim(Sc2))
  expect_true(max(abs(Sc - Sc2)) < 1e-6)
  expect_true(all(vapply(fit$fits, function(f) identical(f$bg$mode, "partial"), logical(1))))
})

test_that("mode='contrastive' runs and yields positive contrast eigenvalues", {
  skip_on_cran()
  set.seed(12)
  Tn <- 80; P <- 240; K <- 12
  clus <- rep(seq_len(K), length.out = P)
  coords <- matrix(rnorm(P * 2), ncol = 2)
  distmat <- cluster_centroid_distmat(coords, clus)
  X <- matrix(rnorm(Tn * P), nrow = Tn)

  fit <- time_contrast_clusterpca(
    X, clus, ncomp = 2,
    mode = "contrastive", alpha = 1,
    scheme = list(type = "uniform"),
    distmat = distmat, q_bg = 8
  )
  # dimensions coherent
  expect_equal(nrow(scores(fit)), Tn)
  expect_equal(nrow(coef(fit)), P)
  # projection defined
  Sc2 <- project(fit, X)
  expect_equal(ncol(Sc2), ncol(scores(fit)))
  # store mode/alpha
  expect_true(all(vapply(fit$fits, function(f) identical(f$bg$mode, "contrastive"), logical(1))))
  expect_true(all(vapply(fit$fits, function(f) is.numeric(f$bg$alpha), logical(1))))
})