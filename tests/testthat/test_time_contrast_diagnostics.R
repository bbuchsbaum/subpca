test_that("whiten mode with q_bg=0 and tau=1 reproduces plain time-PCA per cluster", {
  skip_on_cran()
  set.seed(101)
  Tn <- 90; P <- 360; K <- 18
  clus <- sample(seq_len(K), P, replace = TRUE)
  X <- matrix(rnorm(Tn * P), nrow = Tn)

  # Fit with no background modes and tau=1  -> Y = X (centered)
  fit <- time_contrast_clusterpca(
    X, clus, ncomp = 2,
    scheme = list(type = "uniform"),
    q_bg = 0, tau = 1,               # <- key: no bg, identity whitener
    center_time = TRUE,
    mode = "whiten",
    svd_method = "svd"               # deterministic for test
  )

  # Center X the same way the fitter does
  Xc <- sweep(X, 2L, colMeans(X), "-")

  # For every cluster, compare (U D) and V with base svd(X_cluster)
  for (k in seq_along(fit$fits)) {
    f <- fit$fits[[k]]
    r <- f$ncomp
    if (r == 0L) next

    idx <- f$col_index
    S <- svd(Xc[, idx, drop = FALSE], nu = r, nv = r)

    # Compare time scores (U D); allow sign flips per component
    UD_ref <- S$u %*% diag(S$d[seq_len(r)], r, r)
    TS     <- f$time_scores
    # align signs
    sgn <- sign(colSums(UD_ref * TS)); sgn[sgn == 0] <- 1
    TS_aligned <- sweep(TS, 2L, sgn, "*")
    expect_true(max(abs(UD_ref - TS_aligned)) < 1e-6)

    # Compare right singular vectors V (voxel loadings), up to sign
    V_ref <- S$v
    V     <- f$v
    diagcorr <- diag(abs(t(V_ref) %*% V))
    expect_true(min(diagcorr) > 1 - 1e-6)
  }
})

test_that("contrastive mode with alpha=0 reproduces plain time-PCA per cluster", {
  skip_on_cran()
  set.seed(202)
  Tn <- 80; P <- 320; K <- 16
  clus <- sample(seq_len(K), P, replace = TRUE)
  X <- matrix(rnorm(Tn * P), nrow = Tn)

  fit <- time_contrast_clusterpca(
    X, clus, ncomp = 2,
    scheme = list(type = "uniform"),
    q_bg = 10,                       # any (unused when alpha=0 in contrastive engine)
    alpha = 0,                       # <- key: pure PCA on X_C
    center_time = TRUE,
    mode = "contrastive"
  )

  Xc <- sweep(X, 2L, colMeans(X), "-")

  for (k in seq_along(fit$fits)) {
    f <- fit$fits[[k]]
    r <- f$ncomp
    if (r == 0L) next

    idx <- f$col_index
    S <- svd(Xc[, idx, drop = FALSE], nu = r, nv = r)

    # Compare (U D) with time_scores (sign-aligned)
    UD_ref <- S$u %*% diag(S$d[seq_len(r)], r, r)
    TS     <- f$time_scores
    sgn <- sign(colSums(UD_ref * TS)); sgn[sgn == 0] <- 1
    TS_aligned <- sweep(TS, 2L, sgn, "*")
    expect_true(max(abs(UD_ref - TS_aligned)) < 1e-6)

    # Compare V (right singulars)
    V_ref <- S$v
    V     <- f$v
    diagcorr <- diag(abs(t(V_ref) %*% V))
    expect_true(min(diagcorr) > 1 - 1e-6)
  }
})

test_that("partial mode with ridge=0 annihilates background subspace (orthogonality check)", {
  skip_on_cran()
  set.seed(303)
  Tn <- 100; P <- 400; K <- 20
  clus <- sample(seq_len(K), P, replace = TRUE)
  X <- matrix(rnorm(Tn * P), nrow = Tn)

  fit <- time_contrast_clusterpca(
    X, clus, ncomp = 2,
    scheme = list(type = "uniform"),
    q_bg = 15,                       # ensure a nontrivial bg subspace
    center_time = TRUE,
    mode = "partial",
    ridge_bg = 0,                    # <- key: exact projection (I - U U^T)
    svd_method = "svd"
  )

  # For each cluster, U^T * (time_scores) ≈ 0
  for (k in seq_along(fit$fits)) {
    f <- fit$fits[[k]]
    if (f$ncomp == 0L) next
    U  <- f$bg$U
    TS <- f$time_scores
    if (ncol(U) == 0L) next          # empty bg is vacuously orthogonal
    ortho <- max(abs(crossprod(U, TS)))
    expect_true(ortho < 1e-6)
  }
})