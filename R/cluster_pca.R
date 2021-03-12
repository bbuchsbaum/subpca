cluster_pca <- function(X, clus,
                  weights=NULL,
                  ccomp=2,
                  preproc=center()) {

  assert_that(length(clus) == ncol(X))
  assert_that(min(table(clus)) > 2)

  ngroups <- length(unique(clus))
  groups <- sort(unique(clus))

  if (!is.null(weights)) {
    assert_that(length(weights) == length(groups))
    assert_that(all(weights > 0))
    weights <- weights/sum(weights)
  } else {
    weights <- rep(1, ngroups)
  }

  if (!is.function(ccomp)) {
    if (any(ccomp <= 0)) {
      stop("`ccomp` must be greater than 0")
    }
    if (length(ccomp) == 1) {
      ccomp <- rep(ccomp, ngroups)
    } else {
      assertthat::assert_that(length(ccomp) == ngroups)
    }

    if (ccomp[1] < 1) {
      perc <- ccomp[1]
      ccomp <- function(fit) {
        svar <- cumsum(sdev(fit))
        v <- svar/svar[length(svar)]
        min(which(v > perc))
      }
    }
  }

  proc <- prep(preproc)
  X <- init_transform(proc, X)

  sind <- split(1:ncol(X), clus)
  gsize <- sapply(sind, length)

  fits <- furrr::future_map(1:length(sind), function(i) {
    xb <- X[,sind[[i]]]
    pp <- multivarious::fresh(preproc)

    if (is.function(ccomp)) {
      fit0 <- multivarious::pca(xb, ncomp=min(dim(xb)), preproc=pp)
      n <- ccomp(fit0)
      multivarious::truncate(fit0, n)
    } else {
      multivarious::pca(xb, ncomp=ccomp[i], preproc=pp)
    }
  })


  fits

}
