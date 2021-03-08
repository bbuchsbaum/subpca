#' subpca
#'
#' @param X the data matrix
#' @param clus the cluster index vector
#' @param ncomp the number of components to extract from each block (can vary by block or be fixed)
#' @param preproc the preprocessor
#' @param combine how method to combine cluster-wise PCAs.
#' @param extra args
#' @importFrom furrr future_map
#' @importFrom assertthat assert_that
#'
#' @examples
#'
#' X <- matrix(rnorm(10*100), 10, 100)
#' clus <- rep(1:10, each=10)
#' ncomp=10
#'
#' pres <- subpca(X, clus, ncomp=5)
subpca <- function(X, clus,
                      weights=NULL,
                      ncomp=2,
                      ccomp=2,
                      combine=c("pca", "scaled","MFA"),
                      preproc=center(), ...) {


  combine <- match.arg(combine)
  assert_that(length(clus) == ncol(X))
  assert_that(min(table(clus)) > 2)

  ngroups <- length(unique(clus))
  groups <- sort(unique(clus))

  if (length(ccomp) == 1) {
    ccomp <- rep(ccomp, ngroups)
  } else {
    assertthat::assert_that(length(ccomp) == ngroups)
  }

  if (!is.null(weights)) {
    assert_that(length(weights) == length(groups))
    assert_that(all(weights > 0))
    weights <- weights/sum(weights)
  } else {
    weights <- rep(1, ngroups)
  }

  proc <- prep(preproc)
  X <- init_transform(proc, X)

  sind <- split(1:ncol(X), clus)
  gsize <- sapply(sind, length)

  fits <- furrr::future_map(1:length(sind), function(i) {
      xb <- X[,sind[[i]]]
      pp <- multivarious::fresh(preproc)
      multivarious::pca(xb, ncomp=ccomp[i], preproc=pp)
  })

  # scmat <- do.call(cbind, furrr::future_map(fits, function(fit) {
  #   if (combine == "pca" || combine == "MFA") {
  #     multivarious::scores(fit)
  #   } else {
  #     fit$u
  #   }
  # }))
  #
  # A <- if (combine == "MFA") {
  #   sapply(1:length(fits), function(i) 1/sdev(fits[[i]])[1] * weights[i])
  # } else {
  #   weights
  # }
  #
  # A <- rep(A, sapply(fits, function(fit) ncomp(fit)))
  # final_fit <- genpca::genpca(scmat, A, preproc = pass())

  final_fit <- metapca(fits, ncomp=ncomp, combine=combine)
  attr(final_fit, "class") <- c("subpca", attr(final_fit, "class"))
  #sc <- scores(final_fit)
  #lds <- t(X) %*% sc

}

