
#' sub-block PCA
#'
#'
#' PCA on sublocks or clusters followed again by PCA
#'
#'
#'
#' @param X the data matrix
#' @param clus the cluster index vector
#' @param ncomp the number of final "meta" components.
#' @param ccomp the number of components to extract from each block (can vary by block or be fixed, see details).
#' @param preproc the preprocessor
#' @param combine how method to combine cluster-wise PCAs.
#' @param ... extra args sent to `metapca`
#' @importFrom furrr future_map
#' @importFrom assertthat assert_that
#'
#' @details
#'
#' the `ccomp` argument can either be a scalar integer, a fractional value less than 1,
#' a `vector` with the same length as the number of clusters, or a function taking a `pca` object
#' and returning an `integer`.
#'
#' @examples
#'
#' X <- matrix(rnorm(10*100), 10, 100)
#' clus <- rep(1:10, each=10)
#' ncomp=10
#'
#' pres <- subpca(X, clus, ncomp=5)
#'
#' f <- function(p) {
#'   ev <- sdev(p)^2
#'   v <- cumsum(ev) / sum(ev)
#'   min(which(v > .8))
#' }
#'
#' pres <- subpca(X, clus, ncomp=5, ccomp=f)
#' pres <- subpca(X, clus, ncomp=5, ccomp=.5)
#' @export
subpca <- function(X, clus,
                      weights=NULL,
                      ncomp=2,
                      ccomp=2,
                      combine=c("pca", "scaled","MFA"),
                      preproc=center(), ...) {


  combine <- match.arg(combine)
  assert_that(length(clus) == ncol(X), msg="Length of 'clus' must equal number of columns in X")
  assert_that(min(table(clus)) > 2, msg="Each cluster must have at least 3 observations")

  ngroups <- length(unique(clus))
  groups <- sort(unique(clus))

  if (!is.null(weights)) {
    assert_that(length(weights) == ngroups, msg="Length of 'weights' must equal number of groups")
  }


  cfit <- clusterpca(X, clus=clus, ccomp=ccomp, preproc=preproc)
  outer_block_indices <- split(1:ncol(X), clus)
  final_fit <- metapca(cfit$fits, ncomp=ncomp, combine=combine, weights=weights, outer_block_indices=outer_block_indices, ...)
  attr(final_fit, "class") <- c("subpca", attr(final_fit, "class"))
  attr(final_fit, "nclus") <- ngroups
  final_fit
}

#' @export
print.subpca <- function(x) {
  cat("subpca: ", paste0(class(x)), "\n")
  cat("number of clusters: ", attr(x, "nclus"), "\n")
  cat("input dim: ", nrow(x$v), "\n")
  cat("output dim: ", ncol(x$v), "\n")
}

#' @export
components.subpca <- function(x, ...) {
  x$v
}

