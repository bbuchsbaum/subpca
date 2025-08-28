#' Extract time filters (left singular vectors) for each cluster
#'
#' The "time filters" are the left singular vectors U from the SVD of Y_C = B_C^{-1/2} X_C.
#' They live in sample/time space (T x r).
#'
#' @param x A \code{time_contrast_clusterpca} object.
#' @param by_cluster If TRUE, return a list per cluster; otherwise, column-bind into T x sum(ncomp).
#' @return Either a list of T x r matrices (per cluster) or a single T x sum(ncomp) matrix.
#' @export
time_filters <- function(x, by_cluster = FALSE) UseMethod("time_filters")

#' @export
time_filters.time_contrast_clusterpca <- function(x, by_cluster = FALSE) {
  mats <- lapply(x$fits, function(f) f$u)
  if (isTRUE(by_cluster)) return(mats)
  if (!length(mats)) return(matrix(0, nrow = x$dims$T, ncol = 0L))
  do.call(cbind, mats)
}

#' Contrast "signal vs background" ratios per component
#'
#' For each cluster and component, returns the generalized-eigenvalue proxy
#' corresponding to the SVD of Y_C = B_C^{-1/2} X_C: \eqn{\lambda_i \propto d_i^2}.
#' If \code{normalize=TRUE}, divides by the cluster size |C| to approximate a covariance scaling.
#'
#' @param x A \code{time_contrast_clusterpca} object.
#' @param by_cluster If TRUE, return a list per cluster; otherwise a single numeric vector.
#' @param normalize If TRUE, report \eqn{d_i^2 / |C|}; else \eqn{d_i^2}.
#' @return Numeric vector (or list) of contrast ratios.
#' @export
contrast_ratio <- function(x, by_cluster = FALSE, normalize = TRUE) UseMethod("contrast_ratio")

#' @export
contrast_ratio.time_contrast_clusterpca <- function(x, by_cluster = FALSE, normalize = TRUE) {
  cr_list <- lapply(x$fits, function(f) {
    if (!f$ncomp) return(numeric(0))
    if (normalize) (f$d ^ 2) / max(1L, f$nvars) else (f$d ^ 2)
  })
  if (isTRUE(by_cluster)) return(cr_list)
  unlist(cr_list, use.names = FALSE)
}

#' Background diagnostics per cluster
#'
#' Quick metadata for the background whitener: rank (number of bg modes used),
#' ridge \code{tau}, and how many background clusters contributed.
#'
#' @param x A \code{time_contrast_clusterpca} object.
#' @return Data frame with one row per cluster.
#' @export
background_info <- function(x) {
  K <- length(x$fits)
  data.frame(
    cluster = vapply(x$fits, function(f) as.character(f$cluster), character(1)),
    bg_rank = vapply(x$fits, function(f) length(f$bg$s), integer(1)),
    tau     = vapply(x$fits, function(f) f$bg$tau, numeric(1)),
    used_bg_clusters = vapply(x$fits, function(f) length(f$bg$used), integer(1)),
    nvars   = vapply(x$fits, function(f) f$nvars, integer(1)),
    ncomp   = vapply(x$fits, function(f) f$ncomp, integer(1)),
    stringsAsFactors = FALSE
  )
}

# Optional: make predict() an alias for project()
#' @export
predict.time_contrast_clusterpca <- function(object, newdata, ...) {
  project(object, new_data = newdata, ...)
}