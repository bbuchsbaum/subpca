
#' subpca applied to multiple data blocks with the same variables and clusters
#'
#' multiple subpca
#'
#' @inheritParams subpca
#' @importFrom multidesign is_cstacked
#' @examples
#'
#' X <- multidesign::multiblock(replicate(5, matrix(rnorm(50*50), 50, 50), simplify=FALSE))
#' clus <- rep(1:5, length.out=50)
#' ret <- musubpca(X,clus)
#' ret <- musubpca(X,clus, inner_ccomp=.5, ccomp=.5)
#' @export
musubpca <- function(X,
                    clus,
                    weights=NULL,
                    ncomp=2,
                    inner_ccomp=2,
                    ccomp=2,
                    inner_combine=c("pca", "scaled","MFA"),
                    combine=c("pca", "scaled","MFA"),
                    preproc=center(), ...) {


  assert_that(inherits(X, "multiblock_list"))
  inner_combine <- match.arg(inner_combine)
  combine <- match.arg(combine)

  # Column-stacked: observations x variables, cluster the variables (columns)
  # Row-stacked: variables x observations, need to transpose for column-wise clustering
  if (!is_cstacked(X)) {
    X <- t(X)
  }

  assert_that(min(table(clus)) > 2)
  assert_that(length(clus) == ncol(X[[1]]))

  ngroups <- length(unique(clus))
  groups <- sort(unique(clus))

  if (!is.function(ccomp) && length(ccomp) == 1) {
    ccomp <- rep(ccomp, ngroups)
  } else {
    assertthat::assert_that(length(ccomp) == ngroups)
  }

  if (ccomp[1] < 1) {
    perc <- ccomp[1]
    ccomp <- function(fit) {
      ev <- sdev(fit)^2
      v <- cumsum(ev) / sum(ev)
      min(which(v > perc))
    }
  }

  ## for each block run a cluster_pca
  out <- furrr::future_map(X, function(x) {
    clusterpca(x, clus, ccomp=inner_ccomp, preproc=multivarious::fresh(preproc))
  }, .options = furrr::furrr_options(packages = c("multivarious", "subpca")))


  ## for each cluster, run a metapca
  nfits <- length(out)

  fits2 <- furrr::future_map(1:ngroups, function(i) {
    fs <- lapply(out, function(x) x$fits[[i]])
    pres <- if (is.function(ccomp)) {
      nc <- sum(sapply(fs, function(f) shape(f)[2]))
      fit0 <- metapca(fs, ncomp=nc, combine=inner_combine)
      n <- ccomp(fit0)
      pres <- truncate(fit0, n)
    } else {
      ## combining fits for each cluster ...
      metapca(fs, ncomp=ccomp[i], combine=inner_combine)
    }
    ##...
  }, .options = furrr::furrr_options(packages = c("multivarious", "subpca")))

  pres <- metapca(fits2, ncomp=ncomp, weights=weights, combine=combine)
  return(pres)
}


