
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


  assert_that(inherits(X, "multiblock_list"), msg="X must be a multiblock_list object")
  inner_combine <- match.arg(inner_combine)
  combine <- match.arg(combine)

  # Column-stacked: observations x variables, cluster the variables (columns)
  # Row-stacked: variables x observations, need to transpose for column-wise clustering
  if (!is_cstacked(X)) {
    X <- t(X)
  }

  assert_that(length(X) > 0, msg="X must contain at least one block")
  assert_that(all(vapply(X, is.matrix, logical(1))), msg="Each block in X must be a matrix")

  nvars <- vapply(X, ncol, integer(1))
  nobs <- vapply(X, nrow, integer(1))
  assert_that(length(unique(nvars)) == 1, msg="All blocks in X must have the same number of variables (columns)")
  assert_that(length(unique(nobs)) == 1, msg="All blocks in X must have the same number of observations (rows)")
  assert_that(length(clus) == nvars[[1]], msg="Length of 'clus' must equal number of variables per block")
  assert_that(!anyNA(clus), msg="`clus` must not contain NA values")
  assert_that(min(table(clus)) > 2, msg="Each cluster must have at least 3 variables")

  ngroups <- length(unique(clus))
  groups <- sort(unique(clus))
  nblocks <- length(X)
  dots <- list(...)

  reserved_args <- c("fits", "ncomp", "weights", "combine", "outer_block_indices")
  duplicate_args <- intersect(names(dots), reserved_args)
  if (length(duplicate_args) > 0) {
    stop("Do not pass ", paste(sprintf("`%s`", duplicate_args), collapse = ", "), " via `...`; these are controlled by `musubpca`")
  }

  cluster_indices <- split(seq_along(clus), factor(clus, levels=groups))
  outer_block_indices <- lapply(cluster_indices, function(ind) {
    unlist(lapply(seq_len(nblocks), function(b) ind + (b - 1L) * nvars[[1]]), use.names = FALSE)
  })

  resolve_ncomp <- function(ncomp, fit, i) {
    if (is.function(ncomp)) {
      nc <- if (length(formals(ncomp)) >= 2) ncomp(fit, i) else ncomp(fit)
    } else {
      nc <- ncomp
    }

    if (!is.numeric(nc) || length(nc) != 1 || !is.finite(nc)) {
      stop("`ccomp` function must return one finite numeric value")
    }

    nc <- floor(nc)

    if (nc < 1) {
      stop("`ccomp` function must return values >= 1")
    }

    min(nc, multivarious::ncomp(fit))
  }

  if (!is.function(ccomp)) {
    if (any(ccomp <= 0)) {
      stop("all `ccomp` elements must be greater than 0")
    }

    if (length(ccomp) == 1) {
      ccomp <- rep(ccomp, ngroups)
    } else {
      assertthat::assert_that(length(ccomp) == ngroups, msg="Length of 'ccomp' must equal number of groups")
    }

    if (any(ccomp < 1)) {
      if (any(ccomp > 1)) {
        stop("cannot mix fixed and fractional values of `ccomp`")
      }

      perc <- ccomp
      ccomp <- function(fit, i) {
        ev <- sdev(fit)^2
        v <- cumsum(ev) / sum(ev)
        which(v >= perc[i])[1]
      }
    }
  }

  ## for each block run a clusterpca
  out <- furrr::future_map(X, function(x) {
    clusterpca(x, clus, ccomp=inner_ccomp, preproc=multivarious::fresh(preproc))
  }, .options = furrr::furrr_options(packages = c("multivarious", "subpca")))


  ## for each cluster, run a metapca
  fits2 <- furrr::future_map(1:ngroups, function(i) {
    fs <- lapply(out, function(x) x$fits[[i]])
    pres <- if (is.function(ccomp)) {
      nc <- sum(sapply(fs, function(f) shape(f)[2]))
      fit0 <- metapca(fs, ncomp=nc, combine=inner_combine)
      n <- resolve_ncomp(ccomp, fit0, i)
      pres <- truncate(fit0, n)
    } else {
      ## combining fits for each cluster ...
      metapca(fs, ncomp=ccomp[i], combine=inner_combine)
    }
    ##...
  }, .options = furrr::furrr_options(packages = c("multivarious", "subpca")))

  pres <- metapca(
    fits2,
    ncomp=ncomp,
    weights=weights,
    combine=combine,
    outer_block_indices=outer_block_indices,
    ...
  )
  return(pres)
}
