#' PCA of PCAs
#'
#' metapca computes a principal components analysis over several pca fits
#'
#' @param fits a list of fitted component models of type `bi_projector`
#' @param ncomp number of final components
#' @param weights option non-negative weighting vector -- 1 value per fit
#' @param combine method for weighting input pca fits
#'
#' X1 <- matrix(rnorm(20*10), 20, 10)
#' X2 <- matrix(rnorm(20*20), 20, 20)
#' X3 <- matrix(rnorm(20*30), 20, 30)
#'
#' pc1 <- pca(X1, ncomp=10, preproc=center())
#' pc2 <- pca(X2, ncomp=10, preproc=center())
#' pc3 <- pca(X3, ncomp=10, preproc=center())
#'
#'
#'
#' fits <- list(pc1,pc2,pc3)
#' pfit <- metapca(fits, ncomp=15)
#' ncol(scores(pfit)) == 10
#' @export
#' @import assertthat
metapca <- function(fits, ncomp=2, weights=NULL, combine=c("pca", "scaled","MFA")) {
  assert_that(all(sapply(fits,function(f) inherits(f, "bi_projector"))))

  combine <- match.arg(combine)

  X <- do.call(cbind, lapply(fits, function(x) scores(x)))
  nc <- sapply(fits, function(x) shape(x)[2])

  if (is.null(weights)) {
    weights <- rep(1, sum(nc))
  } else  {
    weights <- rep(weights, nc)
  }

  pres <- if (combine == "scaled") {
    wts <- (1/sqrt(matrixStats::colSds(X))) * weights
    genpca::genpca(unclass(X), A=wts, ncomp=ncomp, preproc=pass())
  } else if (combine == "MFA") {
    wts <- rep(sapply(fits, function(x) sdev(x)[1]), nc) * weights
    genpca::genpca(X, A=wts, ncomp=ncomp, preproc=pass())
  } else {
    genpca::genpca(X, A=weights, ncomp=ncomp, preproc=pass())
  }

  b_ind <- function(nv) {
    offsets <- cumsum(c(1, nv))
    lapply(1:length(nvars), function(i) {
      seq(offsets[i], offsets[i] + nv[i]-1)
    })
  }



  nvars <- sapply(fits, function(x) shape(x)[1])
  tot <- sum(nvars)

  outer_block_indices <- b_ind(nvars)
  inner_block_indices <- b_ind(nc)


  v <- do.call(rbind, lapply(1:length(fits), function(i) {
      coef(fits[[i]]) %*% coef(pres)[inner_block_indices[[i]],]
  }))

  ret <- bi_projector(
    v=v,
    s=scores(pres),
    sdev=sdev(pres),
    metafit=pres,
    fits=fits,
    outer_block_indices=outer_block_indices,
    inner_block_indices=inner_block_indices,
    classes=c("metapca", "pca"))
}

#' @export
project_block.metapca <- function(x, new_data, block, ...) {
  partial_project(x, new_data, colind=x$outer_block_indices[[block]])
}

#' @export
partial_project.metapca <- function(x, new_data, colind, ...) {
  x0 <- do.call(cbind, lapply(1:length(x$outer_block_indices), function(i) {
    ind <- x$outer_block_indices[[i]]
    keep <- colind %in% ind
    if (sum(keep) > 0) {
      idx <- which(colind %in% ind)
      nd <- new_data[,idx,drop=FALSE]
      project(x$fits[[i]], nd,colind=which(keep))
    } else {
      matrix(0, nrow(scores(x$fits[[i]])), x$fits[[i]]$ncomp)
    }
  }))

  project(x$metafit,x0, comp)
}


#' @export
project.metapca <- function(x, new_data) {
  #browser()
  assert_that(ncol(new_data) == sum(sapply(x$fits, function(f) shape(f)[1])))
  x0 <- do.call(cbind, lapply(1:length(x$outer_block_indices), function(i) {
    ind <- x$outer_block_indices[[i]]
    nd <- new_data[,ind,drop=FALSE]
    project(x$fits[[i]], nd)
  }))

  project(x$metafit,x0)

}

#' @export
reconstruct.metapca <- function(x, comp, rowind=1:nrow(scores(x)), colind=1:nrow(coefficients(x))) {

  recon <- reconstruct(x$metafit, comp=comp)


  ## map colind to blocks.
  ## reconstruct all blocks
  ## subset colind
  blocks <- seq_along(x$outer_block_indices)

  inblock <- sapply(1:length(x$outer_block_indices), function(i) {
    ind <- x$outer_block_indices[[i]]
    any(colind %in% ind)
  })

  blocks <- blocks[inblock]
  ## reconstruct here means reconstruct the original data
  ret <- lapply(blocks, function(i) {
    ## reconstruct the inner matrix
    ## recon <- reconstruct(x$metafit, newdata, comp=comp, colind=x$block_indices[[i]], reverse_pre_process=TRUE)
    ## reconstruct the outer matrix
    as.matrix(reconstruct(x$fits[[i]], new_data=recon[,x$inner_block_indices[[i]]]))
  })

  block_matrix(ret)

}



