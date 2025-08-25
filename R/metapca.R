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
#' truncate(pfit,2)
#' @export
#' @import assertthat
metapca <- function(fits, ncomp=2, weights=NULL, combine=c("pca", "scaled","MFA"), outer_block_indices=NULL) {
  assert_that(all(sapply(fits,function(f) inherits(f, "bi_projector"))))

  combine <- match.arg(combine)

  X <- do.call(cbind, lapply(fits, function(x) scores(x)))
  nc <- sapply(fits, function(x) shape(x)[2])

  if (is.null(weights)) {
    weights <- rep(1, sum(nc))
  } else  {
    assert_that(length(weights) == length(fits), msg="weights must have one value per fit")
    weights <- rep(weights, nc)
  }

  pres <- if (combine == "scaled") {
    wts <- (1/sqrt(matrixStats::colSds(X))) * weights
    genpca::genpca(unclass(X), A=wts, ncomp=ncomp, preproc=pass())
  } else if (combine == "MFA") {
    wts <- rep(sapply(fits, function(x) 1/sdev(x)[1]), nc) * weights
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

  if (is.null(outer_block_indices)) {
    outer_block_indices <- b_ind(nvars)
  } else {
    vcount <- sapply(outer_block_indices, length)
    assert_that(all(nvars == vcount), msg="`outer_block_indices` does not correpsond to input shape of `fits`")
  }

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
      # Get positions of colind that are in this block
      idx <- which(colind %in% ind)
      nd <- new_data[,idx,drop=FALSE]
      # For individual block projection, use sequential indices 1:ncol(nd)
      # not the original colind values
      project(x$fits[[i]], nd)
    } else {
      matrix(0, nrow(new_data), ncomp(x$fits[[i]]))
    }
  }))

  project(x$metafit, x0)
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
truncate.metapca <- function(x, ncomp) {
  # Preserve original class hierarchy
  original_classes <- class(x)
  ret <- bi_projector(
    v=coef(x)[,1:ncomp,drop=FALSE],
    s=scores(x)[,1:ncomp,drop=FALSE],
    sdev=sdev(x)[1:ncomp],
    metafit=truncate(x$metafit, ncomp),
    fits=x$fits,
    outer_block_indices=x$outer_block_indices,
    inner_block_indices=x$inner_block_indices,
    classes=original_classes)
  
  # Preserve any extra attributes
  if ("nclus" %in% names(attributes(x))) {
    attr(ret, "nclus") <- attr(x, "nclus")
  }
  ret
}

#' @export
reconstruct.metapca <- function(x, comp=seq_len(ncomp(x)), rowind=1:nrow(scores(x)), colind=seq_len(sum(sapply(x$outer_block_indices, length))), ...) {
  # Check if this is a musubpca context by examining the fit structure
  # In musubpca, each fit is itself a metapca with clustered structure
  is_musubpca <- all(sapply(x$fits, function(f) inherits(f, "metapca")))
  
  if (is_musubpca) {
    # Special handling for musubpca: reconstruct each cluster separately 
    # and assemble back to original variable space
    cluster_recons <- lapply(x$fits, function(fit) {
      reconstruct(fit, rowind=rowind)
    })
    
    # Determine the number of clusters and blocks
    nclus <- length(cluster_recons)
    # Get number of blocks from outer_block_indices of first fit
    nblocks <- length(x$fits[[1]]$outer_block_indices)
    nvars_per_cluster <- ncol(cluster_recons[[1]]) / nblocks  # Each cluster spans nblocks
    
    # Assemble clusters back to original variable structure
    # Total variables = nclus * nvars_per_cluster * nblocks
    total_vars <- nclus * nvars_per_cluster * nblocks
    assembled_recon <- matrix(0, length(rowind), total_vars)
    
    for (i in seq_len(nclus)) {
      # For each block, place the cluster's variables in the right positions
      for (b in seq_len(nblocks)) {
        # Variables for cluster i in block b
        cluster_block_start <- (b-1) * nvars_per_cluster + 1
        cluster_block_end <- b * nvars_per_cluster
        
        # Position in assembled matrix for cluster i, block b
        assembled_start <- (b-1) * (nclus * nvars_per_cluster) + (i-1) * nvars_per_cluster + 1
        assembled_end <- assembled_start + nvars_per_cluster - 1
        
        assembled_recon[, assembled_start:assembled_end] <- 
          cluster_recons[[i]][, cluster_block_start:cluster_block_end]
      }
    }
    
    return(assembled_recon[, colind, drop=FALSE])
    
  } else {
    # Standard metapca reconstruction 
    # Reconstruct meta-level inner matrix (concatenated inner scores)
    S_hat <- reconstruct(x$metafit, comp=comp, rowind=rowind)
    
    # For each outer block, reconstruct original variables from the block's portion of S_hat
    blocks <- seq_along(x$outer_block_indices)
    mats <- lapply(blocks, function(i) {
      inner_idx <- x$inner_block_indices[[i]]
      # Get the block-specific portion of the reconstructed meta scores
      block_scores <- S_hat[, inner_idx, drop=FALSE]
      
      # Manually reconstruct by matrix multiplication: scores %*% t(components)
      # This gives us the reconstruction for just the specified rows
      block_components <- components(x$fits[[i]])
      block_recon <- block_scores %*% t(block_components)
      
      return(block_recon)
    })
    
    # Assemble result respecting colind
    full <- do.call(cbind, mats)
    return(full[, colind, drop=FALSE])
  }
}

#' @export
components.metapca <- function(x, ...) {
  x$v
}



