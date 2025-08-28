#' Cluster-centroid distance matrix from voxel coordinates
#'
#' @param coords Numeric matrix (P x D) of voxel coordinates (e.g., mm in MNI or surface coords).
#' @param clus Factor/integer of length P assigning each voxel to a cluster (K clusters).
#' @param metric One of \code{"euclidean"}, \code{"manhattan"}, or \code{"cosine"} (1 - cosine similarity).
#' @return K x K distance matrix aligned to \code{levels(as.factor(clus))}.
#' @export
cluster_centroid_distmat <- function(coords, clus, metric = c("euclidean","manhattan","cosine")) {
  stopifnot(is.matrix(coords), nrow(coords) == length(clus))
  metric <- match.arg(metric)
  clus <- as.factor(clus)
  K <- nlevels(clus)
  idx <- split(seq_len(nrow(coords)), clus)
  # centroids: K x D
  Ctr <- vapply(idx, function(ii) colMeans(coords[ii, , drop = FALSE]), numeric(ncol(coords)))
  Ctr <- t(Ctr) # K x D
  rownames(Ctr) <- levels(clus)

  if (metric %in% c("euclidean","manhattan")) {
    d <- as.matrix(stats::dist(Ctr, method = metric))
    rownames(d) <- colnames(d) <- levels(clus)
    return(d)
  } else {
    # cosine distance: 1 - (x·y)/(|x||y|)
    nr <- sqrt(rowSums(Ctr * Ctr))
    Cn <- Ctr / pmax(nr, .Machine$double.eps)
    sim <- Cn %*% t(Cn)
    d <- pmax(0, 1 - sim)
    rownames(d) <- colnames(d) <- levels(clus)
    return(d)
  }
}

#' Heat-kernel diffusion weights from a cluster distance matrix
#'
#' Builds a dense K x K diffusion weight matrix W ≈ exp(-beta * L), where
#' L is the graph Laplacian built from a Gaussian KNN affinity on the provided distmat.
#' For large K, consider constructing W elsewhere and passing it via scheme$W to avoid O(K^3).
#'
#' @param distmat K x K distances between clusters (nonnegative, symmetric).
#' @param k Number of neighbors to keep in the affinity (default 20).
#' @param sigma Gaussian kernel bandwidth (default: median of nonzero distances).
#' @param beta Diffusion time/scale (default 1).
#' @return K x K diffusion weights matrix (row-stochastic).
#' @export
cluster_diffusion_weights <- function(distmat, k = 20L, sigma = NULL, beta = 1) {
  stopifnot(is.matrix(distmat), nrow(distmat) == ncol(distmat))
  K <- nrow(distmat)
  if (is.null(sigma) || !is.finite(sigma) || sigma <= 0)
    sigma <- stats::median(distmat[distmat > 0], na.rm = TRUE)
  # Build KNN Gaussian affinity
  A <- matrix(0, K, K)
  for (i in seq_len(K)) {
    d <- distmat[i, ]
    ord <- order(d)
    nn <- setdiff(ord[seq_len(min(k + 1L, K))], i) # drop self
    A[i, nn] <- exp(-(d[nn]^2) / (2 * sigma * sigma))
  }
  # Symmetrize
  A <- pmax(A, t(A))
  # Laplacian L = D - A
  dvec <- rowSums(A)
  L <- diag(dvec) - A
  # Dense eigen (small/moderate K). For large K, construct W externally.
  ev <- eigen(L, symmetric = TRUE, only.values = FALSE)
  U <- ev$vectors
  lam <- pmax(ev$values, 0)
  W <- U %*% diag(exp(-beta * lam), K, K) %*% t(U)
  # Row-normalize to get weights
  rs <- rowSums(W)
  rs[rs == 0] <- 1
  W <- W / rs
  W
}