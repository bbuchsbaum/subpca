#' Time-contrast cluster PCA (sample-space generalized PCA with background)
#'
#' For each cluster of columns (voxels), compute a PCA in \emph{time/sample} space
#' that emphasizes variance in that cluster while downweighting variance captured
#' by a \emph{background} built from the \emph{mean time-series of other clusters}.
#' The background can be weighted (e.g., by spatial distance or similarity) and
#' is handled efficiently via a low-rank whitener:
#' \deqn{B_C = L_C L_C^\top + \tau I,\quad L_C = [\sqrt{w_k}(m_k - \mu^{(C)})]_{k\ne C}}
#' so that we whiten the data as \eqn{Y_C = (B_C)^{-1/2} X_C} and perform a thin SVD.
#'
#' @param X Numeric matrix (T x P): rows are timepoints/samples, columns are voxels/features.
#' @param clus Integer/factor of length P assigning each column to a cluster (K clusters).
#' @param ncomp Integer (>=1), number of components per cluster (can be a scalar or length-K vector).
#' @param scheme List describing how to weight background clusters. Basic fields:
#'   \itemize{
#'     \item \code{type}: one of \code{"uniform"}, \code{"gaussian"}, \code{"knn_near"}, \code{"knn_far"},
#'           \code{"rings"}, \code{"similarity"}.
#'     \item \code{sigma}: bandwidth for \code{"gaussian"} (on distances between clusters).
#'     \item \code{k}: number of neighbors for \code{"knn_*"}.
#'     \item \code{rings}: numeric vector of distance breakpoints for \code{"rings"}.
#'       Note: The method gives uniform weight to all background clusters in the
#'       outermost ring (most distant from target) and zero weight to all others.
#'       Example: \code{c(5, 15)} creates rings at distances 0-5, 5-15, >15.
#'     \item \code{phi}: for \code{"similarity"}, how to convert correlation to weights; one of
#'           \code{"abs"}, \code{"pos"}, \code{"neg"}.
#'   }
#'   Missing fields fall back to sensible defaults.
#' @param distmat Optional (K x K) matrix of nonnegative distances between clusters (e.g., centroids).
#'   Required for distance-based schemes (\code{"gaussian"}, \code{"knn_*"}, \code{"rings"}).
#' @param q_bg Optional integer; if set, truncate the background rank to \code{q_bg} (top modes of \eqn{L_C}).
#' @param tau Nonnegative ridge added to \eqn{B_C}; default is data-driven (small fraction of background energy).
#' @param center_time Logical; if TRUE (default) center each column (voxel) over time before fitting.
#' @param mode Engine: \code{"whiten"} (generalized PCA via background whitener),
#'   \code{"partial"} (residualize cluster on background means w/ ridge),
#'   or \code{"contrastive"} (cPCA: eigen of A_C - alpha B_C).
#' @param alpha Contrastive strength for \code{mode="contrastive"}.
#'   If \code{NULL}, uses \code{alpha_grid} to auto-select.
#' @param alpha_grid Optional numeric vector of candidate alphas for auto-selection (default: 10^seq(-2,2,length=7)).
#' @param ridge_bg Ridge parameter for background regression in \code{mode="partial"} (default: 1e-2).
#' @param svd_method SVD backend when decomposing \eqn{Y_C}: \code{"irlba"} (default) or \code{"svd"}.
#' @param verbose Logical; print per-cluster progress.
#'
#' @return An object of class \code{c("time_contrast_clusterpca","clusterpca")}, with:
#'   \itemize{
#'     \item \code{clus}: factor of cluster assignments for columns.
#'     \item \code{clusters}: cluster labels (levels).
#'     \item \code{indices}: list of column indices per cluster.
#'     \item \code{fits}: list per cluster with fields:
#'           \code{u} (T x r), \code{v} (|C| x r), \code{d} (length r),
#'           \code{bg} whitener parameters (\code{U}, \code{s}, \code{tau}),
#'           and \code{time_scores} (T x r).
#'     \item \code{scheme}, \code{q_bg}, \code{tau}, \code{ncomp}: config
#'     \item \code{dims}: list with \code{T}, \code{P}, \code{K}.
#'   }
#'
#' @details
#' The method constructs a low-rank background for each cluster C using the mean
#' time-series of the \emph{other} clusters. Weighting schemes let you emphasize
#' near vs. far clusters, rings/annuli, or time-course similarity.
#'
#' Projection to new data uses the stored whitener (\eqn{U,s,\tau}) and voxel
#' loadings (\code{v}) to compute time scores in the same contrastive subspace.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' T <- 200; P <- 1000; K <- 50
#' clus <- sample(1:K, P, replace=TRUE)
#' X <- matrix(rnorm(T*P), T, P)
#' # Fake cluster centroid distances (for demo only)
#' distmat <- as.matrix(stats::dist(matrix(rnorm(K*3), K, 3)))
#' fit <- time_contrast_clusterpca(
#'   X, clus, ncomp=3,
#'   scheme=list(type="gaussian", sigma=median(distmat)),
#'   distmat=distmat, q_bg=20, verbose=TRUE
#' )
#' dim(scores(fit))  # T x sum(ncomp)
#' dim(coef(fit))    # P x sum(ncomp)
#' }
#'
#' @export
time_contrast_clusterpca <- function(
  X, clus,
  ncomp = 3L,
  scheme = list(type = "uniform", sigma = NULL, k = 32L, rings = NULL, phi = "abs"),
  distmat = NULL,
  q_bg = NULL,
  tau = NULL,
  center_time = TRUE,
  mode = c("whiten", "partial", "contrastive"),
  alpha = NULL, alpha_grid = 10^seq(-2, 2, length.out = 7L),
  ridge_bg = 1e-2,
  svd_method = c("irlba", "svd"),
  verbose = FALSE
) {
  # Constants
  DENSE_EIG_FALLBACK_LIMIT <- 800L  # Max T for dense eigendecomposition in contrastive mode
  SMALL_MATRIX_THRESHOLD <- 50L     # Below this, use base::svd instead of irlba
  
  stopifnot(is.matrix(X), is.numeric(X))
  P <- ncol(X); Tn <- nrow(X)
  if (length(clus) != P) stop("length(clus) must equal ncol(X)")
  svd_method <- match.arg(svd_method)
  mode <- match.arg(mode)

  # standardize cluster ids
  clus <- as.factor(clus)
  clusters <- levels(clus)
  K <- length(clusters)
  indices <- split(seq_len(P), clus)

  # center over time if requested
  if (isTRUE(center_time)) {
    X <- sweep(X, 2L, colMeans(X), FUN = "-")
  }

  # compute cluster mean time-series: M (T x K)
  M <- .tc_cluster_means(X, indices)
  # precompute time-centering of M just to be safe (already centered by X centering)
  # but do not re-center columns of M independently, to keep alignment with X.

  # allow scalar or vector ncomp
  if (length(ncomp) == 1L) ncomp <- rep.int(as.integer(ncomp), K)
  if (length(ncomp) != K) stop("ncomp must be length 1 or K")

  # distance-based schemes need distmat
  need_dist <- isTRUE(scheme$type %in% c("gaussian", "knn_near", "knn_far", "rings"))
  if (need_dist && (is.null(distmat) || !is.matrix(distmat) || nrow(distmat) != K || ncol(distmat) != K)) {
    stop("distmat (K x K) is required for scheme type '", scheme$type, "'.")
  }

  # default tau: small fraction of mean background variance (computed per cluster below if NULL)
  auto_tau <- is.null(tau)
  
  # Check for irlba availability once
  use_irlba <- svd_method == "irlba" && requireNamespace("irlba", quietly = TRUE)

  fits <- vector("list", K)

  for (ci in seq_len(K)) {
    C_name <- clusters[[ci]]
    C_idx  <- indices[[ci]]
    # build weights over other clusters
    w <- .tc_weights_for_cluster(
      ci, M, distmat, scheme = scheme
    )
    # low-rank background factors L_C = [sqrt(w_k) (m_k - mu)]_{k != C}
    bg <- .tc_build_background(M, w, ci, q_bg = q_bg)

    # choose tau per cluster if auto (used in mode="whiten")
    tau_ci <- if (auto_tau) .tc_default_tau(bg) else tau

    Xc <- X[, C_idx, drop = FALSE]
    r  <- min(ncomp[[ci]], min(nrow(Xc), ncol(Xc)))
    if (r == 0L) {
      fits[[ci]] <- .tc_empty_fit(Tn, length(C_idx), bg, tau_ci, C_name)
      next
    }

    # --- compute per-mode representation ---
    fit_result <- switch(
      mode,
      whiten = .tc_fit_cluster_whiten(Xc, r, bg, tau_ci, use_irlba, SMALL_MATRIX_THRESHOLD),
      partial = .tc_fit_cluster_partial(Xc, r, bg, ridge_bg, use_irlba, SMALL_MATRIX_THRESHOLD),
      contrastive = .tc_fit_cluster_contrastive(Xc, r, bg, alpha, alpha_grid, tau_ci, C_name, Tn, C_idx, DENSE_EIG_FALLBACK_LIMIT)
    )
    
    # Handle NULL result from contrastive mode (empty fit)
    if (is.null(fit_result)) {
      next  # Empty fit was already stored by helper function
    }
    
    u <- fit_result$u
    v <- fit_result$v
    d <- fit_result$d
    time_scores <- fit_result$time_scores
    bg_mode <- fit_result$bg_mode

    fits[[ci]] <- list(
      cluster = C_name,
      ncomp   = r,
      u = u, v = v, d = d,
      time_scores = time_scores,
      bg = bg_mode,
      nobs = Tn, nvars = length(C_idx), col_index = C_idx
    )

    if (isTRUE(verbose)) {
      msg <- sprintf("[time_contrast] cluster %s: |C|=%d, r=%d, bg_rank=%d, tau=%.3e",
                     C_name, length(C_idx), r, length(bg$s), tau_ci)
      message(msg)
    }
  }

  structure(list(
    clus = clus,
    clusters = clusters,
    indices = indices,
    fits = fits,
    scheme = scheme,
    q_bg = q_bg,
    tau = tau,
    ncomp = ncomp,
    dims = list(T = Tn, P = P, K = K),
    center_time = center_time,
    mode = mode,
    alpha = alpha,
    ridge_bg = ridge_bg,
    svd_method = svd_method,
    version = "0.2.0"
  ), class = c("time_contrast_clusterpca", "clusterpca"))
}

# --- Helpers ---------------------------------------------------------------

# Compute cluster means M (T x K)
.tc_cluster_means <- function(X, indices) {
  Tn <- nrow(X); K <- length(indices)
  M <- matrix(0, nrow = Tn, ncol = K)
  for (k in seq_len(K)) {
    idx <- indices[[k]]
    if (length(idx)) {
      # mean across columns in the cluster
      M[, k] <- rowMeans(X[, idx, drop = FALSE])
    } else {
      M[, k] <- 0
    }
  }
  colnames(M) <- as.character(seq_len(K))
  M
}

# Weight builder wrapper (returns length-K vector with w[C]=0)
.tc_weights_for_cluster <- function(ci, M, distmat, scheme) {
  K <- ncol(M)
  type <- scheme$type %||% "uniform"
  w <- switch(
    type,
    graphdiff = {
      W <- scheme$W
      if (is.null(W) || !is.matrix(W) || nrow(W) != ncol(W) || nrow(W) != ncol(M)) {
        stop("scheme$W must be a K x K diffusion weight matrix for type 'graphdiff'.")
      }
      as.numeric(W[ci, ])
    },
    uniform = .tc_w_uniform(K, ci),
    gaussian = .tc_w_gaussian(K, ci, distmat, sigma = scheme$sigma),
    knn_near = .tc_w_knn(K, ci, distmat, k = scheme$k %||% 32L, near = TRUE),
    knn_far  = .tc_w_knn(K, ci, distmat, k = scheme$k %||% 32L, near = FALSE),
    rings    = .tc_w_rings(K, ci, distmat, rings = scheme$rings),
    similarity = .tc_w_similarity(K, ci, M, phi = scheme$phi %||% "abs"),
    stop("Unknown scheme$type: ", type)
  )
  w[ci] <- 0
  # normalize if any positive weights remain
  if (sum(w) > 0) w <- w / sum(w)
  w
}

# Uniform weights on complement
.tc_w_uniform <- function(K, ci) {
  w <- rep(1 / (K - 1), K)
  w[ci] <- 0
  w
}

# Gaussian of distance
.tc_w_gaussian <- function(K, ci, distmat, sigma) {
  if (is.null(sigma) || !is.finite(sigma) || sigma <= 0)
    sigma <- median(distmat[distmat > 0], na.rm = TRUE)
  d <- distmat[ci, ]
  w <- exp(- (d * d) / (2 * sigma * sigma))
  w[ci] <- 0
  w
}

# k-nearest or k-farthest neighbors
.tc_w_knn <- function(K, ci, distmat, k = 32L, near = TRUE) {
  d <- distmat[ci, ]
  ord <- order(d, decreasing = !near)
  ord <- ord[ord != ci]
  take <- head(ord, n = min(k, K - 1L))
  w <- numeric(K); w[take] <- 1
  w
}

# Rings (annuli) – simple on/off using provided breakpoints
# Weights only the outermost ring (most distant clusters) to emphasize 
# local patterns that differ from the global context
.tc_w_rings <- function(K, ci, distmat, rings) {
  if (is.null(rings) || !length(rings)) stop("scheme$rings must be provided for type 'rings'.")
  d <- distmat[ci, ]
  # put mass on the outermost ring by default
  cuts <- sort(unique(rings))
  # assign to bins [0, cuts1], (cuts1, cuts2], ..., (cuts_m, Inf)
  bin <- findInterval(d, vec = cuts, rightmost.closed = FALSE, all.inside = FALSE)
  # give weight 1 to the outermost bin (largest index); others 0
  w <- as.numeric(bin == max(bin) & seq_along(bin) != ci)
  w
}

# Similarity by mean time-course correlation
.tc_w_similarity <- function(K, ci, M, phi = c("abs", "pos", "neg")) {
  phi <- match.arg(phi)
  mc <- M[, ci]
  cval <- suppressWarnings(stats::cor(mc, M))
  cval[is.na(cval)] <- 0
  w <- switch(
    phi,
    abs = abs(cval),
    pos = pmax(0, cval),
    neg = pmax(0, -cval)
  )
  w[ci] <- 0
  w
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Build low-rank background L_C and return its SVD pieces + metadata
.tc_build_background <- function(M, w, ci, q_bg = NULL) {
  # weighted mean of background cluster means
  if (sum(w) <= 0) {
    # no background selected
    return(list(U = matrix(0, nrow = nrow(M), ncol = 0L),
                s = numeric(0), mu = rep(0, nrow(M)), used = integer(0), w = w))
  }
  mu <- drop(M %*% w) # weighted average across columns
  used <- which(w > 0 & seq_len(ncol(M)) != ci)
  if (!length(used)) {
    return(list(U = matrix(0, nrow = nrow(M), ncol = 0L),
                s = numeric(0), mu = mu, used = integer(0), w = w))
  }
  # Form L = [sqrt(w_k) (m_k - mu)]_{k in used}
  W12 <- sqrt(w[used])
  L <- sweep(M[, used, drop = FALSE], 1L, mu, FUN = "-")
  L <- sweep(L, 2L, W12, FUN = "*")
  # Truncated SVD if requested
  rmax <- min(nrow(M), ncol(L))
  if (!is.null(q_bg)) rmax <- min(rmax, as.integer(q_bg))
  if (rmax == 0L) {
    return(list(U = matrix(0, nrow = nrow(M), ncol = 0L),
                s = numeric(0), mu = mu, used = used, w = w))
  }
  # Use SVD fallback for small matrices to avoid irlba warnings
  min_dim <- min(nrow(L), ncol(L))
  svdL <- if (requireNamespace("irlba", quietly = TRUE) && min_dim > 50L && rmax < min_dim * 0.5) {
    irlba::irlba(L, nu = rmax, nv = 0L)
  } else {
    # base svd, compute only U and d
    s <- base::svd(L, nu = rmax, nv = 0L)
    list(u = s$u, d = s$d)
  }
  list(U = svdL$u, s = svdL$d, mu = mu, used = used, w = w)
}

# Default tau based on background energy
.tc_default_tau <- function(bg) {
  if (length(bg$s) == 0L) return(1e-3)
  # fraction of average eigenvalue of L L^T (i.e., mean(s^2)/T)
  ev_mean <- mean(bg$s^2)
  tau <- max(1e-8, 1e-3 * ev_mean / max(1L, nrow(bg$U)))
  tau
}

# Apply (B)^{-1/2} to a matrix Z using U, s for L and ridge tau: B = L L^T + tau I
.tc_whiten_apply <- function(U, s, tau, Z) {
  if (!length(s) || ncol(U) == 0) {
    # B = tau I -> scale
    return(Z / sqrt(tau))
  }
  # Ensure U and s dimensions are consistent
  r <- min(ncol(U), length(s))
  if (r == 0) {
    return(Z / sqrt(tau))
  }
  U <- U[, seq_len(r), drop = FALSE]
  s <- s[seq_len(r)]
  
  # Z_white = (1/sqrt(tau)) Z + U * diag(1/sqrt(s^2+tau) - 1/sqrt(tau)) * (U' Z)
  a <- 1 / sqrt(tau)
  delta <- (1 / sqrt(s * s + tau)) - a
  UtZ <- crossprod(U, Z)               # r x q
  UtZ <- sweep(UtZ, 1L, delta, "*")    # row-scale by delta
  Zw  <- a * Z + U %*% UtZ
  Zw
}

# Fallback SVD returning u, v, d of rank r (thin)
.tc_svd <- function(Y, r) {
  s <- base::svd(Y, nu = r, nv = r)
  list(u = s$u, v = s$v, d = s$d[seq_len(r)])
}

# Compute SVD with automatic fallback for small matrices
.tc_compute_svd <- function(M, r, use_irlba, small_threshold = 50L) {
  # Use base::svd for small matrices to avoid irlba warnings
  min_dim <- min(nrow(M), ncol(M))
  if (use_irlba && min_dim > small_threshold && r < min_dim * 0.5) {
    # Use irlba only for larger matrices where we want few components
    irlba::irlba(M, nu = r, nv = r)
  } else {
    # Fall back to base::svd for small matrices or when irlba unavailable
    .tc_svd(M, r)
  }
}

# Helper function for whiten mode computation
.tc_fit_cluster_whiten <- function(Xc, r, bg, tau_ci, use_irlba, small_threshold) {
  # Y = (LL^T + tau I)^(-1/2) X_C
  Y <- .tc_whiten_apply(bg$U, bg$s, tau_ci, Xc)
  svd_res <- .tc_compute_svd(Y, r, use_irlba, small_threshold)
  
  list(
    u = svd_res$u,
    v = svd_res$v,
    d = svd_res$d,
    time_scores = svd_res$u %*% diag(svd_res$d, nrow = length(svd_res$d), ncol = length(svd_res$d)),
    bg_mode = list(mode = "whiten", U = bg$U, s = bg$s, tau = tau_ci, mu = bg$mu, used = bg$used, w = bg$w)
  )
}

# Helper function for partial mode computation  
.tc_fit_cluster_partial <- function(Xc, r, bg, ridge_bg, use_irlba, small_threshold) {
  # R = (I - H_ridge) X_C with H = U diag(s^2/(s^2+rho)) U^T
  R <- .tc_residualize_apply(bg$U, bg$s, ridge_bg, Xc)
  svd_res <- .tc_compute_svd(R, r, use_irlba, small_threshold)
  
  list(
    u = svd_res$u,
    v = svd_res$v,
    d = svd_res$d,
    time_scores = svd_res$u %*% diag(svd_res$d, nrow = length(svd_res$d), ncol = length(svd_res$d)),
    bg_mode = list(mode = "partial", U = bg$U, s = bg$s, rho = ridge_bg, mu = bg$mu, used = bg$used, w = bg$w)
  )
}

# Helper function for contrastive mode computation
.tc_fit_cluster_contrastive <- function(Xc, r, bg, alpha, alpha_grid, tau_ci, C_name, Tn, C_idx, dense_limit) {
  # Find leading eigenvectors of Delta = A_C - alpha * (L L^T)
  alpha_ci <- if (is.null(alpha)) {
    .tc_alpha_select(Xc, bg$U, bg$s, alpha_grid, r_try = max(3L, r), dense_limit = dense_limit)
  } else {
    alpha
  }
  
  eig <- .tc_contrastive_eigs(Xc, bg$U, bg$s, alpha_ci, r, dense_limit = dense_limit)
  
  if (eig$r == 0L) {
    # Return empty fit with contrastive metadata
    empty_fit <- .tc_empty_fit(Tn, length(C_idx), bg, tau_ci, C_name)
    empty_fit$bg$mode <- "contrastive"
    empty_fit$bg$alpha <- alpha_ci
    return(NULL)  # Signal empty result
  }
  
  u <- eig$vec    # T x r_pos
  lam <- eig$val  # positive eigenvalues
  
  # Map to SVD-like triplet for consistency
  Vu <- crossprod(Xc, u)  # |C| x r
  d <- sqrt(colSums(Vu * Vu))
  keep <- which(d > 0)
  
  if (!length(keep)) {
    empty_fit <- .tc_empty_fit(Tn, length(C_idx), bg, tau_ci, C_name)
    empty_fit$bg$mode <- "contrastive"
    empty_fit$bg$alpha <- alpha_ci
    return(NULL)  # Signal empty result
  }
  
  u <- u[, keep, drop = FALSE]
  v <- sweep(Vu[, keep, drop = FALSE], 2L, d[keep], "/")
  d <- d[keep]
  
  list(
    u = u,
    v = v,
    d = d,
    time_scores = Xc %*% v,  # project onto eigenvectors
    bg_mode = list(mode = "contrastive", U = bg$U, s = bg$s, alpha = alpha_ci, mu = bg$mu, used = bg$used, w = bg$w),
    empty_metadata = list(alpha = alpha_ci)  # For storing in empty fits
  )
}

# Residualize on background means with ridge: R = (I - H_ridge) Z, H = U diag(s^2/(s^2+rho)) U^T
.tc_residualize_apply <- function(U, s, rho, Z) {
  if (!length(s)) return(Z)
  if (!is.finite(rho) || rho < 0) rho <- 0
  # Z_res = Z - U diag( s^2/(s^2+rho) ) U^T Z  = U diag(rho/(s^2+rho)) U^T Z + (I - U U^T) Z
  UtZ <- crossprod(U, Z)                           # r x q
  scale <- (s * s) / pmax(s * s + rho, .Machine$double.eps)
  Z_hat <- U %*% (scale * UtZ)                    # projection onto bg subspace (ridge-shrunk)
  Z - Z_hat
}

# Leading positive eigenpairs of Delta = (1/|C|) X X^T - alpha * (L L^T); returns up to k positives
.tc_contrastive_eigs <- function(Xc, U, s, alpha, k, dense_limit = 800L) {
  Tn <- nrow(Xc)
  # Afun * v = (X X^T v)/|C| - alpha * U diag(s^2) U^T v
  Afun <- function(v, args) {
    v <- as.numeric(v)
    xv <- Xc %*% (crossprod(Xc, v) / ncol(Xc))
    if (length(s)) {
      Utv <- crossprod(U, v)
      bgv <- U %*% ((s * s) * Utv)
      xv - alpha * bgv
    } else {
      xv
    }
  }
  # Try RSpectra if available; otherwise form small dense Delta when T is small
  r_try <- min(k + 5L, Tn - 1L)
  if (requireNamespace("RSpectra", quietly = TRUE)) {
    es <- try(RSpectra::eigs_sym(Afun, k = r_try, which = "LA", n = Tn), silent = TRUE)
    if (inherits(es, "try-error")) {
      # fallback dense for small T
      if (Tn <= dense_limit) {
        Delta <- (Xc %*% t(Xc)) / ncol(Xc)
        if (length(s)) Delta <- Delta - alpha * (U %*% ( (s * s) * t(U) ))
        ev <- eigen(Delta, symmetric = TRUE)
        lam <- ev$values; vec <- ev$vectors
      } else {
        stop("RSpectra failed and T is large; install RSpectra for mode='contrastive'.")
      }
    } else {
      lam <- es$values; vec <- es$vectors
    }
  } else {
    if (Tn > dense_limit) stop("RSpectra not installed; cannot efficiently compute contrastive eigen for large T.")
    Delta <- (Xc %*% t(Xc)) / ncol(Xc)
    if (length(s)) Delta <- Delta - alpha * (U %*% ( (s * s) * t(U) ))
    ev <- eigen(Delta, symmetric = TRUE)
    lam <- ev$values; vec <- ev$vectors
  }
  # keep positive eigenvalues (contrastive directions)
  pos <- which(lam > 0)
  if (!length(pos)) return(list(r = 0L, val = numeric(0), vec = matrix(0, nrow = Tn, ncol = 0)))
  ord <- pos[order(lam[pos], decreasing = TRUE)]
  take <- head(ord, k)
  list(r = length(take), val = lam[take], vec = vec[, take, drop = FALSE])
}

# Simple alpha auto-selection: pick alpha maximizing sum of top-r_try positive eigenvalues
.tc_alpha_select <- function(Xc, U, s, alpha_grid, r_try = 3L, dense_limit = 800L) {
  best <- alpha_grid[1]; best_score <- -Inf
  for (a in alpha_grid) {
    eig <- .tc_contrastive_eigs(Xc, U, s, a, k = r_try, dense_limit = dense_limit)
    scr <- if (eig$r > 0) sum(pmax(eig$val, 0)) else -Inf
    if (scr > best_score) { best_score <- scr; best <- a }
  }
  best
}

# Empty fit when r=0 (edge case)
.tc_empty_fit <- function(Tn, Pk, bg, tau_ci, C_name) {
  list(
    cluster = C_name,
    ncomp = 0L,
    u = matrix(0, nrow = Tn, ncol = 0L),
    v = matrix(0, nrow = Pk, ncol = 0L),
    d = numeric(0),
    time_scores = matrix(0, nrow = Tn, ncol = 0L),
    bg = list(U = bg$U, s = bg$s, tau = tau_ci, mu = bg$mu, used = bg$used, w = bg$w),
    nobs = Tn, nvars = Pk, col_index = integer(0)
  )
}

# --- S3 methods ------------------------------------------------------------

#' @export
scores.time_contrast_clusterpca <- function(x, ...) {
  # T x sum(ncomp) : concatenate per-cluster time scores
  mats <- lapply(x$fits, function(f) f$time_scores)
  if (!length(mats)) return(matrix(0, nrow = x$dims$T, ncol = 0L))
  do.call(cbind, mats)
}

#' @export
coef.time_contrast_clusterpca <- function(object, ...) {
  # P x sum(ncomp) block-sparse voxel loadings
  K <- length(object$fits)
  ncols <- sum(vapply(object$fits, function(f) f$ncomp, integer(1L)))
  if (ncols == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(object$dims$P, 0L)))
  }
  i_list <- vector("list", K)
  j_list <- vector("list", K)
  x_list <- vector("list", K)

  col_offset <- 0L
  for (k in seq_len(K)) {
    f <- object$fits[[k]]
    r <- f$ncomp
    if (r == 0L) next
    idx <- f$col_index
    # fill block: rows are original voxel indices, columns are global component indices
    i_block <- rep(idx, times = r)
    j_block <- rep(seq.int(col_offset + 1L, col_offset + r), each = length(idx))
    x_block <- as.vector(f$v) # column-major: matches (|C| x r)

    i_list[[k]] <- i_block
    j_list[[k]] <- j_block
    x_list[[k]] <- x_block

    col_offset <- col_offset + r
  }
  i <- unlist(i_list, use.names = FALSE)
  j <- unlist(j_list, use.names = FALSE)
  x <- unlist(x_list, use.names = FALSE)
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(object$dims$P, ncols))
}

#' @export
sdev.time_contrast_clusterpca <- function(x, ...) {
  # singular values are column norms of time scores
  unlist(lapply(x$fits, function(f) f$d), use.names = FALSE)
}

#' Number of components per cluster or in total
#' @export
ncomp.time_contrast_clusterpca <- function(x, by_cluster = FALSE, ...) {
  nc <- vapply(x$fits, function(f) f$ncomp, integer(1L))
  if (isTRUE(by_cluster)) nc else sum(nc)
}

#' Project new data into the time-contrast subspaces
#'
#' @param object A \code{time_contrast_clusterpca} fit.
#' @param new_data Numeric matrix (T_new x P_new). If \code{colind} is missing,
#'   columns must align to the original P; otherwise, \code{colind} maps the columns
#'   of \code{new_data} to the original feature indices (length P_new integer).
#' @param colind Optional integer vector giving original column indices corresponding
#'   to \code{new_data}'s columns; useful when projecting subsets/reordered columns.
#' @param center_time If TRUE, center each column over time before projection (default: TRUE).
#' @return Matrix of size \code{T_new x sum(ncomp)} with time scores.
#' @export
project.time_contrast_clusterpca <- function(object, new_data, colind = NULL, center_time = TRUE, ...) {
  stopifnot(is.matrix(new_data), is.numeric(new_data))
  Xn <- new_data
  if (isTRUE(center_time)) Xn <- sweep(Xn, 2L, colMeans(Xn), FUN = "-")

  Tn <- nrow(Xn)
  out_list <- vector("list", length(object$fits))

  col_offset <- 0L
  for (k in seq_along(object$fits)) {
    f <- object$fits[[k]]
    r <- f$ncomp
    if (r == 0L) { out_list[[k]] <- matrix(0, nrow = Tn, ncol = 0L); next }

    idx <- f$col_index
    if (!is.null(colind)) {
      pos <- match(idx, colind, nomatch = 0L)
      keep <- pos[pos > 0L]
      if (!length(keep)) {
        out_list[[k]] <- matrix(0, nrow = Tn, ncol = r)
        next
      }
      Xc <- Xn[, keep, drop = FALSE]
      # If some cluster columns are missing, we still project with the available ones.
    } else {
      # assume alignment
      Xc <- Xn[, idx, drop = FALSE]
    }

    # per-mode transform
    if (is.null(f$bg$mode) || f$bg$mode == "whiten") {
      Ynew <- .tc_whiten_apply(f$bg$U, f$bg$s, f$bg$tau, Xc)
    } else if (f$bg$mode == "partial") {
      Ynew <- .tc_residualize_apply(f$bg$U, f$bg$s, f$bg$rho, Xc)
    } else { # "contrastive" uses identity transform
      Ynew <- Xc
    }
    # project onto voxel loadings: scores = Y * V
    # (T x |C|) %*% (|C| x r) => T x r
    scores_k <- Ynew %*% f$v
    out_list[[k]] <- scores_k
    col_offset <- col_offset + r
  }
  do.call(cbind, out_list)
}

#' @export
print.time_contrast_clusterpca <- function(x, ...) {
  cat("time_contrast_clusterpca\n")
  cat(sprintf("  T = %d, P = %d, K = %d\n", x$dims$T, x$dims$P, x$dims$K))
  cat(sprintf("  mode = %s\n", x$mode))
  nc <- vapply(x$fits, function(f) f$ncomp, integer(1L))
  cat(sprintf("  total components = %d (mean per cluster = %.2f)\n", sum(nc), mean(nc)))
  cat("  scheme: "); str(x$scheme, no.list = TRUE)
  invisible(x)
}