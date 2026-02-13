# 100 + 32*8 + 128*4 + 256*2 + 520

## orthogonalize
## M=t(p1$s) %*% p2$s
## o2s=(p2$s %*% corpcor::pseudoinverse(M))
## crossprod(p1$s, o2s)

## need set of sparse matrices, 1 per tree level.
## for every cluster, get the indices and the matrix one level up.
## orthogonalize with respect to submatrix at posiiton one level up.


#' Hierarchical PCA
#'
#' A multiscale PCA that defines sub-scales using a provided hierarchical clustering.
#'
#' @inheritParams clusterpca
#' @importFrom dendextend cutree
#' @examples
#'
#' grid <- expand.grid(1:20, 1:10)
#'
#' ## 100 images each with 50 features
#' X <- matrix(rnorm(200*100), 200, 100)
#' cuts <- c(4, 8, 16)
#' hclus <- hclust(dist(grid), method = "ward.D2")
#' hres1 <- hcluspca(X, hclus, cuts, est_method="standard", ccomp=c(4,2,2,2))
#' hres2 <- hcluspca(X, hclus, cuts, skip_global=TRUE, est_method="standard", ccomp=c(1,1,1))
#'
#' f <- function(fit,i) {  max(1, round(log(shape(fit)[1]))) }
#' hres3 <- hcluspca(X, hclus, cuts, skip_global=TRUE, est_method="standard", ccomp=f)
#'
#'
#' ncomp(hres1) == (sum(cuts) +4)
#'
#' @importFrom dendextend cutree
#' @export
hcluspca <- function(X, hclus, cuts,
                             est_method=c("standard","smooth"),
                             skip_global=FALSE,
                             ccomp=1,
                             spat_smooth=rep(0, length(cuts)+1),
                             cds=NULL,
                             intercept=FALSE,
                             svd_method=c("fast", "base", "irlba", "propack", "rsvd", "svds"),
                             orthogonalize=FALSE,
                             preproc=center()) {

  est_method <- match.arg(est_method)
  svd_method <- match.arg(svd_method)
  #preproc <- pre_processor(X, center=center, scale=scale)
  #Xp <- pre_process(preproc, X)

  nlevs <- if (skip_global) length(cuts) else length(cuts) + 1
  assertthat::assert_that(length(spat_smooth) >= nlevs, msg="`spat_smooth` must have at least one entry per fitted level")
  assertthat::assert_that(is.logical(intercept) && length(intercept) == 1 && !is.na(intercept), msg="`intercept` must be TRUE or FALSE")

  if (is.function(ccomp) && length(ccomp) == 1) {
    ccomp <- lapply(1:nlevs, function(i) ccomp)
  }

  if (!is.function(ccomp) && length(ccomp) == 1) {
    ccomp <- rep(ccomp, nlevs)
  }

  assertthat::assert_that(all(cuts > 1), msg="all `cuts` must be greater than 1")
  assertthat::assert_that(all(diff(cuts) > 0), msg="`cuts` must be strictly increasing from coarse to fine levels")

  if (!is.function(ccomp) && length(ccomp) != nlevs) {
    stop("if `ccomp` is a vector it must have an entry for every level, including the global level if skip_global=FALSE")
  }

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

  pfun_by_level <- function(level) {
    force(level)
    function(X, ncomp, preproc, ind) {
      if (est_method == "standard") {
        # If ncomp is a function, we need to fit with max components first, then truncate
        if (is.function(ncomp)) {
          message("fitting pca, method ", svd_method, " ncomp = function")
          fit0 <- pca(X, ncomp=min(dim(X)), preproc=preproc, method=svd_method)
          nc <- resolve_ncomp(ncomp, fit0, level)
          truncate(fit0, nc)
        } else {
          message("fitting pca, method ", svd_method, " ncomp = ", ncomp)
          pca(X, ncomp=ncomp, preproc=preproc, method=svd_method)
        }
      } else if (est_method == "smooth") {
        stopifnot(!is.null(cds))
        S <- spatial_constraints(cds[ind,,drop=FALSE], sigma_within=spat_smooth[level])
        S <- neighborweights::make_doubly_stochastic(S)
        genpca(X, M=S, ncomp=ncomp, preproc=preproc)
      } else {
        stop("unknown est_method")
      }
    }
  }

  fits <- vector(nlevs, mode="list")
  cutset <- lapply(cuts, function(ci) dendextend::cutree(hclus, ci))
  #kind <- dendextend::cutree(hclus, cuts[i])

  if (!skip_global) {
    ## outer fit
    pfun <- pfun_by_level(1)
    fit0 <- pfun(X, ccomp[[1]], preproc, 1:nrow(X))

    Xresid0 <- residuals(fit0, ncomp=multivarious::ncomp(fit0), xorig=X)
    Xresid <- Xresid0

    fits[[1]] <- fit0
    fi <- 1
    proc <- fit0$preproc
    S <- scores(fit0)
  } else {
    ## skip outer fit
    fit_preproc <- multivarious::fit_transform(preproc, X)
    proc <- fit_preproc$preproc
    Xresid0 <- fit_preproc$transformed
    Xresid <- Xresid0
    fi <- 0
    S <- NULL
  }

  for (i in seq_along(cuts)) {
    #print(i)
    message("residuals for level", i, " = ", sum(Xresid^2))
    kind <- cutset[[i]]
    #kind <- dendextend::cutree(hclus, cuts[i])
    message('fitting clusterpca, level ', i)
    ## no intercept...
    if (orthogonalize && !is.null(S)) {
      sind <- split(1:length(kind), kind)
      for (j in 1:length(sind)) {
        Z <- S[sind[[j]],,drop=FALSE]
        keep <- which(colSums(abs(Z)) > 0)
        if (length(keep) == 0) {
          next
        }

        Z <- as.matrix(Z[,keep,drop=FALSE])
        if (intercept) {
          Z <- cbind(`(Intercept)`=1, Z)
        }
        Y <- Xresid[sind[[j]],]
        rfit <- lm.fit(as.matrix(Z),as.matrix(Y))
        Xresid[sind[[j]],] <- resid(rfit)
      }
    }

    fit <- tryCatch(
      clusterpca(Xresid, clus = kind,
                 ccomp=ccomp[[fi+i]],
                 preproc=pass(),
                 colwise=FALSE, pcafun=pfun_by_level(fi+i)),
      error = function(e) {
        stop(
          sprintf("hcluspca failed at level %d (cut=%d): %s", i, cuts[i], conditionMessage(e)),
          call. = FALSE
        )
      }
    )


    message('computing residuals, level ', i)

    Xresid <- residuals.clusterpca(fit, ncomp=fit$ncomp, xorig=Xresid)
    fits[[fi+i]] <- fit

    if (is.null(S)) {
      S <- multivarious::scores(fit)
    } else {
      S <- cbind(S, multivarious::scores(fit))
    }

    message("score matrix has ", ncol(S), " columns.")
  }


  v <- do.call(cbind, lapply(fits, multivarious::components))
  s <- lapply(fits, scores)


  # if (orthogonalize) {
  #   for (i in 2:length(fits)) {
  #     message("i:", i)
  #     fit_i <- fits[[i]]
  #     for (j in seq_along(fit_i$block_indices)) {
  #       message("j:", j)
  #
  #       bind <- fit_i$block_indices[[j]]
  #       cind <- fit_i$comp_indices[[j]]
  #
  #       print(bind)
  #       print(cind)
  #       s_i <- s[[i]][bind,cind]
  #       s_parent <- s[[i-1]][bind,]
  #
  #
  #       M <- Matrix::crossprod(s_parent, s_i)
  #
  #
  #       s_io <- (s_i %*% corpcor::pseudoinverse(M))
  #       s[[i]][bind,cind] <- s_io
  #     }
  #   }
  # }

  s <- as.matrix(do.call(cbind, s))

  bi_projector(v,s, sdev=unlist(lapply(fits,sdev)), preproc=proc, classes="hcluspca",
               cuts=cuts, levels=length(cuts), hclus=hclus, spat_smooth=spat_smooth, cds=cds,
               ccomp=ccomp,
               svd_method=svd_method)


  # u_final <- do.call(cbind, fits)
  # reg <- regress(u_final, X, method="linear", intercept=TRUE)
  # class(reg) <- c("hpca", class(reg))
  # reg$levels <- length(cuts)
  # reg$hclus <- hclus
  # reg$cuts <- cuts
  # reg$spat_smooth <- spat_smooth
  # reg
}
