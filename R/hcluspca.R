# 100 + 32*8 + 128*4 + 256*2 + 520

#' Hierarchical PCA
#'
#' A multiscale PCA that defines sub-scales using a provided hierarchical clustering.
#'
#' @inheritParams clusterpca
#' @importFrom dendextend cutree
#' @examples
#'
#' grid <- expand.grid(1:10, 1:10)
#'
#' ## 100 images each with 50 features
#' X <- matrix(rnorm(100*50), 100, 50)
#' cuts <- c(4, 8, 16)
#' hclus <- dclust::dclust(grid, nstart=10)
#' hres1 <- hcluspca(X, hclus, cuts, est_method="standard", ccomp=c(4,1,1,1))
#' hres2 <- hcluspca(X, hclus, cuts, skip_global=TRUE, est_method="standard", ccomp=c(1,1,1))
#' ncomp(hres1) == (sum(cuts) +4)
#'
#' @importFrom dendextend cutree
#' @import dclust
#' @export
hcluspca <- function(X, hclus, cuts,
                             est_method=c("standard","smooth"),
                             skip_global=FALSE,
                             ccomp=1,
                             spat_smooth=rep(0, length(cuts)+1),
                             cds=NULL,
                             intercept=FALSE,
                             svd_method=c("fast", "base", "irlba", "propack", "rsvd", "svds"),
                             preproc=center()) {

  est_method <- match.arg(est_method)
  svd_method <- match.arg(svd_method)
  #preproc <- pre_processor(X, center=center, scale=scale)
  #Xp <- pre_process(preproc, X)

  nlevs <- if (skip_global) length(cuts) else length(cuts) + 1


  if (!is.function(ccomp) && length(ccomp) == 1) {
    ccomp <- rep(ccomp, nlevs)
  }

  assertthat::assert_that(all(cuts > 1), msg="all `cuts` must be greater than 1")

  if (!is.function(ccomp) && length(ccomp) != nlevs) {
    stop("if `ccomp` is a vector it must have an entry for every level, including the global level if skip_global=FALSE")
  }

  pfun <- function(X, ncomp, preproc, ind) {
    if (est_method == "standard") {
      pca(X, ncomp=ncomp, preproc=preproc, method=svd_method)
    } else if (est_method == "smooth") {
      coord <- cds[ind,]
      S <- spatial_constraints(cds[ind,,drop=FALSE], sigma_within=spat_smooth[level])
      S <- neighborweights::make_doubly_stochastic(S)
      genpca(x, M=S, ncomp=ncomp, preproc=preproc)
    } else {
      stop()
    }

  }

  if (!skip_global) {
    ## outer fit
    fit0 <- pfun(X, ccomp[[1]], preproc, 1:nrow(X))


    Xresid0 <- residuals(fit0, ncomp=multivarious::ncomp(fit0), xorig=X)
    Xresid <- Xresid0

    fits <- vector(nlevs, mode="list")
    fits[[1]] <- fit0
    fi <- 1
  } else {
    ## skip outer fit
    proc <- prep(preproc)
    Xresid0 <- init_transform(proc, X)
    Xresid <- Xresid0
    fi <- 0
  }

  for (i in 1:length(cuts)) {
    #print(i)
    message("residuals for level", i, " = ", sum(Xresid^2))
    kind <- dendextend::cutree(hclus, cuts[i])
    ## no intercept...
    fit <- clusterpca(Xresid,clus = kind, ccomp=ccomp[[fi+i]], preproc=pass(), colwise=FALSE, pcafun=pfun)
    Xresid <- residuals.clusterpca(fit, ncomp=ncomp(fit),xorig=Xresid)
    fits[[fi+i]] <- fit
  }

  v <- do.call(cbind, lapply(fits,coef))
  s <- do.call(cbind, lapply(fits, scores))

  bi_projector(v,s, sdev=unlist(lapply(fits,sdev)), preproc=fit0$preproc, classes="hcluspca",
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

