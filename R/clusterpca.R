
truncpca <- function(fit, ccomp, i) {
  if (is.function(ccomp)) {
    n <- ccomp(fit,i)
    multivarious::truncate(fit, n)
  } else if (multivarious::ncomp(fit) != ccomp[i]) {
    nc <- min(ccomp[i], ncomp(fit))
    truncate(fit, nc)
  } else {
    fit
  }
}

#' clusterwise pca
#'
#' Fit a separate PCA model for each cluster or block of data.
#'
#' @inheritParams subpca
#'
#'
#' @param colwise clusters are columnwise if `TRUE` otherwise they are rowwise
#' @param pcafun an optional function that computes PCA fit that returns a `bi_projector` instance. see details.
#' @examples
#'
#' X <- matrix(rnorm(10*100), 10, 100)
#' clus <- rep(1:10, each=10)
#'
#'
#' pres <- clusterpca(X, clus, ccomp=5, colwise=TRUE)
#' cf <- coef(pres)
#'
#' X <- matrix(rnorm(200*100), 200, 100)
#' clus <- rep(1:10, each=20)
#' pres <- clusterpca(X, clus, ccomp=5, colwise=FALSE)
#' @import multivarious
#' @export
clusterpca <- function(X, clus,
                  ccomp=2,
                  preproc=center(),
                  colwise=TRUE,
                  pcafun=NULL) {

  if (colwise) {
    assert_that(length(clus) == ncol(X))
  } else {
    assert_that(length(clus) == nrow(X))
  }

  assert_that(min(table(clus)) > 2)

  ngroups <- length(unique(clus))
  groups <- sort(unique(clus))

  if (is.null(pcafun)) {
    pcafun <- function(X, ncomp, preproc, ind=NULL) {
      multivarious::pca(X, ncomp=ncomp, preproc=preproc)
    }
  }

  if (!is.function(ccomp)) {
    if (any(ccomp <= 0)) {
      stop("all `ccomp` elements must be greater than 0")
    }

    if (length(ccomp) == 1) {
      ccomp <- rep(ccomp, ngroups)
    } else {
      assertthat::assert_that(length(ccomp) == ngroups)
    }

    if (any(ccomp < 1)) {
      if (any(ccomp > 1)) {
        stop("cannot mix fixed and fractional values of `ccomp`")
      }

      if (length(ccomp) != ngroups) {
        ccomp <- rep(ccomp, length.out=ngroups)
      }
    }

    if (ccomp[1] < 1) {
      perc <- ccomp
      ccomp <- function(fit, i) {
        svar <- cumsum(sdev(fit))
        v <- svar/svar[length(svar)]
        min(which(v > perc[i]))
      }
    }
  }

  proc <- prep(preproc)
  X <- init_transform(proc, X)

  sind <- if (colwise) {
    split(1:ncol(X), clus)
  } else {
    split(1:nrow(X), clus)
  }

  gsize <- sapply(sind, length)

  fits <- furrr::future_map(1:length(sind), function(i) {
    xb <- if (colwise) X[,sind[[i]],drop=FALSE] else X[sind[[i]],,drop=FALSE]
    pp <- multivarious::fresh(preproc)
    nc <- if (is.function(ccomp)) min(dim(xb)) else ccomp[i]
    fit0 <- pcafun(xb, ncomp=nc, preproc=pp, sind[[i]])
    truncpca(fit0, ccomp,i)
  })

  nc <- sapply(fits,ncomp)
  cunc <- c(0, cumsum(nc))

  comp_indices <- lapply(1:(length(cunc)-1), function(i) {
    seq(cunc[i]+1,cunc[i+1])
  })

  ret <- list(fits=fits,
              ncomp=nc,
              comp_indices=comp_indices,
              preproc=proc,
              groups=groups,
              clus=clus,
              colwise=colwise,
              block_indices=sind)

  class(ret) <- c("clusterpca", "bi_projector", "projector")
  ret

}

#' @export
#' @importFrom multivarious scores
#' @importFrom Matrix sparseMatrix
coef.clusterpca <- function(object) {

  v <- if (object$colwise) {
    nv <- sapply(object$fits, function(f) shape(f)[2])
    offsets <- cumsum(c(1, nv))
    comp_indices <- lapply(1:length(nv), function(i) {
      seq(offsets[i], offsets[i] + nv[i]-1)
    })

    Reduce("+", lapply(seq_along(object$fits), function(i) {
      fit <- object$fits[[i]]
      v <- coef(fit)
      sc <- scores(fit)

      sparseMatrix(i=rep(object$block_indices[[i]], ncol(v)), j=rep(comp_indices[[i]], each=nrow(v)), x=as.vector(v),
                   dims=c(length(object$clus),sum(nv)))

    }))
  } else {
    do.call(cbind, lapply(object$fits, function(fit) {
      coef(fit)
    }))
  }

  v


}

#' @export
scores.clusterpca <- function(x, ...) {

  sc <- if (!x$colwise) {
    nv <- sapply(x$fits, function(f) shape(f)[2])
    offsets <- cumsum(c(1, nv))
    comp_indices <- lapply(1:length(nv), function(i) {
      seq(offsets[i], offsets[i] + nv[i]-1)
    })

    Reduce("+", lapply(seq_along(x$fits), function(i) {
      fit <- x$fits[[i]]
      v <- scores(fit)
      #sc <- scores(fit)

      sparseMatrix(i=rep(x$block_indices[[i]], ncol(v)), j=rep(comp_indices[[i]], each=nrow(v)),
                   x=as.vector(v),
                   dims=c(length(x$clus),sum(nv)))

    }))
  } else {
    do.call(cbind, lapply(x$fits, function(fit) {
      scores(fit)
    }))
  }

  sc


}

#' @export
sdev.clusterpca <- function(x) {
  sc <- scores(x)
  apply(sc,2,function(x) sqrt(sum(x^2)))
}

#' @export
ncomp.clusterpca <- function(x) {
  x$ncomp
}


#' @export
shape.clusterpca <- function(x) {
  v <- coef(x)
  c(nrow(v), ncol(v))
}

#' @export
residuals.clusterpca <- function(x, ncomp=x$ncomp, xorig, ...) {
  if (length(ncomp) == 1) {
    ncomp <- rep(ncomp, length(x$fits))
  }

  res <- lapply(seq_along(x$fits), function(i) {
    fit <- x$fits[[i]]
    if (x$colwise) {
      residuals(fit, ncomp=ncomp[i], xorig=xorig[,x$block_indices[[i]], drop=FALSE])
    } else {
      residuals(fit, ncomp=ncomp[i], xorig=xorig[x$block_indices[[i]], , drop=FALSE])
    }
  })

  xres <- matrix(0, nrow(xorig), ncol(xorig))

  for (i in seq_along(x$fits)) {
    if (x$colwise) {
      xres[,x$block_indices[[i]]] <- res[[i]]
    } else {
      xres[x$block_indices[[i]],] <- res[[i]]
    }
  }

  xres

}

