#' Plot diagnostics for time-contrast cluster PCA
#'
#' @param x A time_contrast_clusterpca object.
#' @param cluster Cluster to show (name or index).
#' @param type One of "spectrum", "time", "bg".
#'   \itemize{
#'     \item \strong{"spectrum"}: barplot of d^2 (contrast ratio proxy).
#'     \item \strong{"time"}: overlay of selected time filters (U columns).
#'     \item \strong{"bg"}: barplot of background weights w^{(C)} across clusters.
#'   }
#' @param comps Integer vector of component indices for type="time".
#' @param reorder_bg If TRUE, reorder bg bars by weight (desc).
#' @param ... Passed to base plot primitives.
#' @export
plot.time_contrast_clusterpca <- function(x, cluster = 1, type = c("spectrum","time","bg"),
                                          comps = 1:3, reorder_bg = TRUE, ...) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  type <- match.arg(type)
  ci <- if (is.numeric(cluster)) as.integer(cluster) else match(as.character(cluster), x$clusters)
  if (is.na(ci) || ci < 1L || ci > length(x$fits)) stop("Invalid 'cluster'.")

  f <- x$fits[[ci]]
  if (type == "spectrum") {
    vals <- (f$d)^2
    if (!length(vals)) { plot.new(); title("No components"); return(invisible()) }
    ttl <- paste0("Spectrum (", f$bg$mode %||% "whiten", "): ", f$cluster)
    barplot(vals, border = NA, ylab = expression(d[i]^2), xlab = "component", main = ttl, ...)
    abline(h = 0, col = "gray80")
    return(invisible())
  }

  if (type == "time") {
    if (!length(f$d)) { plot.new(); title("No components"); return(invisible()) }
    comps <- comps[comps >= 1 & comps <= ncol(f$u)]
    if (!length(comps)) comps <- 1
    mat <- f$u[, comps, drop = FALSE]
    Tn <- nrow(mat)
    ylim <- range(mat)
    plot(seq_len(Tn), mat[, 1], type = "l", lwd = 2,
         xlab = "time", ylab = "filter value", ylim = ylim,
         main = paste("Time filters:", f$cluster), ...)
    if (ncol(mat) > 1L) {
      for (j in 2:ncol(mat)) lines(seq_len(Tn), mat[, j], lwd = 2)
    }
    abline(h = 0, col = "gray80", lty = 2)
    legend("topright", legend = paste0("comp ", comps), lwd = 2, bty = "n")
    return(invisible())
  }

  if (type == "bg") {
    w <- f$bg$w
    if (is.null(w) || !length(w)) { plot.new(); title("No background"); return(invisible()) }
    # hide self-weight
    w_self <- w
    w_self[ci] <- 0
    ord <- if (reorder_bg) order(w_self, decreasing = TRUE) else seq_along(w_self)
    labs <- x$clusters[ord]
    barplot(w_self[ord], names.arg = labs, las = 2,
            ylab = "weight", main = paste("Background weights:", f$cluster), border = NA, ...)
    return(invisible())
  }
}

#' Summarize background ranks/tau per cluster
#' @inherit background_info
#' @export
summary.time_contrast_clusterpca <- function(object, ...) {
  print(background_info(object))
  invisible(object)
}