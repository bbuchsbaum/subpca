# library(Matrix)
# library(CVXR)
#
# X <- matrix(rnorm(100*100), 100,100)
# B <- matrix(rnorm(100*10), 100,10)
# g <- neighborweights::graph_weights(X, k=5, neighbor_mode="knn", weight_mode="normalized", sigma=.7)
# A <- neighborweights::adjacency(g)
# D <- rowSums(A)
# L <- D - A
#
# lap_reg <- function(beta, lambda, alpha) {
#   pen1 <- alpha * CVXR::matrix_trace(quad_form(beta, L))
#   pen2 <- lambda * CVXR::norm2(beta)
#   pen1 + pen2
#   #pen2
# }
# alpha=.1
# lambda=.1
# n <- nrow(X)
# beta <- Variable(ncol(B), ncol(X))
#
# #loss <- CVXR::norm2(X - B %*% beta)
# loss <- sum((X - B %*% beta)^2) / (2 * n)
# obj <- loss + lap_reg(beta,alpha, lambda)
# prob <- Problem(Minimize(obj))
# result <- solve(prob)
# result$getValue(beta)
#
