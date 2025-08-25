# Debug script for hcluspca test failures

library(subpca)
library(multivarious)

# ==== TEST 1: Function-based ccomp issue ====
cat("=== DEBUG TEST 1: Function-based ccomp ===\n")

set.seed(345)
X <- matrix(rnorm(60*40), 60, 40)
d <- dist(1:60)
hclus <- hclust(d)
cuts <- c(3, 6)

# The problematic function
f <- function(fit, i) {
  max(1, round(log(shape(fit)[1])))
}

# Debug: What is ccomp when it's a function?
print("ccomp type:")
print(typeof(f))
print("ccomp class:")
print(class(f))

# The issue is in line 72 of hcluspca.R:
# message("fitting pca, method ", svd_method, " ncomp = ", ncomp)
# When ncomp is a function, the message function tries to convert it to character

# Let's see what happens when we try to convert function to character in message
try({
  message("ncomp = ", f)
}, silent = FALSE)

cat("\n=== DEBUG TEST 2: Orthogonalization issue ===\n")

set.seed(678)
X2 <- matrix(rnorm(60*40), 60, 40)
d2 <- dist(1:60)
hclus2 <- hclust(d2)
cuts2 <- c(3, 6)

# The issue is in the orthogonalization code at line 120:
# Z <- S[sind[[j]],]
# keep <- which(Z[1,] != 0)

# Let's simulate what might be happening
cat("Simulating orthogonalization issue:\n")
# When S is a single column vector, Z[1,] fails because Z is 1D
S_example <- matrix(1:10, ncol=1)  # Single column
print("S dimensions:")
print(dim(S_example))
print("Trying Z[1,] on single column matrix:")
Z <- S_example[1:3,]  # This becomes a vector when ncol=1
print("Z type after subsetting:")
print(class(Z))
print("Z dimensions:")
print(dim(Z))

# This is the problematic line that would fail:
try({
  keep <- which(Z[1,] != 0)
}, silent = FALSE)

cat("\n=== DEBUG TEST 3: Reconstruction issue ===\n")

set.seed(789)
X3 <- matrix(rnorm(40*30), 40, 30)
d3 <- dist(1:40)
hclus3 <- hclust(d3)
cuts3 <- c(2, 4)

# Run the hcluspca to see what the result looks like
result <- hcluspca(X3, hclus3, cuts3, ccomp = 2)

cat("Result class:")
print(class(result))
cat("Result structure:\n")
str(result, max.level = 1)

# The issue is in reconstruct - let's see what happens
cat("Attempting reconstruct:\n")
try({
  X_recon <- reconstruct(result)
}, silent = FALSE)

# Let's examine the components that might be causing the issue
cat("Components structure:\n")
comps <- components(result)
print(class(comps))
print(dim(comps))
print(str(comps))