# CCMMR
CCMMR implements convex clustering using the minimization algorithm presented in the paper _Convex Clustering through MM: An Efficient Algorithm to Perform Hierarchical Clustering_ by D.J.W. Touw, P.J.F. Groenen, and Y. Terada. For issues, please use [Github Issues](https://github.com/djwtouw/CCMMR/issues).

There is also a [Python package](https://github.com/djwtouw/CCMMPy) available.

## Contents
[Installation](#installation)

[Examples](#examples)

## Installation
CCMMR has the following dependencies:
- r2r
- RANN
- Rcpp
- RcppEigen

To install CCMMR, clone the repository, open `CCMMR.Rproj` in RStudio, and press install in the build panel.

## Examples
The following examples are also part of the documentation of the package. Start with loading the package.
```R
library(CCMMR)
```
### Example 1: Computation of a clusterpath
After loading the data, a sparse weight matrix is constructed based on the `k = 5` nearest neighbors. This means that nonzero weights are computed only for pairs of objects that are _k_ nearest neighbors of each other. By default, the weight matrix is constructed so that every observation is (in)directly connected to all other observations via nonzero weights. This ensures that the minimum number of clusters is one. To turn off this behavior, set `connected = FALSE`. 
```R
# Load data
data(two_half_moons)
data = as.matrix(two_half_moons)
X = data[, -3]
y = data[, 3]

# Get sparse weights in dictionary of keys format with k = 5 and phi = 8
W = sparse_weights(X, 5, 8.0)

# Set a sequence for lambda
lambdas = seq(0, 2400, 1)

# Compute clusterpath
res = convex_clusterpath(X, W, lambdas)

# Get cluster labels for two clusters
labels = clusters(res, 2)

# Plot the clusterpath with colors based on the cluster labels
plot(res, col = labels)
```
### Example 2: Searching for a number of clusters
In the previous example, the choice for $\lambda$ has determined what the number of clusters was going to be. However, it can be difficult to guess in advance what value for $\lambda$ corresponds to a particular number of clusters. The following code looks for clusterings in a specified range. If no upper bound is specified, just a single number of clusters (equal to `target_low`) is looked for.
```Python
# Load data
data(two_half_moons)
data = as.matrix(two_half_moons)
X = data[, -3]
y = data[, 3]

# Get sparse weights in dictionary of keys format with k = 5 and phi = 8
W = sparse_weights(X, 5, 8.0)

# Perform convex clustering with a target number of clusters
res1 = convex_clustering(X, W, target_low = 2, target_high = 5)

# Plot the clustering for 2 to 5 clusters
par(mfrow=c(2, 2))
plot(X, col = clusters(res1, 2), main = "2 clusters", pch = 19)
plot(X, col = clusters(res1, 3), main = "3 clusters", pch = 19)
plot(X, col = clusters(res1, 4), main = "4 clusters", pch = 19)
plot(X, col = clusters(res1, 5), main = "5 clusters", pch = 19)

# A more generalized approach to plotting the results of a range of clusters
res2 = convex_clustering(X, W, target_low = 2, target_high = 7)

# Plot the clusterings
k = length(res2$num_clusters)
par(mfrow=c(ceiling(k / ceiling(sqrt(k))), ceiling(sqrt(k))))

for (i in 1:k) {
    labels = clusters(res2, res2$num_clusters[i])
    c = length(unique(labels))

    plot(X, col = labels, main = paste(c, "clusters"), pch = 19)
}
```
### Example 3: Alternative visualizations
As an alternative to the clusterpath, convex clustering results can also be visualized using a dendrogram. In the following example, convex clustering is applied to a small generated data set, after which the `as.hclust()` function transforms the output into a `hclust` object. Consequently, the standard `plot()` can be used to plot a dendrogram. Note that `hclust` objects require the clusterpath to terminate in a single cluster.
```R
# Demonstration of converting a clusterpath into a dendrogram, first generate
# data
set.seed(6)
X = matrix(rnorm(14), ncol = 2)
y = rep(1, nrow(X))

# Get sparse distances in dictionary of keys format with k = 3
W = sparse_weights(X, 3, 4.0)

# Sequence for lambda
lambdas = seq(0, 45, 0.02)

# Compute results
res = convex_clusterpath(X, W, lambdas)

# Generate hclust object
hcl = as.hclust(res)
hcl$height = sqrt(hcl$height)

# Plot clusterpath and dendrogram
par(mfrow=c(1, 2))
plot(res, y, label = c(1:7))
plot(hcl, ylab = expression(sqrt(lambda)), xlab = NA, sub = NA, main = NA,
     hang = -1)
```
