#' Conversion of a \code{cvxclust} object into an \code{hclust} object
#'
#' @description Converts the output of \link{convex_clustering} or
#' \link{convex_clusterpath} into a hclust object. Note that a step in the
#' clusterpath from one value for lambda to the next may cause the number of
#' clusters to decrease by more than one. It is a hard requirement that the
#' clusterpath ends in a single cluster, as standard dendrogram plotting
#' methods fail if this is not the case.
#'
#' @param x A \code{cvxclust} object.
#' @param ... Unused.
#'
#' @return A \code{hclust} object.
#'
#' @examples
#' # Demonstration of converting a clusterpath into a dendrogram, first generate
#' # data
#' set.seed(6)
#' X = matrix(rnorm(14), ncol = 2)
#' y = rep(1, nrow(X))
#'
#' # Get sparse distances in dictionary of keys format with k = 3
#' W = sparse_weights(X, 3, 4.0)
#'
#' # Sequence for lambda
#' lambdas = seq(0, 45, 0.02)
#'
#' # Compute results
#' res = convex_clusterpath(X, W, lambdas)
#'
#' # Generate hclust object
#' hcl = as.hclust(res)
#' hcl$height = sqrt(hcl$height)
#'
#' # Plot clusterpath and dendrogram
#' par(mfrow=c(1, 2))
#' plot(res, y, label = c(1:7))
#' plot(hcl, ylab = expression(sqrt(lambda)), xlab = NA, sub = NA, main = NA,
#'      hang = -1)
#'
#' @seealso \link{hclust}
#'
#' @export
as.hclust.cvxclust <- function(x, ...)
{
    # Input checks
    .check_cvxclust(x, "obj")

    if (!(1 %in% x$num_clusters)) {
        message = paste("The clusterpath does not terminate in a single",
                        "cluster. Consider setting connected = TRUE in",
                        "sparse_weights(...) or choose a higher maximum value",
                        "for lambda")
        stop(message)
    }

    # Create hclust object
    result = list()
    result$merge = x$merge
    result$height = x$height
    result$order = x$order
    result$labels = c(1:x$n)

    result$call = list()
    result$call[[1]] = "Convex Clustering"
    result$call$d[[1]] = "dist"
    result$call$d[[2]] = "X"

    result$dist.method = "Euclidean"

    # Set class
    class(result) = "hclust"

    return(result)
}
