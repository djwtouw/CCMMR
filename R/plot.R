#' Plot 2D clusterpath
#'
#' @description Plot a clusterpath for two-dimensional data.
#'
#' @param x A \code{cvxclust} object.
#' @param col A vector containing cluster membership information. Default is
#' \code{NULL}.
#' @param labels A vector containing labels for each object. Default is
#' \code{NULL}.
#' @param ... Further graphical parameters.
#'
#' @examples
#' # Load data
#' data(two_half_moons)
#' data = as.matrix(two_half_moons)
#' X = data[, -3]
#' y = data[, 3]
#'
#' # Get sparse distances in dictionary of keys format with k = 5 and phi = 8
#' W = sparse_weights(X, 5, 8.0)
#'
#' # Set a sequence for lambda
#' lambdas = seq(0, 2400, 1)
#'
#' # Compute results CMM
#' res = convex_clusterpath(X, W, lambdas)
#' plot(res, y + 1)
#'
#' @export
plot.cvxclust <- function(x, col = NULL, labels = NULL, ...)
{
    # Input checks
    .check_cvxclust(x, "x")

    if (is.null(x$coordinates)) {
        message = paste("The clusterpath coordinates were not saved. Make sure",
                        "to set save_clusterpath = TRUE in",
                        "convex_clustering(...) or convex_clusterpath(...)")
        stop(message)
    }

    if (ncol(x$coordinates) != 2) {
        stop(paste("plot.clusterpath is only implemented for two-dimensional",
                   "clusterpaths"))
    }

    # Get the number of observations and the number of lambdas used for the path
    n_obs = x$n
    n_lam = length(x$lambdas)

    # X is the first set of observations
    X = x$coordinates[c(1:n_obs), ]
    plot(X, pch = 1, las = 1, xlab = expression(X[1]), ylab = expression(X[2]),
         asp = 1, ...)

    # Draw the paths
    for (i in 1:n_obs) {
        graphics::lines(x$coordinates[seq(i, n_obs * n_lam, n_obs), ],
                        col = "grey", cex=0.5)
    }

    # Give the dots colors to denote cluster membership
    if (!is.null(col)) {
        graphics::points(X[, 1:2], pch = 19, col = col)
    } else {
        graphics::points(X[, 1:2], pch = 19, col = "black")
    }

    # Add labels
    if (!is.null(labels)) {
        graphics::text(X[, 1], X[, 2], labels = labels, cex= 0.7, pos=3)
    }
}
