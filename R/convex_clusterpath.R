#' Minimize the convex clustering loss function
#'
#' @description Minimizes the convex clustering loss function for a given set of
#' values for lambda.
#'
#' @param X An \eqn{n} x \eqn{p} numeric matrix. This function assumes that each
#' row represents an object with \eqn{p} attributes.
#' @param W A \code{sparseweights} object, see \link{sparse_weights}.
#' @param lambdas A vector containing the values for the penalty parameter.
#' @param tau Parameter to compute the threshold to fuse clusters. Default is
#' 0.001.
#' @param center If \code{TRUE}, center \code{X} so that each column has mean
#' zero. Default is \code{TRUE}.
#' @param scale If \code{TRUE}, scale the loss function to ensure that the
#' cluster solution is invariant to the scale of \code{X}. Default is
#' \code{TRUE}. Not recommended to set to \code{FALSE} unless comparing to
#' algorithms that minimize the unscaled convex clustering loss function.
#' @param eps_conv Parameter for determining convergence of the minimization.
#' Default is 1e-6.
#' @param burnin_iter Number of updates of the loss function that are done
#' without step doubling. Default is 25.
#' @param max_iter_conv Maximum number of iterations for minimizing the loss
#' function. Default is 5000.
#' @param save_clusterpath If \code{TRUE}, store the solution that minimized
#' the loss function for each lambda. Is required for drawing the clusterpath.
#' Default is \code{FALSE}. To store the clusterpath coordinates, \eqn{n} x
#' \eqn{p} x \eqn{no. lambdas} have to be stored, this may require too much
#' memory for large data sets.
#' @param target_losses The values of the loss function that are used to
#' determine convergence of the algorithm (tested as: loss - target <=
#' \code{eps_conv} * target). If the input is not \code{NULL}, it should be a
#' vector with the same length as \code{lambdas}. Great care should be exercised
#' to make sure that the target losses correspond to attainable values for the
#' minimization. The inputs (\code{X}, \code{W}, \code{lambdas}) should be the
#' same, but also the same version of the loss function (centered, scaled)
#' should be used. Default is \code{NULL}.
#'
#' @return A \code{cvxclust} object containing the following
#' \item{\code{info}}{A dataframe containing for each value for lambda: the
#' number of different clusters, and the value of the loss function at the
#' minimum.}
#' \item{\code{merge}}{The merge table containing the order at which the
#' observations in \code{X} are clustered.}
#' \item{\code{height}}{The value for lambda at which each reduction in the
#' number of clusters occurs.}
#' \item{\code{order}}{The order of the observations in \code{X} in order to
#' draw a dendrogram without conflicting branches.}
#' \item{\code{elapsed_time}}{The number of seconds that elapsed while
#' running the code. Note that this does not include the time required for
#' input checking and possibly scaling and centering \code{X}.}
#' \item{\code{coordinates}}{The clusterpath coordinates. Only part of the
#' output in case that \code{save_clusterpath=TRUE}.}
#' \item{\code{lambdas}}{The values for lambda for which a clustering was
#' found.}
#' \item{\code{eps_fusions}}{The threshold for cluster fusions that was used by
#' the algorithm.}
#' \item{\code{num_clusters}}{The different numbers of clusters that have been
#' found.}
#' \item{\code{n}}{The number of observations in \code{X}.}
#'
#' @examples
#' # Load data
#' data(two_half_moons)
#' data = as.matrix(two_half_moons)
#' X = data[, -3]
#' y = data[, 3]
#'
#' # Get sparse weights in dictionary of keys format with k = 5 and phi = 8
#' W = sparse_weights(X, 5, 8.0)
#'
#' # Set a sequence for lambda
#' lambdas = seq(0, 2400, 1)
#'
#' # Compute clusterpath
#' res = convex_clusterpath(X, W, lambdas)
#'
#' # Get cluster labels for two clusters
#' labels = clusters(res, 2)
#'
#' # Plot the clusterpath with colors based on the cluster labels
#' plot(res, col = labels)
#'
#' @seealso \link{convex_clustering}, \link{sparse_weights}
#'
#' @export
convex_clusterpath <- function(X, W, lambdas, tau = 1e-3, center = TRUE,
                               scale = TRUE, eps_conv = 1e-6, burnin_iter = 25,
                               max_iter_conv = 5000, save_clusterpath = TRUE,
                               target_losses = NULL)
{
    # Input checks
    .check_array(X, 2, "X")
    .check_weights(W)
    .check_lambdas(lambdas)
    .check_scalar(tau, TRUE, "tau", upper_bound = 1)
    .check_boolean(center, "center")
    .check_boolean(scale, "scale")
    .check_scalar(eps_conv, TRUE, "eps_conv", upper_bound = 1)
    .check_int(burnin_iter, FALSE, "burnin_iter")
    .check_int(max_iter_conv, FALSE, "max_iter_conv")
    .check_boolean(save_clusterpath, "save_clusterpath")

    # Check the vector of target losses
    if (!is.null(target_losses)) {
        .check_array(target_losses, 1, "target_losses")

        if (length(lambdas) != length(target_losses)) {
            message = "target_losses should have the same length as lambdas"
            stop(message)
        }

        use_target = TRUE
    } else {
        target_losses = rep(-1, length(lambdas))
        use_target = FALSE
    }

    # Preliminaries
    n = nrow(X)

    # Set the means of each column of X to zero
    if (center) {
        X_ = X - matrix(apply(X, 2, mean), byrow = TRUE, ncol = ncol(X),
                        nrow = nrow(X))
    } else {
        X_ = X
    }

    # Transpose X
    X_ = t(X_)

    # Separate the weights into keys and values
    W_idx = t(W$keys) - 1
    W_val = W$values

    # Compute fusion threshold
    eps_fusions = .fusion_threshold(X_, tau)

    t_start = Sys.time()
    clust = .convex_clusterpath(X_, W_idx, W_val, lambdas, target_losses,
                                eps_conv, eps_fusions, scale, save_clusterpath,
                                use_target, burnin_iter, max_iter_conv)
    elapsed_time = difftime(Sys.time(), t_start, units = "secs")

    # Construct result
    result = list()
    result$info = data.frame(
        clust$info_d[1, ],
        clust$info_i[2, ],
        clust$info_d[2, ],
        clust$info_i[1, ]
    )
    names(result$info) = c("lambda", "clusters", "loss", "iterations")

    # Merge table and height vector
    result$merge = t(clust$merge)
    result$height = clust$height

    # Determine order of the observations for a dendrogram
    # Start with an entry in a hashmap for each observation
    d = hashmap()
    for (i in 1:n) {
        d[[-i]] = i
    }

    # Work through the merge table to make sure that everything that is
    # merged is next to each other
    if (nrow(result$merge) >= 1) {
        for (i in 1:nrow(result$merge)) {
            d[[i]] = c(d[[result$merge[i, 1]]], d[[result$merge[i, 2]]])
            delete(d, result$merge[i, 1])
            delete(d, result$merge[i, 2])
        }
    }

    # Finally, create a vector with the order of the observations
    result$order = c()
    keys = keys(d)
    for (key in keys) {
        result$order = c(result$order, d[[key]])
    }

    # Add elapsed time
    result$elapsed_time = elapsed_time

    # Add clusterpath coordinates
    if (save_clusterpath) {
        result$coordinates = t(clust$clusterpath)
    }

    # Add lambdas
    result$lambdas = result$info$lambda

    # Add fusion threshold
    result$eps_fusions = eps_fusions

    # Add vector of possible cluster counts
    result$num_clusters = unique(result$info$clusters)

    # Add the number of observations
    result$n = nrow(X)

    # Give the result a class
    class(result) = "cvxclust"

    return(result)
}
