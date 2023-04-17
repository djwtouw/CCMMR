#' Find a target number of clusters in the data using convex clustering
#'
#' @description \code{convex_clustering} attempts to find the number of clusters
#' specified by the user by means of convex clustering. The algorithm looks for
#' each number of clusters between \code{target_low} and \code{target_high}. If
#' \code{target_low} = \code{target_high}, the algorithm searches for a single
#' clustering. It is recommended to specify a range around the desired number of
#' clusters, as not each number of clusters between 1 and \code{nrow(X)} may be
#' attainable due to numerical inaccuracies.
#'
#' @param X An \eqn{n} x \eqn{p} numeric matrix. This function assumes that each
#' row represents an object with \eqn{p} attributes.
#' @param W A \code{sparseweights} object, see \link{sparse_weights}.
#' @param target_low Lower bound on the number of clusters that should be
#' searched for. If \code{target_high = NULL}, this is the exact number of
#' clusters that is searched for.
#' @param target_high Upper bound on the number of clusters that should be
#' searched for. Default is \code{NULL}, in that case, it is set equal to
#' \code{target_low}.
#' @param max_iter_phase_1 Maximum number of iterations to find an upper and
#' lower bound for the value for lambda for which the desired number of clusters
#' is attained. Default is 2000.
#' @param max_iter_phase_2 Maximum number of iterations to to refine the upper
#' and lower bounds for lambda. Default is 20.
#' @param lambda_init The first value for lambda other than 0 to use for convex
#' clustering. Default is 0.01.
#' @param factor The percentage by which to increase lambda in each step.
#' Default is 0.025.
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
#' \eqn{p} x \eqn{no. lambdas} values have to be stored, this may require too
#' much memory for large data sets.
#' @param verbose Verbosity of the information printed during clustering.
#' Default is 0, no output.
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
#' \item{\code{phase_1_instances}}{The number of instances of the loss function
#' that were minimized while finding an upper and lower bound for lambda. The
#' sum \code{phase_1_iterations + phase_2_iterations} gives the total number of
#' instances solved.}
#' \item{\code{phase_2_instances}}{The number of instances of the loss function
#' that were minimized while refining the value for lambda. The sum
#' \code{phase_1_iterations + phase_2_iterations} gives the total number of
#' instances solved.}
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
#' # Perform convex clustering with a target number of clusters
#' res1 = convex_clustering(X, W, target_low = 2, target_high = 5)
#'
#' # Plot the clustering for 2 to 5 clusters
#' oldpar = par(mfrow=c(2, 2))
#' plot(X, col = clusters(res1, 2), main = "2 clusters", pch = 19)
#' plot(X, col = clusters(res1, 3), main = "3 clusters", pch = 19)
#' plot(X, col = clusters(res1, 4), main = "4 clusters", pch = 19)
#' plot(X, col = clusters(res1, 5), main = "5 clusters", pch = 19)
#'
#' # A more generalized approach to plotting the results of a range of clusters
#' res2 = convex_clustering(X, W, target_low = 2, target_high = 7)
#'
#' # Plot the clusterings
#' k = length(res2$num_clusters)
#' par(mfrow=c(ceiling(k / ceiling(sqrt(k))), ceiling(sqrt(k))))
#'
#' for (i in 1:k) {
#'     labels = clusters(res2, res2$num_clusters[i])
#'     c = length(unique(labels))
#'
#'     plot(X, col = labels, main = paste(c, "clusters"), pch = 19)
#' }
#' par(oldpar)
#'
#' @seealso \link{convex_clusterpath}, \link{sparse_weights}
#'
#' @export
convex_clustering <- function(X, W, target_low, target_high = NULL,
                              max_iter_phase_1 = 2000, max_iter_phase_2 = 20,
                              lambda_init = 0.01, factor = 0.025, tau = 1e-3,
                              center = TRUE, scale = TRUE, eps_conv = 1e-6,
                              burnin_iter = 25, max_iter_conv = 5000,
                              save_clusterpath = FALSE, verbose = 0)
{
    # Input checks
    .check_array(X, 2, "X")
    .check_weights(W)
    .check_int(max_iter_phase_1, FALSE, "max_iter_phase_1")
    .check_int(max_iter_phase_2, FALSE, "max_iter_phase_2")
    .check_scalar(lambda_init, TRUE, "lambda_init")
    .check_scalar(factor, TRUE, "factor")
    .check_scalar(tau, TRUE, "tau", upper_bound = 1)
    .check_boolean(center, "center")
    .check_boolean(scale, "scale")
    .check_scalar(eps_conv, TRUE, "eps_conv", upper_bound = 1)
    .check_int(burnin_iter, FALSE, "burnin_iter")
    .check_int(max_iter_conv, FALSE, "max_iter_conv")
    .check_boolean(save_clusterpath, "save_clusterpath")
    .check_int(verbose, FALSE, "verbose")

    if (is.null(target_high)) {
        target_high = target_low
    }
    .check_cluster_targets(target_low, target_high, nrow(X))

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
    clust = .convex_clustering(X_, W_idx, W_val, eps_conv, eps_fusions, scale,
                               save_clusterpath, burnin_iter, max_iter_conv,
                               target_low, target_high, max_iter_phase_1,
                               max_iter_phase_2, verbose, lambda_init, factor)
    elapsed_time = difftime(Sys.time(), t_start, units = "secs")

    # Stop the program if no clusterings in the range [low, high] have been
    # found
    if (clust$targets_found < 1) {
        message = paste("No clusterings within the range [target_low,",
                        "target_high] were found. Consider a wider range or",
                        "set connected = TRUE in sparse_weights(...)")
        stop(message)
    }

    # Construct result
    result = list()
    result$info = data.frame(
        clust$info_d[1, ],
        clust$info_i[2, ],
        clust$info_d[2, ]
    )
    result$info = result$info[1:clust$targets_found, ]
    names(result$info) = c("lambda", "clusters", "loss")

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

    # Add the number of instances solved
    result$phase_1_instances = clust$phase_1_instances
    result$phase_2_instances = clust$phase_2_instances

    # Add vector of possible cluster counts
    result$num_clusters = result$info$clusters

    # Add the number of observations
    result$n = nrow(X)

    # Give the result a class
    class(result) = "cvxclust"

    return(result)
}
