#' Obtain clustering from a clusterpath
#'
#' @description Get a particular clustering of the data. If there is a
#' clustering for \code{n_clusters}, it is returned, otherwise the function will
#' stop with a message. To know whether a query is going to be successful
#' beforehand, check the \code{num_clusters} attribute of the \code{cvxclust}
#' object, this lists all possible options for the number of clusters.
#'
#' @param obj A \code{cvxclust} object.
#' @param n_clusters An integer that specifies the number of clusters that
#' should be returned.
#'
#' @return A vector with the cluster labels for each object in the data.
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
#'
#' # Get labels for three clusters
#' labels = clusters(res, 3)
#'
#' @export
clusters <- function(obj, n_clusters)
{
    # Input checks
    .check_cvxclust(obj)
    .check_int(n_clusters, TRUE, "n_clusters")

    if (!(n_clusters %in% obj$num_clusters)) {
        message = paste(n_clusters, "is not among the possible number of",
                        "clusters")
        stop(message)
    }

    # Start with an entry in a hashmap for each observation
    D = list()
    for (i in 1:obj$n) {
        D[as.character(-i)] = i
    }

    # Work through the merge table to reduce the number of clusters until
    # the desired number is reached
    for (i in 1:nrow(obj$merge)) {
        if (length(D) == n_clusters) {
            break
        }

        D[[as.character(i)]] = c(D[[as.character(obj$merge[i, 1])]],
                                 D[[as.character(obj$merge[i, 2])]])
        D[as.character(obj$merge[i, 1])] = NULL
        D[as.character(obj$merge[i, 2])] = NULL
    }

    # Create cluster labels
    result = rep(0, obj$n)
    label = 1

    for (key in names(D)) {
        result[D[[key]]] = label

        # Increment label
        label = label + 1
    }

    return(result)
}
