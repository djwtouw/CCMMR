#' Computation of sparse weight matrix
#'
#' @description Construct a sparse weight matrix in a dictionary-of-keys format.
#' Each nonzero weight is computed as \eqn{exp(-phi * ||x_i - x_j||^2)}, where
#' the squared Euclidean distance may be scaled by the average squared Euclidean
#' distance, depending on the argument \code{scale}. Sparsity is achieved by
#' only setting weights to nonzero values that correspond to two objects that
#' are among each other's \eqn{k} nearest neighbors.
#'
#' @param X An \eqn{n} x \eqn{p} numeric matrix. This function assumes that each
#' row represents an object with \eqn{p} attributes.
#' @param k The number of nearest neighbors to be used for non-zero weights.
#' @param phi Tuning parameter of the Gaussian weights. Input should be a
#' nonnegative value.
#' @param connected If \code{TRUE}, guarantee a connected structure of the
#' weight matrix. This ensures that groups of observations that would not be
#' connected through weights that are based only on the \code{k} nearest
#' neighbors are (indirectly) connected anyway. The method is determined by
#' the argument \code{connection_type}. Default is \code{TRUE}.
#' @param scale If \code{TRUE}, scale each squared l2-norm by the mean squared
#' l2-norm to ensure scale invariance of the weights. Default is \code{TRUE}.
#' @param connection_type Determines the method to ensure a connected weight
#' matrix if \code{connected} is \code{TRUE}. Should be one of
#' \code{c("SC", "MST")}. SC stands for the method using a symmetric circulant
#' matrix, connecting objects \eqn{i} with objects \eqn{i+1} (and \eqn{n} with
#' \eqn{1}). MST stands for minimum spanning tree. The graph that results from
#' the nonzero weights determined by the \eqn{k} nearest neighbors is divided
#' into \eqn{c} subgraphs and a minimum spanning tree algorithm is used to add
#' \eqn{c-1} nonzero weights to ensure that all objects are indirectly
#' connected. Default is \code{"SC"}.
#'
#' @return A \code{sparseweights} object containing the nonzero weights in
#' dictionary-of-keys format.
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
#' @export
sparse_weights <- function(X, k, phi, connected = TRUE, scale = TRUE,
                           connection_type = "SC")
{
    # Input checks
    CCMMR:::.check_array(X, 2, "X")
    CCMMR:::.check_int(k, TRUE, "k")
    CCMMR:::.check_scalar(phi, FALSE, "phi")
    CCMMR:::.check_boolean(connected, "connected")
    CCMMR:::.check_boolean(scale, "scale")
    if (!(connection_type %in% c("SC", "MST"))) {
        message = "Expected one of 'SC' and 'MST' for connection_type"
        stop(message)
    }

    # Preliminaries
    n = nrow(X)

    # Get the k nearest neighbors
    nn_res = RANN::nn2(X, X, k + 1)
    nn_idx = nn_res$nn.idx - 1
    nn_dists = nn_res$nn.dists

    # Transform the indices of the k-nn into a dictionary of keys sparse matrix
    res = CCMMR:::.sparse_weights(t(X), t(nn_idx), t(nn_dists), phi, k,
                                  (connection_type == "SC") && connected, scale)

    # Unique keys and value pairs
    keys = t(res$keys)
    values = res$values
    u_idx = !duplicated(keys)

    keys = keys[u_idx, ]
    values = values[u_idx]

    if (connection_type == "MST" && connected) {
        # Use the keys of the sparse weight matrix to find clusters in the data
        id = CCMMR:::.find_subgraphs(t(keys), nrow(X)) + 1

        # Number of disconnected parts of the graph
        n_clusters = max(id)

        if (n_clusters > 1) {
            # List to keep track of which objects are responsible for those
            # distances
            closest_objects = list()

            # Matrix to keep track of the shortest distance between the clusters
            D_between = matrix(0, nrow = n_clusters, ncol = n_clusters)

            for (a in 2:nrow(D_between)) {
                for (b in 1:(a - 1)) {
                    # Find out which member of cluster a is closest to each of
                    # the members of cluster b
                    nn_between = RANN::nn2(X[id == a, ], X[id == b, ], 1)

                    # Get the indices of the objects with respect to their
                    # cluster
                    b_idx = which.min(nn_between$nn.dists)
                    a_idx = nn_between$nn.idx[b_idx]

                    # Save the distance
                    D_between[a, b] = nn_between$nn.dists[b_idx]
                    D_between[b, a] = nn_between$nn.dists[b_idx]

                    # Get the original indices of the objects
                    a_idx = c(1:nrow(X))[id == a][a_idx] - 1
                    b_idx = c(1:nrow(X))[id == b][b_idx] - 1

                    # Store the objects
                    idx = (a - 1) * (a - 2) / 2 + b
                    closest_objects[[idx]] = c(a_idx, b_idx)
                }
            }

            # Find minimum spanning tree for D_between
            mst_keys = CCMMR:::.find_mst(D_between) + 1

            # Vector for the weights
            mst_values = rep(0, nrow(mst_keys))

            # Get the true object indices from the closest_objects list
            for (i in 1:nrow(mst_keys)) {
                # Get the index for the closest_objects list
                ii = min(mst_keys[i, 1], mst_keys[i, 2]) - 1
                jj = max(mst_keys[i, 1], mst_keys[i, 2]) - 1
                idx = jj * (jj - 1) / 2 + ii + 1

                # Fill in the weights in the values vector
                mst_values[i] = D_between[mst_keys[i, 1], mst_keys[i, 2]]

                # Replace the cluster ids by the object ids
                mst_keys[i, ] = closest_objects[[idx]]
            }

            # Because both the upper and lower part of the weight matrix are
            # stored, the key and value pairs are duplicated
            mst_keys = rbind(mst_keys, mst_keys[, c(2, 1)])
            mst_values = c(mst_values, mst_values)

            # Compute the weights from the distances
            if (scale) {
                mst_values = exp(-phi * mst_values^2 / res$msd)
            } else {
                mst_values = exp(-phi * mst_values^2)
            }

            # Append everything to the existing keys and values, ensuring a
            # connected weight matrix based on a minimum spanning tree
            keys = rbind(keys, mst_keys)
            values = c(values, mst_values)
        }
    }

    # Store the keys in column major format
    sorted_idx = order(keys[, 2], keys[, 1])
    keys = keys[sorted_idx, ]
    values = values[sorted_idx]

    # In R, counting starts at 1
    keys = keys + 1

    # Prepare result
    result = list()
    result$keys = keys
    result$values = values
    class(result) = "sparseweights"

    return(result)
}
