#' Computation of sparse weight matrix
#'
#' @description Construct a sparse weight matrix in a dictionary-of-keys format.
#' Each nonzero weight is computed as \eqn{exp(-phi * ||x_i - x_j||^2)}, where
#' the squared Euclidean distance may be scaled by the average squared Euclidean
#' distance, depending on the argument \code{scale}.
#'
#' @param X An \eqn{n} x \eqn{p} numeric matrix. This function assumes that each
#' row represents an object with \eqn{p} attributes.
#' @param k The number of nearest neighbors to be used for non-zero weights.
#' @param phi Tuning parameter of the Gaussian weights. Input should be a
#' nonnegative value.
#' @param connected If \code{TRUE}, guarantee a connected structure of the
#' weight matrix by using a symmetric circulant matrix to add nonzero weights.
#' This ensures that groups of observations that would not be connected through
#' weights that are based only on the k nearest neighbors are (indirectly)
#' connected anyway. Default is \code{TRUE}.
#' @param scale If \code{TRUE}, scale each squared l2-norm by the mean squared
#' l2-norm to ensure scale invariance of the weights. Default is \code{TRUE}.
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
sparse_weights <- function(X, k, phi, connected = TRUE, scale = TRUE)
{
    # Input checks
    .check_array(X, 2, "X")
    .check_int(k, TRUE, "k")
    .check_scalar(phi, FALSE, "phi")
    .check_boolean(connected, "connected")
    .check_boolean(scale, "scale")

    # Preliminaries
    n = nrow(X)

    # Get the k nearest neighbors
    nn_res = nn2(X, X, k + 1)
    nn_idx = nn_res$nn.idx - 1
    nn_dists = nn_res$nn.dists

    # Transform the indices of the k-nn into a dictionary of keys sparse matrix
    res = .sparse_weights(t(X), t(nn_idx), t(nn_dists), phi, k, connected,
                          scale)

    # Unique keys and value pairs
    keys = t(res$keys)
    values = res$values
    u_idx = !duplicated(keys)

    keys = keys[u_idx, ]
    values = values[u_idx]

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
