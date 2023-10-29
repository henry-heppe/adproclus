# Visualization functions for an ADPROCLUS model

#' Network plot of a (low dimensional) ADPROCLUS solution
#'
#' Produce a representation of a (low dimensional) ADPROCLUS solution,
#' where each cluster is a vertex and the edge between two vertices represents
#' the overlap between the corresponding clusters.
#' The size of a vertex corresponds to the cluster size.
#' The overlap is represented through color, width and numerical label
#' of the edge.
#' The numerical edge labels can be relative
#' (number of overlap observations / total observations)
#' or absolute (number of observations in both clusters).
#' \strong{NOTE:} This function can be called through the
#' \code{plot(model, type = "Network")} function with model an
#' object of class \code{adpc}.
#'
#' @param model ADPROCLUS solution (class: \code{adpc}). Low dimensional model
#' possible.
#' @param relative_overlap Logical. If \code{TRUE} (default), the number of
#' observations belonging to two clusters
#' is divided by the total number of observations. If \code{FALSE}
#' the number of observations in a cluster overlap will be displayed on the
#' edges.
#' @param title String. Default: " Cluster network of ADPROCLUS solution"
#' @param filetype Optional. Choose type of file to save the plot.
#' Possible choices: \code{"R", "pdf", "svg", "tex", "jpg", "tiff", "png", ""}
#' Default: \code{NULL} does not create a file.
#'
#' @return Invisibly returns the input model.
#' @export
#'
#' @examples
#' # Loading a test dataset into the global environment
#' x <- stackloss
#'
#' # Quick low dimensional clustering with K = 3 clusters and S = 1 dimensions
#' clust <- adproclus_low_dim(x, 3, 1)
#'
#' # Plot the overlapping the clusters
#' plot_cluster_network(clust)
plot_cluster_network <- function(model,
                                 title = "Cluster network of ADPROCLUS
                                 solution",
                                 relative_overlap = TRUE,
                                 filetype = NULL) {
  checkmate::assertClass(model, "adpc")
  checkmate::assertString(title, null.ok = TRUE)
  checkmate::assertFlag(relative_overlap)
  checkmate::assertChoice(filetype, c(
    "R", "pdf", "svg", "tex",
    "jpg", "tiff", "png", ""
  ),
  null.ok = TRUE
  )

  if (is.null(title)) {
    title <- "Cluster network of ADPROCLUS solution"
  }


  withr::local_seed(1)

  A <- model$A
  k <- ncol(A)

  sizes <- colSums(A)

  # Create adjacency matrix for graph, where each cluster is a node
  # All nodes are connected to all nodes except themselves
  adjacency_matrix <- matrix(rep(1, k^2), k, k)
  diag(adjacency_matrix) <- rep(0, k)

  network <- igraph::graph_from_adjacency_matrix(adjacency_matrix,
    mode = "undirected"
  )
  edgelist <- data.frame(t(igraph::as_edgelist(network)))
  weights <- sapply(edgelist, .extract_overlap, A = A)

  # Compute cluster sizes and then add to string of node label
  labels <- c()
  for (i in 1:k) {
    labels[i] <- paste(colnames(A)[i], "\n ", "obs: ", colSums(A)[i], sep = "")
  }

  weights_internal <- weights
  # Fr layout algorithm cannot deal with zero-weights
  weights_internal[which(weights == 0)] <- 0.9

  # Display edges, when there is no overlap at all
  if (sum(weights) == 0) {
    weights_internal <- weights_internal + 1
  }

  if (relative_overlap) {
    qgraph::qgraph(
      input = cbind(igraph::as_edgelist(network), weights_internal),
      mar = c(7, 7, 7, 7),
      layout = "spring",
      minimum = 0.9,
      theme = "TeamFortress",
      colFactor = 0.6,
      filetype = filetype,
      node.width = 1 + 1.5 * (sizes / sqrt(sum(sizes^2))),
      title = title,
      labels = labels,
      label.scale = TRUE,
      edge.labels = round(weights / nrow(A), digits = 4),
      edge.label.cex = 1.5,
      edge.label.color = "black",
      edge.label.bg = FALSE,
      edge.label.position = 0.3,
      edge.width = 1.5 + (weights / sqrt(sum(weights^2) + 0.001))^0.01,
      directed = FALSE,
      nNodes = k,
      weighted = TRUE,
      edgelist = TRUE
    )
  } else {
    qgraph::qgraph(
      input = cbind(igraph::as_edgelist(network), weights_internal),
      mar = c(7, 7, 7, 7),
      layout = "spring",
      minimum = 0.9,
      theme = "TeamFortress",
      colFactor = 0.6,
      filetype = filetype,
      node.width = 1 + 1.5 * (sizes / sqrt(sum(sizes^2))),
      title = title,
      labels = labels,
      label.scale = TRUE,
      edge.labels = weights,
      edge.label.cex = 1.5,
      edge.label.color = "black",
      edge.label.bg = FALSE,
      edge.label.position = 0.3,
      edge.width = 1.5 + (weights / sqrt(sum(weights^2) + 0.001))^0.01,
      directed = FALSE,
      nNodes = k,
      weighted = TRUE,
      edgelist = TRUE
    )
  }
  invisible(model)
}

#' Plot profile matrix of ADPROCLUS solution
#'
#' Produce a representation of profile matrix \eqn{\boldsymbol{P}} of a
#' (low dimensional) ADPROCLUS solution of class \code{adpc}.
#' The plot displays the profiles in the style of a correlation plot.
#' \strong{NOTE:} This function can also be called through the
#' \code{plot(model, type = "Profiles")} function with model an object of
#' class \code{adpc}.
#'
#' @param model Object of class \code{adpc}. (Low dimensional) ADPROCLUS
#' solution
#' @param title String. Default: "Profiles of ADPROCLUS solution"
#'
#' @return Invisibly returns the input model.
#' @export
#'
#' @examples
#' # Loading a test dataset into the global environment
#' x <- stackloss
#'
#' # Quick clustering with K = 3 clusters
#' clust <- adproclus(x, 3)
#'
#' # Plot the profile scores of each cluster
#' plot_profiles(clust)
plot_profiles <- function(model, title = "Profiles of ADPROCLUS solution") {
  checkmate::assertClass(model, "adpc")
  checkmate::assertString(title, null.ok = TRUE)

  if (is.null(title)) {
    title <- "Profiles of ADPROCLUS solution"
  }
  if (is.null(model$C)) {
    corrplot::corrplot(model$P,
      is.corr = FALSE, title = title,
      mar = c(0, 0, 2, 0)
    )
  } else {
    if (title == "Profiles of ADPROCLUS solution") {
      title <- "Profiles of Low dimensional ADPROCLUS solution"
    }
    corrplot::corrplot(model$C,
      is.corr = FALSE, title = title,
      mar = c(0, 0, 2, 0)
    )
  }
  invisible(model)
}

#' Plot variable to component matrix of ADPROCLUS solution
#'
#' Produce a representation of variable to component matrix
#' \eqn{\boldsymbol{B'}} of a \strong{low dimensional} ADPROCLUS solution
#' of class \code{adpc}. The plot displays the scores in the style of a
#' correlation plot.
#' \strong{NOTE:} This function can be called through the
#' \code{plot(model, type = "VarsByComp")} function
#' with model an object of class \code{adpc}.
#'
#' @param model Object of class \code{adpc}. Must be \strong{Low dimensional}
#' ADPROCLUS solution
#' @param title String. Default: "B' of Low Dimensional ADPROCLUS Solution"
#'
#' @return Invisibly returns the input model.
#' @export
#'
#' @examples
#' # Loading a test dataset into the global environment
#' x <- stackloss
#'
#' # Quick low dimensional clustering with K = 3 clusters and S = 1 dimensions
#' clust <- adproclus_low_dim(x, 3, 1)
#'
#' # Plot the matrix B', connecting components with variables
#' plot_vars_by_comp(clust)
plot_vars_by_comp <- function(model,
                              title = "B' of Low Dimensional
                              ADPROCLUS Solution") {
  checkmate::assertClass(model, "adpc")
  checkmate::assertString(title, null.ok = TRUE)

  if (is.null(title)) {
    title <- "B' of Low Dimensional ADPROCLUS Solution"
  }
  if (is.null(model$C)) {
    stop("Model must be a low dimensional ADPROCLUS solution.")
  }
  corrplot::corrplot(t(model$B),
    is.corr = FALSE, title = title,
    mar = c(0, 0, 2, 0)
  )
  invisible(model)
}

# Helper function
#' Calculate the number of observations in the overlap of two clusters
#'
#' @param edge Pair of cluster numbers.
#' @param A Cluster membership matrix.
#'
#' @return Number of observations that are in both clusters simultaneously.
#'
#' @noRd
.extract_overlap <- function(edge, A) {
  cluster1 <- edge[1]
  cluster2 <- edge[2]
  # no. of rows in which both clusters (cols) are 1
  overlap <- nrow(A[A[, cluster1] == 1 & A[, cluster2] == 1, ])
  overlap
}
