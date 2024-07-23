# Visualization functions for an ADPROCLUS model

if(getRversion() >= "2.15.1")  utils::globalVariables(c("clusters", "components"))

#' Scree plot of (low dimensional) ADPROCLUS models
#'
#' Used for scree-plot based model selection. Visualizes a set of ADPROClUS models
#' in terms of their number of clusters and model fit (SSE or unexplained variance).
#' For low dimensional ADPROCLUS models plots are made with the number of
#' components on the x-axis for each given number of clusters. One can then
#' choose to have them displayed all in one plot (\code{grid = FALSE}) or next
#' to each other in separate plots (\code{grid = TRUE}).
#'
#'
#' @param model_fit Matrix of SSE or unexplained variance scores as given by the
#' output of \code{\link{mselect_adproclus}} or
#' \code{\link{mselect_adproclus_low_dim}}.
#' @param title String. Title to be displayed in plot.
#' @param grid Boolean. \code{FALSE} means for low dimensional ADPROCLUS all
#' lines will be in one plot. \code{TRUE} means separate plots.
#' @param digits Integer. The number of decimal places to display.
#'
#' @return Invisibly returns the \code{ggplot2} object.
#' @export
#'
#' @examples
#' # Loading a test dataset into the global environment
#' x <- stackloss
#'
#' # Estimating models with cluster parameter values ranging from 1 to 4
#' model_fits <- mselect_adproclus(data = x, min_nclusters = 1, max_nclusters = 4, seed = 1)
#'
#' # Plot the results as a scree plot to select the appropriate number of clusters
#' plot_scree_adpc(model_fits)
#'
#' # Estimating models with cluster parameter values ranging from 1 to 4
#' # and component parameter values also ranging from 1 to 4
#' model_fits <- mselect_adproclus_low_dim(data = x, 1, 4, 1, 4, seed = 1)
#'
#' # Plot the results as a scree plot to select the appropriate number of clusters
#' plot_scree_adpc(model_fits)
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{mselect_adproclus}}}{to obtain the \code{model_fit} input from the possible ADPROCLUS models}
#'   \item{\code{\link{mselect_adproclus_low_dim}}}{to obtain the \code{model_fit} input from the possible low dimensional ADPROCLUS models}
#'   \item{\code{\link{select_by_CHull}}}{for automatic model selection via CHull method}
#' }
plot_scree_adpc <- function(model_fit,
                          title = "Scree plot of ADPROCLUS models",
                          grid = FALSE,
                          digits = 3) {
        checkmate::assert_matrix(model_fit)
        checkmate::assert_string(title)
        checkmate::assert_count(digits, positive = TRUE, coerce = TRUE)

        model_fit <- round(model_fit, digits)

        #option for full dimensional model
        if (ncol(model_fit) == 1) {
                data <- data.frame(clusters = strtoi(rownames(model_fit)), model_fit[, 1])
                colnames(data)[2] <- colnames(model_fit)
                fit_var <- colnames(model_fit)
                scree_plot <- ggplot2::ggplot(data, ggplot2::aes(x = clusters, y = !!(rlang::ensym(fit_var)))) +
                        ggplot2::geom_line() +
                        ggplot2::geom_point() +
                        ggplot2::labs(x = "Number of Clusters", y = gsub("_", " ", fit_var), title = title) +
                        ggplot2::scale_x_continuous(breaks = scales::breaks_extended(nrow(model_fit))) +
                        ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 10^(-digits))) +
                        ggplot2::theme_classic()
        } else {
                #option for low dimensional model
                if (substring(colnames(model_fit)[[1]], 1, 1) == "S") {
                        fit_var <- "SSE"
                } else {
                        fit_var <- "Unexplained_Variance"
                }

                data_wide <- data.frame(clusters = strtoi(rownames(model_fit)), model_fit)
                data <- tidyr::pivot_longer(data = data_wide,
                                                 cols = !clusters,
                                                 names_to = "components",
                                                 names_transform = readr::parse_number,
                                                 values_to = fit_var)
                data$clusters <- as.factor(data$clusters)

                if (grid) {
                        scree_plot <- ggplot2::ggplot(data, ggplot2::aes(x = components, y = !!(rlang::ensym(fit_var)), group = 1)) +
                                ggplot2::geom_line() +
                                ggplot2::geom_point() +
                                ggplot2::labs(x = "Number of Components", y = gsub("_", " ", fit_var), title = title) +
                                ggplot2::scale_x_continuous(breaks = scales::breaks_extended(nrow(model_fit))) +
                                ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 10^(-digits))) +
                                ggplot2::facet_wrap(ggplot2::vars(clusters), labeller = "label_both") +
                                ggplot2::theme_classic()
                } else {
                        scree_plot <- ggplot2::ggplot(data, ggplot2::aes(x = components, y = !!(rlang::ensym(fit_var)), color = clusters)) +
                                ggplot2::geom_line() +
                                ggplot2::geom_point() +
                                ggplot2::scale_color_manual(values = c("blue", "red", "green", "black")) +
                                ggplot2::labs(x = "Number of Components", y = gsub("_", " ", fit_var), title = title) +
                                ggplot2::scale_x_continuous(breaks = scales::breaks_extended(nrow(model_fit))) +
                                ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 10^(-digits))) +
                                ggplot2::theme_classic()
                }
        }
        suppressWarnings(print(scree_plot))
        invisible(scree_plot)
}

#' Scree plot of a pre-selection of low dimensional ADPROCLUS models
#'
#' To be used when one has selected a number of components for each number
#' of clusters. Plots the remaining sets of models to compare SSE or unexplained
#' variances. The input \code{model_fit} is supposed to be the output from the
#' \code{\link{select_by_CHull}} function applied to the output from
#' the \code{\link{mselect_adproclus_low_dim}} function.
#'
#' @param model_fit Matrix with SSE or unexplained variance values.
#' Can be obtained from \code{\link{select_by_CHull}}.
#' @param title String. Title to be displayed in plot.
#' @param digits Integer. The number of decimal places to display.
#'
#' @return Returns the \code{ggplot2} object.
#' @export
#'
#' @examples
#' # Loading a test dataset into the global environment
#' x <- stackloss
#'
#' # Estimating models with cluster parameter values ranging from 1 to 4
#' # and component parameter values also ranging from 1 to 4
#' model_fits <- mselect_adproclus_low_dim(data = x, 1, 4, 1, 4, seed = 1)
#'
#' # Choosing for each number of cluster the best number of components
#' model_fits_preselected <- select_by_CHull(model_fits)
#'
#' # Plot the results as a scree plot to select the appropriate number of clusters
#' plot_scree_adpc_preselected(model_fits_preselected)
plot_scree_adpc_preselected <- function(model_fit,
                             title = "Scree plot of low dim ADPROCLUS pre-selected models",
                             digits = 3) {
        checkmate::assert_matrix(model_fit)
        checkmate::assert_string(title)
        checkmate::assert_count(digits, positive = TRUE, coerce = TRUE)

        model_fit <- round(model_fit, digits)

        fit_var <- colnames(model_fit)[3]
        scree_plot <- ggplot2::ggplot(data.frame(model_fit), ggplot2::aes(x = clusters, y = !!(rlang::ensym(fit_var)), label = rownames(model_fit))) +
                ggplot2::geom_line() +
                ggplot2::geom_point() +
                ggrepel::geom_label_repel(box.padding = 0.5, segment.linetype = 0) +
                ggplot2::labs(x = "Number of Clusters", y = gsub("_", " ", fit_var), title = title) +
                ggplot2::scale_x_continuous(breaks = scales::breaks_extended(nrow(model_fit))) +
                ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 10^(-digits))) +
                ggplot2::theme_classic()
        scree_plot
}

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
#' @param filename Optional. Name of the file without extension.
#' @param ... Additional arguments passing to the
#' \code{qgraph::qgraph()} function, to customize the graph visualization.
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
                                 title = "Cluster network of ADPROCLUS solution",
                                 relative_overlap = TRUE,
                                 filetype = NULL,
                                 filename = NULL,
                                 ...) {
  checkmate::assertClass(model, "adpc")
  checkmate::assertString(title, null.ok = TRUE)
  checkmate::assertFlag(relative_overlap)
  checkmate::assertChoice(filetype, c(
          "R", "pdf", "svg", "tex",
          "jpg", "tiff", "png", ""
  ),
  null.ok = TRUE
  )
  checkmate::assertString(filename, null.ok = TRUE)

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
      filename = filename,
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
      edgelist = TRUE,
      ...
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
      filename = filename,
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
      edgelist = TRUE,
      ...
    )
  }
  invisible(model)
}

#' Plot profile matrix of ADPROCLUS solution
#'
#' Produce a representation of profile matrix \eqn{P}
#' (or \eqn{C} for low dimensional solution) of an ADPROCLUS
#' solution of class \code{adpc}.
#' The plot displays the profiles in the style of a correlation plot.
#' \strong{NOTE:} This function can also be called through the
#' \code{plot(model, type = "Profiles")} function with model an object of
#' class \code{adpc}.
#'
#' @param model Object of class \code{adpc}. (Low dimensional) ADPROCLUS
#' solution
#' @param title String. Default: "Profiles of ADPROCLUS solution"
#' @param ... Additional arguments passing to the
#' \code{corrplot::corrplot()} function, to customize the plot.
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
plot_profiles <- function(model,
                          title = "Profiles of ADPROCLUS solution",
                          ...) {
        checkmate::assertClass(model, "adpc")
        checkmate::assertString(title, null.ok = TRUE)

        if (is.null(title)) {
                title <- "Profiles of ADPROCLUS solution"
        }
        if (is.null(model$C)) {
                corrplot::corrplot(model$P,
                                   is.corr = FALSE, title = title,
                                   mar = c(0, 0, 2, 0),
                                   ...
                )
        } else {
                if (title == "Profiles of ADPROCLUS solution") {
                        title <- "Low dim Profiles C of ADPROCLUS solution"
                }
                corrplot::corrplot(model$C,
                                   is.corr = FALSE, title = title,
                                   mar = c(0, 0, 2, 0),
                                   ...
                )
        }
        invisible(model)
}

#' Plot variable to component matrix of ADPROCLUS solution
#'
#' Produce a representation of variable to component matrix
#' \eqn{B'} of a \strong{low dimensional} ADPROCLUS solution
#' of class \code{adpc}. The plot displays the scores in the style of a
#' correlation plot.
#' \strong{NOTE:} This function can be called through the
#' \code{plot(model, type = "vars_by_comp")} function
#' with model an object of class \code{adpc}.
#'
#' @param model Object of class \code{adpc}. Must be \strong{Low dimensional}
#' ADPROCLUS solution
#' @param title String. Default: "B' of Low Dimensional ADPROCLUS Solution"
#' @param ... Additional arguments passing to the
#' \code{corrplot::corrplot()} function, to customize the plot
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
                              title = "B' of Low Dimensional ADPROCLUS Solution",
                              ...) {
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
                           mar = c(0, 0, 2, 0),
                           ...
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
        overlap <- nrow(A[A[, cluster1] == 1 & A[, cluster2] == 1, , drop = FALSE])
        overlap
}
