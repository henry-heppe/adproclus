#Visualization functions for an ADPROCLUS model

#helper function
.extract_overlap <- function(edge, A) {
        cluster1 <- edge[1] #vertex on one side of the edge
        cluster2 <- edge[2] #vertex on the other side of the edge
        overlap <- nrow(A[A[,cluster1] == 1 & A[,cluster2] == 1,]) #no. of rows in which both clusters (cols) are 1

}

#' Network plot of a (low dimensional) ADPROCLUS solution
#'
#' Produce a representation of a (low dimensional) ADPROCLUS solution, where each cluster is a vertex and
#' the edge between two vertices represents the overlap between the corresponding clusters.
#' The size of a vertex corresponds to the cluster size. The overlap is represented
#' through color, width and numerical label of the edge.
#' The numerical edge-labels can be relative (number of overlap observations / total observations)
#' or absolute (number of observations in both clusters).
#' \strong{NOTE:} This function can be called through the \code{plot(model, type = "Network")} function.
#'
#' @param model ADPROCLUS solution (class: \code{adpc}). Low dimensional model possible.
#' @param cluster_names Names of the clusters as list. OPTIONAL
#' @param component_names Names of the variables (components if low dim model) as list. OPTIONAL
#' @param relative_overlap Boolean value. If \code{TRUE} (default), the number of observations belonging to two clusters
#' is divided by the total number of observations. If \code{FALSE}
#' the number of observations in a cluster overlap will be displayed on the edges.
#' @param vertex_color Color of the vertices as string. OPTIONAL
#' @param edge_color_low Lower end of color spectrum for edges as string. OPTIONAL
#' @param edge_color_high Upper end of color spectrum for edges as string. OPTIONAL
#'
#' @return Invisibly returns the input model.
#' @export
#'
#' @examples
#' # Loading a test dataset into the global environment
#' x <- ADPROCLUS::CGdata[1:100,]
#'
#' # Quick low dimensional clustering with K = 3 clusters and S = 1 dimensions
#' clust <- adproclusLD(x, 3, 1)
#'
#' # Plot the overlapping the clusters
#' plotClusterNetwork(clust)
plotClusterNetwork <- function(model, cluster_names = NULL, component_names = NULL,
                                     relative_overlap = TRUE,
                                     vertex_color = "antiquewhite",
                                     edge_color_low = "lavenderblush2",
                                     edge_color_high = "lavenderblush4") {

        checkmate::assertClass(model, "adpc")
        checkmate::assertFlag(relative_overlap)
        checkmate::assertString(vertex_color)
        checkmate::assertString(edge_color_low)
        checkmate::assertString(edge_color_high)

        model_temp <- model

        if(!is.null(cluster_names)) {
                model_temp <- name_clusters_adpc(model, cluster_names)

        }

        if(!is.null(component_names)) {
                model_temp <- name_components_adpc(model, component_names)
        }

        withr::local_seed(1)

        data <- model_temp$Model
        A <- model_temp$A
        C <- model_temp$C
        B <- model_temp$B

        k <- ncol(A)

        adjacency_matrix <- matrix(rep(1,k^2),k,k)
        diag(adjacency_matrix) <- rep(0,k)
        network <- igraph::as.undirected(igraph::graph_from_adjacency_matrix(adjacency_matrix))
        edgelist <- data.frame(t(igraph::as_edgelist(network)))

        weights <- unlist(lapply(edgelist, .extract_overlap, A = A))


        edge_colors <- grDevices::colorRampPalette(c(edge_color_low, edge_color_high))(length(weights))
        edge_colors <- edge_colors[rank(weights)]

        sizes <- colSums(A)
        sizes <- 100 * (sizes / sum(sizes))

        igraph::V(network)$size <- sizes

        if(relative_overlap) {#relative overlap: overlapping obs / totals observations
                weights <- round(weights / nrow(A), digits = 3)
                igraph::E(network)$weights <- weights
                plot(network,
                     layout = igraph::layout_with_fr(network, weights = weights * 10000),
                     vertex.label = cluster_names,
                     vertex.color = vertex_color,
                     edge.label = weights,
                     edge.label.cex = 1,
                     edge.label.dist = 1,
                     edge.width = (weights^1.2)*100,
                     edge.color = edge_colors
                )

        } else {
                igraph::E(network)$weights <- weights
                plot(network,
                     layout = igraph::layout_with_fr(network, weights = weights * 10000 /sum(weights)),
                     vertex.label = cluster_names,
                     vertex.color = vertex_color,
                     edge.label = weights,
                     edge.label.cex = 1,
                     edge.label.dist = 1,
                     edge.width = ((weights/sum(weights))^1.2)*100,
                     edge.color = edge_colors
                )
        }

}

#' Plot profile matrix of ADPROCLUS solution
#'
#' @param model Object of class \code{adpc}. (Low dimensional) ADPROCLUS solution
#' @param title String. Default: "Profiles of ADPROCLUS solution"
#'
#' @return Invisibly returns the input model.
#' @export
#'
#' @examples
#' #add exmple
plotProfiles <- function(model, title = "Profiles of ADPROCLUS solution") {
        if(is.null(model$C)) {
                corrplot::corrplot(model$P, is.corr = FALSE, title = title)
        } else {
                if (title == "Profiles of ADPROCLUS solution") {
                        title = "Profiles of Low dimensional ADPROCLUS solution"
                }
                corrplot::corrplot(model$C, is.corr = FALSE, title = title)
        }
        invisible(model)

}

#' Plot variable to component matrix of ADPROCLUS solution
#'
#' @param model Object of class \code{adpc}. Must be \strong{Low dimensional} ADPROCLUS solution
#' @param title String. Default: "B' of Low Dimensional ADPROCLUS Solution"
#'
#' @return Invisibly returns the input model.
#' @export
#'
#' @examples
#' #add examples
plotVarsByComp <- function(model, title = "B' of Low Dimensional ADPROCLUS Solution") {
        if (is.null(model$C)) {
                stop("Model must be a low dimensional ADPROCLUS solution.")
        }
        corrplot::corrplot(t(model$B), is.corr = FALSE, title = title)
        invisible(model)
}
