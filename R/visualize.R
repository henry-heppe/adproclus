#Visualization functions for an ADPROCLUS model

#' Plots a low dimensional ADPROCLUS solution as a network
#'
#' @param model the model to visualized
#' @param cluster_names optional: the names of the clusters
#' @param relative_overlap optional: choose whether you want to absolute number of observations on the edges
#' @param vertex_color optional: color of vertices
#' @param edge_color_low optional: lower end of color spectrum for edges
#' @param edge_color_high optional: upper end of color spectrum for edges
#'
#' @return nothing, outputs a plot
#' @export
#'
#' @examples
#' # Loading a test dataset into the global environment
#' x <- ADPROCLUS::CGdata
#'
#' # Quick low dimensional clustering with K = 3 clusters and S = 1 dimensions
#' clust <- adproclusLD(x, 3, 1)
#'
#' # Plot the overlapping the clusters
#' plot_clusters_as_network(clust)
plot_clusters_as_network <- function(model, cluster_names = NULL,
                                     relative_overlap = TRUE,
                                     vertex_color = "antiquewhite",
                                     edge_color_low = "lavenderblush2",
                                     edge_color_high = "lavenderblush4") {

        withr::local_seed(1)

        data <- model$Model
        A <- model$A
        C <- model$C
        B <- model$B

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

.extract_overlap <- function(edge, A) {
        cluster1 <- edge[1] #vertex on one side of the edge
        cluster2 <- edge[2] #vertex on the other side of the edge
        overlap <- nrow(A[A[,cluster1] == 1 & A[,cluster2] == 1,]) #no. of rows in which both clusters (cols) are 1

}
