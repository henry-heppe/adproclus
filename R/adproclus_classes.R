#S3 methods for ADPROCLUS solution representation, printing and plotting

#' Constructor for a (low dimensional) ADPROCLUS solution object
#'
#' Yields an object of class \code{adpc}, which can be printed, plotted and summarized by the corresponding methods.
#'
#' @param model model
#' @param A A
#' @param P P
#' @param sse sse
#' @param totvar totvar
#' @param explvar explvar
#' @param iterations iterations
#' @param timer timer
#' @param initialStart initialstart
#' @param C C
#' @param B B
#' @param runs runs
#'
#' @return Object of class \code{adpc}
#' @export
#'
#' @examples
#' #Create the information needed for an object of class adpc
#' x <- ADPROCLUS::CGdata[1:100,]
#' result <- adproclus(x, 3)
#' model <- result$model
#' A <- result$A
#' P <- result$P
#' totvar <- result$totvar
#' explvar <- result$explvar
#' iterations <- result$iterations
#' timer <- result$timer
#' initialStart <- result$ initialStart
#'
#' #Create object of class adpc
#' result_object <- adpc(model, A, P, sse, totvar, explvar, iterations, timer, initialStart)
#'
adpc <- function(model, A, P, sse, totvar, explvar, iterations, timer, initialStart, C = NULL, B = NULL, runs = NULL) {
        stopifnot(is.matrix(model))
        stopifnot(is.matrix(A))
        stopifnot(is.matrix(P))
        object <- list(model = model, A = A, P = P)
        class(object) <- "adpc"
        return(object)
        #issue: not done
}

#' Summarize important information on ADPROCLUS solution
#'
#' @param object ADPROCLUS solution (class: \code{adpc}). Low dimensional model possible.
#' @param title String. Default: "ADPROCLUS solution"
#' @param digits Integer. The number of digits that all decimal numbers will be rounded to.
#' @param matrix_rows Integer. The number of matrix rows to display. OPTIONAL
#' @param matrix_cols Integer. The number of matrix columns to display. OPTIONAL
#' @param ... ignored
#'
#' @return invisibly returns object of class \code{summary.adpc}
#' @export
#'
#' @examples
#' x <- ADPROCLUS::CGdata[1:100,]
#' model <- adproclus(x, 3)
#' summary(model)
summary.adpc <- function(object, title = "ADPROCLUS solution", digits = 3, matrix_rows = 10, matrix_cols = 15, ...) {
        #members per cluster
        #cluster overlaps
        #mean, st.dev, min, max, percentiles per cluster (use summary() on vector and select the relevant numbers)
        #inter-cluster variances?
        A <- object$A
        k <- ncol(A)
        cluster_sizes <- colSums(object$A)
        cluster_sizes_overlaps <- matrix(rep(0, k^2), k, k)

        for (i in 1:k) {
                for (j in 1:k) {
                        cluster_sizes_overlaps[i,j] <- nrow(A[A[,i] == 1 & A[,j] == 1,])
                }
        }
        cluster_characteristics <- list()
        if (is.null(object$C)) {
                n_var_true <- ncol(object$model) #limit number of variables to include in summary statistics
                n_var_inc <- min(matrix_cols, n_var_true)
                for (i in 1:k) {
                        members <- which(as.logical(A[,i]))
                        cluster_characteristics <- append(cluster_characteristics, list(summary(data.frame(object$model[members,1:n_var_inc]))[c(1,4,6),]))
                        names(cluster_characteristics)[i] <- paste("cluster_", i, sep = "")#issue: put cluster name here
                }
        } else {
                n_var_true <- ncol(object$model_lowdim) #limit number of variables to include in summary statistics
                n_var_inc <- min(matrix_cols, n_var_true)
                for (i in 1:k) {
                        members <- which(as.logical(A[,i]))
                        cluster_characteristics <- append(cluster_characteristics, list(summary(data.frame(object$model_lowdim[members,1:n_var_inc]))[c(1,4,6),]))
                        names(cluster_characteristics)[i] <- paste("cluster_", i, sep = "")#issue: put cluster name here
                }
        }


        summary_res <- list(model_complete = object,
                            cluster_sizes_overlaps = cluster_sizes_overlaps,
                            cluster_characteristics = cluster_characteristics)
        print_settings <- list(digits = digits, matrix_rows = matrix_rows, matrix_cols = matrix_cols, title = title)
        summary_res <- append(summary_res, list(print_settings = print_settings))
        class(summary_res) <- "summary.adpc"
        return(summary_res)
}

#' Printing (low dim) ADPROCLUS summary of class \code{summary.adpc}
#'
#' Prints a set of summary statistics for a (low dimensional) ADPROCLUS solution as provided by the function \code{summary.adpc()}.
#' A number of parameters for how the results should be printed can be passed as an argument to \code{summary.adpc()}
#' which then passes it on to this method.
#'
#' @param x Object of class \code{summary.adpc}
#' @param ... ignored
#'
#' @return \code{summary.adpc} object invisibly
#' @export
#'
#' @examples
#' x <- ADPROCLUS::CGdata[1:100,]
#' model <- adproclus(x, 3)
#' summary(model)
print.summary.adpc <- function(x, ...) {
        n_var_true <- ncol(x$model_complete$model) #limit number of variables to include in summary statistics
        n_var_inc <- min(x$print_settings$matrix_cols, n_var_true)
        print(x$model_complete,
              digits = x$print_settings$digits,
              matrix_rows = x$print_settings$matrix_rows,
              matrix_cols = x$print_settings$matrix_cols)
        cat("Cluster sizes and overlaps:\n")
        print(x$cluster_sizes_overlaps)
        cat(" [diagonal entries: number of observations in a cluster]\n")
        cat(" [off-diagonal entry [i,j]:  number of observations both in cluster i and j]\n\n")
        if(is.null(x$model_complete$C)) {
                cat("Summary statistics of model variables per cluster:\n")
                if (n_var_true > n_var_inc) {
                        cat("[", n_var_true - n_var_inc, "variables were omitted ]\n")
                }
        } else {
                cat("Summary statistics of low dimensional components per cluster:\n")
                if (n_var_true > n_var_inc) {
                        cat("[", n_var_true - n_var_inc, "components were omitted ]\n")
                }
        }

        lapply(1:ncol(x$model_complete$A),
               function(i) {cat(names(x$cluster_characteristics)[i], "\n"); print(x$cluster_characteristics[[i]])})
        invisible(x)


}


#' Plotting a (low dimensional) ADPROCLUS solution
#'
#' When passing a (low dimensional) ADPROCLUS solution of class \code{adpc} to the generic \code{plot()},
#' this method plots the solution in one of the following three ways:
#' \describe{
#' \item{\code{Network}}{Each cluster is a vertex and
#' the edge between two vertices represents the overlap between the corresponding clusters.
#' The size of a vertex corresponds to the cluster size. The overlap is represented
#' through color, width and numerical label of the edge.
#' The numerical edge-labels can be relative (number of overlap observations / total observations)
#' or absolute (number of observations in both clusters).}
#' \item{\code{Profiles}}{Plot the profile matrix (\eqn{\boldsymbol{P}} for full dimensional model,
#' \eqn{\boldsymbol{C}} for low dimensional model)in the style of a correlation plot to visualize the
#' relation of each cluster with each variable.}
#' \item{code{Low Dimensional Components}}{Plot the low dimensional component-by-variable matrix \eqn{\boldsymbol{B'}}
#'  in the style of a correlation plot to visualize the relation of each component with each original variable.
#'  \strong{NOTE:} Only works for low dimensional ADPROCLUS.}
#' }
#'
#' @param x Object of class \code{adpc}. (Low dimensional) ADPROCLUS solution
#' @param type Choice for type of plot: one of "Network", "Profiles", "LowDimComponents". Default: "Network".
#' @param title String. OPTIONAL.
#' @param ... ignored
#'
#' @return Invisibly returns the input model.
#' @export
#'
#' @examples
#'#add example
plot.adpc <- function(x, type = "Network", title = NULL, ...) {
        checkmate::assertChoice(type, c("Network", "Profiles", "LowDimComponents"))
        #check for illegal choice of LowDimComponentes for full dim model is done in plotVarsByComp()
        if (type == "LowDimComponents") {
                if (!is.null(title)) {
                        plotVarsByComp(x, title)
                } else {
                        plotVarsByComp(x)
                }

        } else if (type == "Profiles") {
                if (!is.null(title)) {
                        plotProfiles(x, title)
                } else {
                        plotProfiles(x)
                }

        } else {
                if (!is.null(title)) {
                        plotClusterNetwork(x, title)
                } else {
                        plotClusterNetwork(x)
                }

        }

        invisible(x)
}



#' Print basic information on ADPROCLUS solution
#'
#' For an object of class adpc as input, this method prints basic information
#' about the ADPROCLUS solution represented by the object. Works for both full and low dimensional solutions.
#'
#' @param x ADPROCLUS solution (class: \code{adpc})
#' @param title String. Default: "ADPROCLUS solution"
#' @param digits Integer. The number of digits that all decimal numbers will be rounded to.
#' @param matrix_rows Integer. The number of matrix rows to display. OPTIONAL
#' @param matrix_cols Integer. The number of matrix columns to display. OPTIONAL
#' @param ... ignored
#'
#' @export
#'
#' @examples
#' x <- ADPROCLUS::CGdata[1:100,]
#' model <- adproclus(x, 3)
#' print(model)
print.adpc <- function(x, title = "ADPROCLUS solution", digits = 3, matrix_rows = 10, matrix_cols = 15, ...) {
        n_obs_true <- nrow(x$model)
        n_obs_inc <- min(matrix_rows, n_obs_true)
        n_clust_true <- ncol(x$A)
        n_clust_inc_row <- min(matrix_rows, n_clust_true)
        n_clust_inc_col <- min(matrix_cols, n_clust_true)
        n_var_true <- ncol(x$model)
        n_var_inc <- min(matrix_cols, n_var_true)

        if(!is.null(x$C)) {
                cat("Low Dimensional", title, "\n")
                cat("   number of clusters:", ncol(x$A), "\n")
                cat("   data format: ", nrow(x$model), "x", ncol(x$model), "\n")
                cat("   number of (semi-) random starts:", x$parameters$nrandomstart + x$parameters$nsemirandomstart, "\n")
                if (!is.null(x$parameters$start_allocation)) {
                        cat("   A rational start was also included.\n")
                }
                cat("Results:", "\n")
                cat("   explained variance:", round(x$explvar, digits), "\n")
                cat("   processing time:", round(x$timer, digits), "s", "\n")
                cat("   iterations:", x$iterations, "\n")
                cat("A (cluster membership matrix):", "\n")
                print(round(x$A[1:n_obs_inc,1:n_clust_inc_col], digits))
                if (n_obs_true - n_obs_inc > 0) {
                        cat("[", n_obs_true - n_obs_inc, " rows were omitted ]\n")
                }
                if (n_clust_true - n_clust_inc_col > 0) {
                        cat("[", n_clust_true - n_clust_inc_col, " columns were omitted ]\n")
                }
                cat("C (profiles in terms of components - cluster by component):", "\n")
                print(round(x$C[1:n_clust_inc_row,], digits))
                if (n_clust_true - n_clust_inc_row > 0) {
                        cat("[", n_clust_true - n_clust_inc_row, " rows were omitted ]\n")
                }
                cat("B' (components by variables): ", "\n")
                print(round(t(x$B)[,1:n_var_inc], digits))
                if (n_var_true - n_var_inc > 0) {
                        cat("[", n_var_true - n_var_inc, " columns were omitted ]\n")
                }
        } else {
                cat(title, "\n")
                cat("Setup:", "\n")
                cat("   number of clusters:", ncol(x$A), "\n")
                cat("   data format: ", nrow(x$model), "x", ncol(x$model), "\n")
                cat("   number of (semi-) random starts:", x$parameters$nrandomstart + x$parameters$nsemirandomstart, "\n")
                if (!is.null(x$parameters$start_allocation)) {
                        cat("   A rational start was also included.\n")
                }
                cat("Results:", "\n")
                cat("   explained variance:", round(x$explvar, digits), "\n")
                cat("   processing time:", round(x$timer, digits), "s", "\n")
                cat("   iterations:", x$iterations, "\n")
                cat("A (cluster membership matrix):", "\n")
                print(round(x$A[1:n_obs_inc,1:n_clust_inc_col], digits))
                if (n_obs_true - n_obs_inc > 0) {
                        cat("[", n_obs_true - n_obs_inc, " rows were omitted ]\n")
                }
                if (n_clust_true - n_clust_inc_col > 0) {
                        cat("[", n_clust_true - n_clust_inc_col, " columns were omitted ]\n")
                }
                cat("P (profiles):", "\n")
                print(round(x$P[1:n_clust_inc_row,1:n_var_inc], digits))
                if (n_clust_true - n_clust_inc_row > 0) {
                        cat("[", n_clust_true - n_clust_inc_row, " rows were omitted ]\n")
                }
                if (n_var_true - n_var_inc > 0) {
                        cat("[", n_var_true - n_var_inc, " columns were omitted ]\n")
                }
        }
}

#' Naming the clusters of an ADPROCLUS solution
#'
#' @param model ADPROCLUS solution (class: \code{adpc})
#' @param cluster_names Names of the clusters as list.
#'
#' @return The input model  of class \code{adpc} with adjusted cluster names.
#' @export
#'
#' @examples
#' #example
name_clusters_adpc <- function(model, cluster_names) {
        result <- model
        checkmate::assertVector(cluster_names, len = ncol(model$A))
        #issue: implement the naming here
        return(result)
}

#' Naming the components of a low dimensional ADPROCLUS solution
#'
#' @param model ADPROCLUS solution (class: \code{adpc})
#' @param component_names Names of the components of the low dimensional model as list.
#' Can be used to name variables for full dimensional model.
#'
#' @return The input model  of class \code{adpc} with adjusted component names.
#' @export
#'
#' @examples
#' #example
name_components_adpc <- function(model, component_names) {
        result <- model
        if(is.null(model$C)) {
                checkmate::assertVector(component_names, len = ncol(model$P))
        } else {
                checkmate::assertVector(component_names, len = ncol(model$C))
        }
        #issue: implement the naming here
        return(result)
}
