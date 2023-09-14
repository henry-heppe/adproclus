#S3 methods for ADPROCLUS solution

#' Constructor for a ADPROCLUS solution object
#'
#' Yields an object of class \code{adpc}.
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

}

#' Summarize important information on ADPROCLUS solution
#'
#' @param object ADPROCLUS solution (class: \code{adpc})
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
        for (i in 1:k) {
                members <- which(as.logical(A[,i]))
                cluster_characteristics <- append(cluster_characteristics, list(summary(data.frame(object$model[members,]))[c(1,4,6),]))
                names(cluster_characteristics)[i] <- paste("cluster_", i, sep = "")#issue: put cluster name here
        }

        summary_res <- list(model_complete = object,
                            cluster_sizes_overlaps = cluster_sizes_overlaps,
                            cluster_characteristics = cluster_characteristics)
        print_settings <- list(digits = digits, matrix_rows = matrix_rows, matrix_cols = matrix_cols, title = title)
        summary_res <- append(summary_res, list(print_settings = print_settings))
        class(summary_res) <- "summary.adpc"
        return(summary_res)
}

print.summary.adpc <- function(x, ...) {
        # object <- x$model_complete
        # if(!is.null(object$C)) {
        #         res <- paste("Low Dimensional ", title, "\n", sep = "")
        #         res <- paste(res, "Setup:", "\n")
        #         res <- paste(res, "   number of clusters: ", ncol(object$A), "\n")
        #         res <- paste(res, "   data format: ", nrow(object$model), "x", ncol(object$model), "\n")
        #         res <- paste(res, "   number of runs: ", object$parameters$nrandomstarts + object$parameters$nsemirandomstarts, "\n")
        #         res <- paste(res, "Results:", "\n")
        #         res <- paste(res, "   explained variance: ", object$explvar, "\n   time: ", object$timer, "\n   iterations: ", object$iterations, "\n")
        #
        #         res <- paste(res, "First few rows of A: ", object$A[1:5], "\n")
        #         res <- paste(res, "First few profiles in terms of components: ", object$C[1:10], "\n")
        #         res <- paste(res, "Components by variables matrix B: ", object$B, "\n")
        # } else {
        #         res <- title
        #         res <- paste(res, "Setup:", "\n")
        #         res <- paste(res, "   number of clusters:", ncol(object$A), "\n")
        #         res <- paste(res, "   number of dimensions:", ncol(object$C), "\n")
        #         res <- paste(res, "   data format: ", nrow(object$model), "x", ncol(object$model), "\n")
        #         res <- paste(res, "   number of runs: ", object$parameters$nrandomstarts + object$parameters$nsemirandomstarts, "\n")
        #         res <- paste(res, "Results:", "\n")
        #         res <- paste(res, "explained variance:", object$explvar, "time:", object$timer, "iterations:", object$iterations, "\n")
        #         res <- paste(res, "First few rows of A:", object$A[1:5], "\n")
        #         res <- paste(res, "First few profiles:", object$P[1:5], "\n")
        # }

}


plot.adpc <- function() {

}



#' Print basic information on ADPROCLUS solution
#'
#' For an object of class adpc as input, this method prints basic information
#' about the ADPROCLUS solution represented by the object. Works for both full and low dimensional solutions.
#'
#' @param x ADPROCLUS solution (class: \code{adpc})
#' @param title String. Default: "ADPROCLUS solution"
#' @param ... ignored
#' @param digits Integer. The number of digits that all decimal numbers will be rounded to.
#' @param matrix_rows Integer. The number of matrix rows to display. OPTIONAL
#' @param matrix_cols Integer. The number of matrix columns to display. OPTIONAL
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
                cat("   number of runs:", x$parameters$nrandomstarts + x$parameters$nsemirandomstarts, "\n")
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
                cat("   number of runs:", x$parameters$nrandomstarts + x$parameters$nsemirandomstarts, "\n")
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
