# S3 methods for ADPROCLUS solution representation, printing and plotting

#' Constructor for a (low dimensional) ADPROCLUS solution object
#'
#' Yields an object of class \code{adpc}, which can be printed, plotted and
#' summarized by the corresponding methods. Mandatory input are the membership
#' matrix \eqn{\boldsymbol{A}} and the profile matrix \eqn{\boldsymbol{P}}
#' (where the number of columns from \eqn{\boldsymbol{A}} corresponds to
#' the number of rows in \eqn{\boldsymbol{P}}),
#' if the object is to represent a full dimensional ADPROCLUS model.
#' For a low dimensional ADPROCLUS model, the matrices \eqn{\boldsymbol{C}}
#' and \eqn{\boldsymbol{B}} have to be provided and \eqn{\boldsymbol{P}} can
#' be inferred from those. All other inputs are optional but may be included
#' so that the output from the \code{summary(), print(), plot()} is complete.
#' For further details on the (low dimensional) ADPROCLUS model a
#' nd what every element of the objects means
#' see \code{\link{adproclus}} and \code{\link{adproclus_low_dim}}.
#'
#' @param A Membership matrix A.
#' @param P Profile matrix P.
#' @param sse Sum of Squared Error.
#' @param totvar Total variance.
#' @param explvar Explained variance.
#' @param iterations Number of iterations.
#' @param timer Time run the complete algorithm.
#' @param timer_one_run Time to complete this single algorithm start.
#' @param initial_start  List of \code{type} of start and
#' \code{start_allocation} matrix
#' @param C Low dimensional profiles matrix C.
#' @param B Matrix of base vectors connecting low dimensional components with
#' original variables B.
#' @param runs List of suboptimal models.
#' @param parameters List of algorithm parameters.
#'
#' @return Object of class \code{adpc}.
#' @export
#'
#' @examples
#' # Create the information needed for a minimal object of class adpc
#' x <- stackloss
#' result <- adproclus(x, 3)
#' A <- result$A
#' P <- result$P
#'
#' # Use constructor to obtain object of class adpc
#' result_object <- adpc(A, P)
#'
adpc <- function(A, P,
                 sse = NULL, totvar = NULL, explvar = NULL,
                 iterations = NULL, timer = NULL, timer_one_run = NULL,
                 initial_start = NULL,
                 C = NULL, B = NULL,
                 runs = NULL, parameters = NULL) {
  checkmate::assert_matrix(A, any.missing = FALSE)
  checkmate::assert_matrix(P, any.missing = FALSE)
  checkmate::assert_number(sse, null.ok = TRUE)
  checkmate::assert_number(totvar, null.ok = TRUE)
  checkmate::assert_number(explvar, null.ok = TRUE)
  checkmate::assert_count(iterations, null.ok = TRUE)
  checkmate::assert_number(timer, null.ok = TRUE)
  checkmate::assert_number(timer_one_run, null.ok = TRUE)
  checkmate::assert_list(initial_start,
    types = c("character", "matrix"),
    null.ok = TRUE
  )
  checkmate::assert_list(runs, null.ok = TRUE)
  checkmate::assert_list(parameters, null.ok = TRUE)

  model_lowdim <- NULL

  if (!is.null(C) || !is.null(B)) {
    stopifnot(!is.null(C))
    stopifnot(!is.null(B))
    checkmate::assert_matrix(B, any.missing = FALSE)
    checkmate::assert_matrix(C, any.missing = FALSE)
    model_lowdim <- A %*% C

    if (!isTRUE(all.equal(P, C %*% t(B)))) {
      stopifnot(ncol(A) == nrow(C %*% t(B)))
      P <- C %*% t(B)
      warning("Inferred P as CB', since they were not equal.")
    }
  } else {
    checkmate::assert_matrix(P, any.missing = FALSE)
    stopifnot(ncol(A) == nrow(P))
  }


  object <- list(
    model = A %*% P, model_lowdim = model_lowdim, A = A, P = P,
    sse = sse, totvar = totvar, explvar = explvar,
    iterations = iterations, timer = timer, timer_one_run = timer_one_run,
    C = C, B = B,
    runs = runs, parameters = parameters
  )
  class(object) <- "adpc"
  object
}

#' Summary of ADPROCLUS solution
#'
#' For an object of class \code{adpc} as input, this method yields a summary
#' object of class \code{summary.adpc} including group characteristics of the
#' clusters in the solution. Works for both full and low dimensional solutions.
#' Adjust the parameters \code{digits, matrix_rows, matrix_cols} to change the
#' level of detail for the printing of the summary.
#'
#' @param object ADPROCLUS solution (class: \code{adpc}). Low dimensional model
#' possible.
#' @param title String. Default: "ADPROCLUS solution"
#' @param digits Integer. The number of digits that all decimal numbers will be
#' rounded to.
#' @param matrix_rows Integer. The number of matrix rows to display. OPTIONAL
#' @param matrix_cols Integer. The number of matrix columns to display. OPTIONAL
#' @param ... ignored
#'
#' @return Invisibly returns object of class \code{summary.adpc}.
#' @export
#'
#' @examples
#' # Obtain data, compute model, summarize model
#' x <- stackloss
#' model <- adproclus(x, 3)
#' model_summary <- summary(model)
summary.adpc <- function(object,
                         title = "ADPROCLUS solution",
                         digits = 3, matrix_rows = 10, matrix_cols = 5,
                         ...) {
  checkmate::assert_class(object, "adpc")
  checkmate::assert_string(title)
  checkmate::assert_int(digits, lower = 1, coerce = TRUE)
  checkmate::assert_int(matrix_rows, lower = 1, coerce = TRUE)
  checkmate::assert_int(matrix_cols, lower = 1, coerce = TRUE)

  A <- object$A
  k <- ncol(A)
  cluster_sizes_overlaps <- matrix(rep(0, k^2), k, k)

  for (i in 1:k) {
    for (j in 1:k) {
      cluster_sizes_overlaps[i, j] <- length(which(A[, i] == 1 & A[, j] == 1, ))
      cluster_sizes_overlaps[j, i] <- cluster_sizes_overlaps[i, j]
    }
  }
  cluster_characteristics <- list()
  if (is.null(object$C)) {
    for (i in 1:k) {
      members <- which(as.logical(A[, i]))
      cluster_characteristics <- append(
        cluster_characteristics,
        list(summary(object$model[members, , drop = FALSE])[c(1, 4, 6), , drop = FALSE])
      )
      names(cluster_characteristics)[i] <- colnames(A)[i]
    }
  } else {
    for (i in 1:k) {
      members <- which(as.logical(A[, i]))
      cluster_characteristics <- append(
        cluster_characteristics,
        list(summary(object$model_lowdim[members, , drop = FALSE])[c(1, 4, 6), , drop = FALSE])
      )
      names(cluster_characteristics)[i] <- colnames(A)[i]
    }
  }


  summary_res <- list(
    model_complete = object,
    cluster_sizes_overlaps = cluster_sizes_overlaps,
    cluster_characteristics = cluster_characteristics
  )
  print_settings <- list(
    digits = digits,
    matrix_rows = matrix_rows,
    matrix_cols = matrix_cols,
    title = title
  )
  summary_res <- append(summary_res, list(print_settings = print_settings))
  class(summary_res) <- "summary.adpc"
  summary_res
}

#' Print (low dimensional) ADPROCLUS summary
#'
#' Prints an object of class \code{summary.adpc} to represent and summarize a
#' (low dimensional) ADPROCLUS solution. A number of parameters for how the
#' results should be printed can be passed as an argument to
#' \code{summary.adpc()} which then passes it on to this method. This method
#' does not take a model of class \code{adpc} directly as input.
#'
#' @param x Object of class \code{summary.adpc}
#' @param ... ignored
#'
#' @return Invisibly returns object of class \code{summary.adpc}.
#' @export
#'
#' @examples
#' # Obtain data, compute model, print summary of model
#' x <- stackloss
#' model <- adproclus(x, 3)
#' print(summary(model))
print.summary.adpc <- function(x, ...) {
  checkmate::assert_class(x, "summary.adpc")

  # limit number of variables to print for cluster summary stats
  if (is.null(x$model_complete$C)) {
    n_var_true <- ncol(x$model_complete$model)
  } else {
    n_var_true <- ncol(x$model_complete$model_lowdim)
  }


  n_var_inc <- min(x$print_settings$matrix_cols, n_var_true)
  print(x$model_complete,
    digits = x$print_settings$digits,
    matrix_rows = x$print_settings$matrix_rows,
    matrix_cols = x$print_settings$matrix_cols
  )
  cat("Cluster sizes and overlaps:\n")
  print(x$cluster_sizes_overlaps)
  cat(" (diagonal entries: number of observations in a cluster)\n")
  cat(" (off-diagonal entry [i,j]:  number of observations both in cluster i and j)\n\n")
  if (is.null(x$model_complete$C)) {
    cat("Summary statistics of model variables per cluster:\n")
    if (n_var_true > n_var_inc) {
      cat("[", n_var_true - n_var_inc, "variables per cluster were omitted ]\n")
    }
  } else {
    cat("Summary statistics of low dimensional components per cluster:\n")
    if (n_var_true > n_var_inc) {
      cat(
        "[", n_var_true - n_var_inc,
        "components per cluster were omitted ]\n"
      )
    }
  }

  lapply(
    seq_len(ncol(x$model_complete$A)),
    function(i) {
      cat(names(x$cluster_characteristics)[i], "\n")
      print(x$cluster_characteristics[[i]][, 1:n_var_inc, drop = FALSE])
    }
  )


  invisible(x)
}


#' Plotting a (low dimensional) ADPROCLUS solution
#'
#' When passing a (low dimensional) ADPROCLUS solution of class \code{adpc} to
#' the generic \code{plot()}, this method plots the solution in one of the
#' following three ways:
#' \describe{
#' \item{\code{Network}}{Each cluster is a vertex and
#' the edge between two vertices represents the overlap between the
#' corresponding clusters. The size of a vertex corresponds to the cluster size.
#' The overlap is represented through color, width and numerical label of the
#' edge. The numerical edge-labels can be relative
#' (number of overlap observations / total observations)
#' or absolute (number of observations in both clusters).}
#' \item{\code{Profiles}}{Plot the profile matrix (\eqn{\boldsymbol{P}}
#' for full dimensional model, \eqn{\boldsymbol{C}} for low dimensional model)
#' in the style of a correlation plot to visualize the relation of each cluster
#' with each variable.}
#' \item{code{Low Dimensional Components}}{Plot the low dimensional
#' component-by-variable matrix \eqn{\boldsymbol{B'}} in the style of a
#' correlation plot to visualize the relation of each component with each
#' original variable. \strong{NOTE:} Only works for low dimensional ADPROCLUS.}
#' }
#'
#' @param x Object of class \code{adpc}. (Low dimensional) ADPROCLUS solution
#' @param type Choice for type of plot: one of "Network", "Profiles",
#' "VarsByComp". Default: "Network".
#' @param title String. OPTIONAL.
#' @param relative_overlap Logical, only applies to plot of
#' \code{type = "Network"}. If \code{TRUE} (default), the number of observations
#' belonging to two clusters is divided by the total number of observations.
#' @param ... additional arguments will be passed on to the functions
#' \code{plot_cluster_network(), plot_profiles(), plot_vars_by_comp()}
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
#' # Produce three plots of the model
#' plot(clust, type = "Network")
#' plot(clust, type = "Profiles")
#' plot(clust, type = "vars_by_comp")
plot.adpc <- function(x,
                      type = "Network",
                      title = NULL,
                      relative_overlap = TRUE,
                      ...) {
  checkmate::assertClass(x, "adpc")
  checkmate::assertChoice(type, c("Network", "Profiles", "vars_by_comp"))
  checkmate::assertString(title, null.ok = TRUE)
  checkmate::assertFlag(relative_overlap)

  # Check for illegal choice of vars_by_comp for full dim x is in plotVarsByComp()
  if (type == "vars_by_comp") {
    plot_vars_by_comp(model = x, title = title, ...)
  } else if (type == "Profiles") {
    plot_profiles(model = x, title = title, ...)
  } else {
    plot_cluster_network(
      model = x, title = title,
      relative_overlap = relative_overlap,
      ...
    )
  }

  invisible(x)
}



#' Print basic information on ADPROCLUS solution
#'
#' For an object of class \code{adpc} as input, this method prints basic
#' information about the ADPROCLUS solution represented by the object.
#' Works for both full and low dimensional solutions. Adjust the parameters
#' \code{digits, matrix_rows, matrix_cols}
#' to change the level of detail printed.
#'
#' @param x ADPROCLUS solution (class: \code{adpc})
#' @param title String. Default: "ADPROCLUS solution"
#' @param digits Integer. The number of digits that all decimal numbers will
#' be rounded to.
#' @param matrix_rows Integer. The number of matrix rows to display. OPTIONAL
#' @param matrix_cols Integer. The number of matrix columns to display. OPTIONAL
#' @param ... ignored
#'
#' @export
#'
#' @examples
#' # Obtain data, compute model, print model
#' x <- stackloss
#' model <- adproclus(x, 3)
#' print(model)
print.adpc <- function(x,
                       title = "ADPROCLUS solution",
                       digits = 3,
                       matrix_rows = 10, matrix_cols = 15,
                       ...) {
  checkmate::assert_class(x, "adpc")
  checkmate::assert_string(title)
  checkmate::assert_int(digits, lower = 1, coerce = TRUE)
  checkmate::assert_int(matrix_rows, lower = 1, coerce = TRUE)
  checkmate::assert_int(matrix_cols, lower = 1, coerce = TRUE)

  n_obs_true <- nrow(x$model)
  n_obs_inc <- min(matrix_rows, n_obs_true)
  n_clust_true <- ncol(x$A)
  n_clust_inc_row <- min(matrix_rows, n_clust_true)
  n_clust_inc_col <- min(matrix_cols, n_clust_true)
  n_var_true <- ncol(x$model)
  n_var_inc <- min(matrix_cols, n_var_true)
  n_randomstart <- x$parameters$nrandomstart
  n_semirandomstart <- x$parameters$nsemirandomstart
  start_allocation <- x$parameters$start_allocation

  if (!is.null(x$C)) {
    cat("Low Dimensional", title, "\n")
    cat("   number of clusters:", ncol(x$A), "\n")
    cat("   number of components: ", ncol(x$C), "\n")
    cat("   data format: ", nrow(x$model), "x", ncol(x$model), "\n")
    cat("   number of total starts:",
      n_randomstart + n_semirandomstart + 1 * !is.null(start_allocation),
      "\n"
    )
    if (!is.null(start_allocation)) {
      cat("   A rational start was also included.\n")
    }
    cat("Results:", "\n")
    cat("   explained variance:", round(x$explvar, digits), "\n")
    cat("   processing time:", round(x$timer, digits), "s", "\n")
    cat("   iterations:", x$iterations, "\n")
    cat("\n")
    cat("A (cluster membership matrix):", "\n")
    print(x$A[1:n_obs_inc, 1:n_clust_inc_col, drop = FALSE])
    if (n_obs_true - n_obs_inc > 0) {
      cat("[", n_obs_true - n_obs_inc, " rows were omitted ]\n")
    }
    if (n_clust_true - n_clust_inc_col > 0) {
      cat("[", n_clust_true - n_clust_inc_col, " columns were omitted ]\n")
    }
    cat("\n")
    cat("C (profiles in terms of components - cluster by component):", "\n")
    print(round(x$C[1:n_clust_inc_row, , drop = FALSE], digits))
    if (n_clust_true - n_clust_inc_row > 0) {
      cat("[", n_clust_true - n_clust_inc_row, " rows were omitted ]\n")
    }
    cat("\n")
    cat("B' (components by variables): ", "\n")
    print(round(t(x$B)[, 1:n_var_inc, drop = FALSE], digits))
    if (n_var_true - n_var_inc > 0) {
      cat("[", n_var_true - n_var_inc, " columns were omitted ]\n")
    }
  } else {
    cat(title, "\n")
    cat("Setup:", "\n")
    cat("   number of clusters:", ncol(x$A), "\n")
    cat("   data format: ", nrow(x$model), "x", ncol(x$model), "\n")
    cat("   number of total starts:",
      n_randomstart + n_semirandomstart + 1 * !is.null(start_allocation),
      "\n"
    )
    if (!is.null(start_allocation)) {
      cat("   A rational start was also included.\n")
    }
    cat("Results:", "\n")
    cat("   explained variance:", round(x$explvar, digits), "\n")
    cat("   processing time:", round(x$timer, digits), "s", "\n")
    cat("   iterations:", x$iterations, "\n")
    cat("A (cluster membership matrix):", "\n")
    print(x$A[1:n_obs_inc, 1:n_clust_inc_col, drop = FALSE])
    if (n_obs_true - n_obs_inc > 0) {
      cat("[", n_obs_true - n_obs_inc, " rows were omitted ]\n")
    }
    if (n_clust_true - n_clust_inc_col > 0) {
      cat("[", n_clust_true - n_clust_inc_col, " columns were omitted ]\n")
    }
    cat("P (profiles):", "\n")
    print(round(x$P[1:n_clust_inc_row, 1:n_var_inc, drop = FALSE], digits))
    if (n_clust_true - n_clust_inc_row > 0) {
      cat("[", n_clust_true - n_clust_inc_row, " rows were omitted ]\n")
    }
    if (n_var_true - n_var_inc > 0) {
      cat("[", n_var_true - n_var_inc, " columns were omitted ]\n")
    }
  }
}
