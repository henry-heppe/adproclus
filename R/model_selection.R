

#' Model selection helper for ADPROCLUS
#'
#' Performs ADPROCLUS for the number of clusters from \code{min_nclusters} to \code{max_nclusters}.
#' This replaces the need to manually estimate multiple models to select the best
#' number of clusters and returns the results in a format compatible with
#' \code{\link{plot_scree_adpc}} to obtain a scree plot.
#' Output is also compatible with \code{\link{select_by_CHull}} to
#' automatically select a suitable number of clusters.
#' The compatibility with both functions is only given if
#' \code{return_models = FALSE}.
#'
#'
#' @param data Object-by-variable data matrix of class \code{matrix} or
#'   \code{data.frame}.
#' @param min_nclusters Minimum number of clusters to estimate.
#' @param max_nclusters Maximum number of clusters to estimate.
#' @param return_models Boolean. If \code{FALSE} a vector of model fit scores is
#' returned, which is compatible with the \code{\link{plot_scree_adpc}} function.
#' If \code{TRUE} the list of actually estimated models is returned.
#' @param unexplvar Boolean. If \code{TRUE} the model fit is specified in terms
#' of unexplained variance. Otherwise it will be specified in terms of
#' Sum of Squared Errors (SSE). This propagates through to the scree plots.
#' @param start_allocation Optional starting cluster membership matrix to be
#' passed to the ADPROCLUS procedure. See \code{\link{get_rational}} for more information.
#' @param nrandomstart Number of random starts computed for each model.
#' @param nsemirandomstart Number of semi-random starts computed for each model.
#' @param algorithm Character string "\code{ALS1}" or "\code{ALS2}" (default),
#'   denoting the type of alternating least squares algorithm. Can be
#'   abbreviated with "1" or "2".
#' @param save_all_starts Logical. If \code{TRUE} and \code{return_models = TRUE},
#' the results of all algorithm starts are returned.
#' By default, only the best solution is retained.
#' @param seed Integer. Seed for the random number generator.
#' Default: NULL, meaning no reproducibility.
#'
#' @return Matrix with one column of SSE or unexplained variance scores for all estimated
#' models. Row names are the value of the cluster parameter for the relevant model.
#' Depends on the choice of \code{return_models}.
#' If \code{TRUE} a list of estimated models is returned.
#' @export
#'
#' @examples
#' # Loading a test dataset into the global environment
#' x <- stackloss
#'
#' # Estimating models with cluster parameter values ranging from 1 to 4
#' model_fits <- mselect_adproclus(data = x, min_nclusters = 1, max_nclusters = 4, seed = 10)
#'
#' # Plot the results as a scree plot to select the appropriate number of clusters
#' plot_scree_adpc(model_fits)
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{adproclus}}}{for the actual ADPROCLUS procedure}
#'   \item{\code{\link{plot_scree_adpc}}}{for plotting the model fits}
#'   \item{\code{\link{select_by_CHull}}}{for automatic model selection via CHull method}
#' }
mselect_adproclus <- function(data, min_nclusters, max_nclusters,
                              return_models = FALSE,
                              unexplvar = TRUE,
                              start_allocation = NULL,
                              nrandomstart = 1, nsemirandomstart = 1,
                              algorithm = "ALS2",
                              save_all_starts = FALSE,
                              seed = NULL) {

        checkmate::assertCount(min_nclusters, positive = TRUE, coerce = TRUE)
        checkmate::assertCount(max_nclusters, positive = TRUE, coerce = TRUE)
        all_models <- list()
        for (i in min_nclusters:max_nclusters) {
                all_models <- append(all_models, list(adproclus(data = data,
                                                                        nclusters = i,
                                                                        start_allocation = start_allocation,
                                                                        nrandomstart = nrandomstart,
                                                                        nsemirandomstart = nsemirandomstart,
                                                                        algorithm = algorithm,
                                                                        save_all_starts = save_all_starts,
                                                                        seed = seed)))
                names(all_models)[length(all_models)] <- paste("model_", i, sep = "")
        }
        class(all_models) <- "adpclist"

        if (return_models) {
                all_models
        } else {
                results <- matrix(min_nclusters:max_nclusters, max_nclusters + 1 - min_nclusters, 1)
                if (unexplvar) {
                        for (i in 1:length(all_models)) {
                                results[i,1] <- 1 - all_models[[i]]$explvar
                        }
                        colnames(results) <- "Unexplained_Variance"

                } else {
                        for (i in 1:length(all_models)) {
                                results[i,1] <- all_models[[i]]$sse
                        }
                        colnames(results) <- "SSE"
                }
                rownames(results) <- min_nclusters:max_nclusters
                results
        }

}

#' Model selection helper for low dimensional ADPROCLUS
#'
#' Performs low dimensional ADPROCLUS for the number of clusters from
#' \code{min_nclusters} to \code{max_nclusters} and the number of components
#' from \code{min_ncomponents} to \code{max_ncomponents}.
#' This replaces the need to manually estimate multiple models to select the best
#' number of clusters and components and returns the results in a format compatible with
#' \code{\link{plot_scree_adpc}} to obtain a scree plot / multiple scree plots.
#' Output is also compatible with \code{\link{select_by_CHull}} to
#' automatically select a suitable number of components for each number of clusters.
#' The compatibility with both functions is only given if
#' \code{return_models = FALSE}.
#'
#' @param data Object-by-variable data matrix of class \code{matrix} or
#'   \code{data.frame}.
#' @param min_nclusters Minimum number of clusters to estimate.
#' @param max_nclusters Maximum number of clusters to estimate.
#' @param min_ncomponents Minimum number of components to estimate.
#' Must be smaller or equal than \code{min_nclusters}.
#' @param max_ncomponents Maximum number of components to estimate.
#' Must be smaller or equal than \code{max_nclusters}.
#' @param return_models Boolean. If \code{FALSE} a matrix of model fit scores is
#' returned, which is compatible with the \code{\link{plot_scree_adpc}} function.
#' If \code{TRUE} the list of actually estimated models is returned.
#' @param unexplvar Boolean. If \code{TRUE} the model fit is specified in terms
#' of unexplained variance. Otherwise it will be specified in terms of
#' Sum of Squared Errors (SSE). This propagates through to the scree plots.
#' @param start_allocation Optional starting cluster membership matrix to be
#' passed to the low dimensional ADPROCLUS procedure.
#' See \code{\link{get_rational}} for more information.
#' @param nrandomstart Number of random starts computed for each model.
#' @param nsemirandomstart Number of semi-random starts computed for each model.
#' @param save_all_starts Logical. If \code{TRUE} and \code{return_models = TRUE},
#' the results of all algorithm starts are returned.
#' By default, only the best solution is retained.
#' @param seed Integer. Seed for the random number generator.
#' Default: NULL, meaning no reproducibility.
#'
#' @return Number of clusters by number of components matrix
#' where the values are SSE or unexplained variance scores for all estimated
#' models. Row names are the value of the cluster parameter for the relevant
#' model. Column names contain the value of the components parameter.
#' Depends on the choice of \code{return_models}.
#' If \code{TRUE} a list of estimated models is returned.
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
#' # Plot the results as a scree plot to select the appropriate number of clusters
#' plot_scree_adpc(model_fits)
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{adproclus_low_dim}}}{for the actual low dimensional ADPROCLUS procedure}
#'   \item{\code{\link{plot_scree_adpc}}}{for plotting the model fits}
#'   \item{\code{\link{select_by_CHull}}}{for automatic model selection via CHull method}
#' }
mselect_adproclus_low_dim <- function(data,
                                      min_nclusters, max_nclusters,
                                      min_ncomponents, max_ncomponents,
                                      return_models = FALSE,
                                      unexplvar = TRUE,
                                      start_allocation = NULL,
                                      nrandomstart = 1, nsemirandomstart = 1,
                                      save_all_starts = FALSE,
                                      seed = NULL) {

        checkmate::assertCount(min_nclusters, positive = TRUE, coerce = TRUE)
        checkmate::assertCount(max_nclusters, positive = TRUE, coerce = TRUE)
        checkmate::assertCount(min_ncomponents, positive = TRUE, coerce = TRUE)
        checkmate::assertCount(max_ncomponents, positive = TRUE, coerce = TRUE)
        if (max_ncomponents > max_nclusters) {
                stop("max_ncomponents must be smaller or equal than max_nclusters")
        }
        if (min_ncomponents > min_nclusters) {
                stop("min_ncomponents must be smaller or equal than min_nclusters")
        }
        all_models <- list()
        for (i in min_nclusters:max_nclusters) {
                all_models_i <- list()
                for (j in min_ncomponents:min(max_ncomponents, i)) {
                        all_models_i <- append(all_models_i,
                                               list(adproclus_low_dim(data = data,
                                                                      nclusters = i,
                                                                      ncomponents = j,
                                                                      start_allocation = start_allocation,
                                                                      nrandomstart = nrandomstart,
                                                                      nsemirandomstart = nsemirandomstart,
                                                                      save_all_starts = save_all_starts,
                                                                      seed = seed)))
                        names(all_models_i)[length(all_models_i)] <- paste("model_", i, "_", j, sep = "")
                }
                all_models <- append(all_models, list(all_models_i))
                names(all_models)[length(all_models)] <- paste("models_", i, sep = "")
        }
        class(all_models) <- "adpclist"

        if (return_models) {
                all_models
        } else {
                ncluster_param <- max_nclusters - min_nclusters + 1
                ncomponent_param <- max_ncomponents - min_ncomponents + 1
                results <- matrix(1:(ncluster_param * ncomponent_param),
                                  ncluster_param,
                                  ncomponent_param)
                if (unexplvar) {
                        for (i in 1:ncluster_param) {
                                results_i <- c()
                                for (j in 1:ncomponent_param) {
                                        if (j > length(all_models[[i]])) {
                                                results_i <- append(results_i, NA)
                                        } else {
                                                results_i <- append(results_i, 1 - all_models[[i]][[j]]$explvar)
                                        }

                                }
                                results[i,] <- results_i
                        }
                        colnames(results) <- sapply(list(min_ncomponents:max_ncomponents),
                                                    function(x) {paste("Unexplained_Variance_", x, "comp", sep = "")})
                } else {
                        for (i in 1:ncluster_param) {
                                results_i <- c()
                                for (j in 1:ncomponent_param) {
                                        if (j > length(all_models[[i]])) {
                                                results_i <- append(results_i, NA)
                                        } else {
                                                results_i <- append(results_i, all_models[[i]][[j]]$sse)
                                        }

                                }
                                results[i,] <- results_i
                        }
                        colnames(results) <- sapply(list(min_ncomponents:max_ncomponents),
                                                    function(x) {paste("SSE_", x, "comp", sep = "")})
                }
                rownames(results) <- min_nclusters:max_nclusters

                results

        }

}


#' Automatic Model Selection for ADPROCLUS with CHull Method
#'
#' For a set of full dimensional ADPROCLUS models (each with different number of clusters),
#' this function finds the "elbow" in the scree plot by using the
#' CHull procedure (Wilderjans, Ceuleman & Meers, 2013) implemented in
#' the \code{\link[multichull]{multichull}} package.
#' For a matrix of low dimensional ADPROCLUS models
#' (each with different number of cluster and components),
#' this function finds the "elbow" in the scree plot for each
#' number of clusters with the CHull methods.
#' That is, it reduces the number of model to choose from to the number of
#' different cluster parameter values by choosing the "elbow" number of
#' components for a given number of clusters. The resulting list can in turn
#' be visualized with \code{\link{plot_scree_adpc_preselected}}.
#' For this procedure to work, the SSE or unexplained variance values must be
#' decreasing in the number of clusters (components). If that is not the case
#' increasing the number of (semi-) random starts can help.
#'
#' This procedure cannot choose the model with
#' the largest or smallest number of clusters (components), i.e. for a set of
#' three models it will always choose the middle one. If for a given number of
#' clusters exactly two models were estimated, this function chooses the model
#' with the lower SSE/unexplained variance.
#'
#' The name of the model fit criterion is propagated from the input matrix based
#' on the first column name. It is either "SSE" or "Unexplained_Variance".
#'
#'
#'
#' @param model_fit Matrix containing SSEs or unexplained variance of all models
#' as in the output of \code{\link{mselect_adproclus}} or \code{\link{mselect_adproclus_low_dim}}.
#' @param percentage_fit Required proportion of increase in fit of a more complex model.
#' @param ... Additional parameters to be passed on to \code{multichull::CHull()} function.
#'
#' @return For full dimensional ADPROCLUS a \code{CHull} object describing the
#' chosen model.
#' For low dimensional ADPROCLUS a matrix containing the list of chosen models
#' and the relevant model parameter, compatible with
#' \code{\link{plot_scree_adpc_preselected}}.
#' @export
#'
#' @references Wilderjans, T. F., Ceulemans, E., & Meers, K. (2012). CHull: A generic convex hull based model
#' selection method. \emph{Behavior Research Methods}, 45, 1-15
#'
#' @examples
#' # Loading a test dataset into the global environment
#' x <- stackloss
#'
#' # Estimating models with cluster parameter values ranging from 1 to 4
#' model_fits <- mselect_adproclus(data = x, min_nclusters = 1, max_nclusters = 4)
#'
#' # Use and visualize CHull method
#' selected_model <- select_by_CHull(model_fits)
#' selected_model
#' plot(selected_model)
#'
#' # Estimating low dimensional models with cluster parameter values
#' # ranging from 1 to 4 and component parameter values also ranging from 1 to 4
#' model_fits <- mselect_adproclus_low_dim(data = x, 1, 4, 1, 4, nsemirandomstart = 10, seed = 1)
#'
#' # Using the CHull method
#' pre_selection <- select_by_CHull(model_fits)
#'
#' # Visualize pre-selected models
#' plot_scree_adpc_preselected(pre_selection)
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{mselect_adproclus}}}{to obtain the \code{model_fit} input from the possible ADPROCLUS models}
#'   \item{\code{\link{mselect_adproclus_low_dim}}}{to obtain the \code{model_fit} input from the possible low dimensional ADPROCLUS models}
#'   \item{\code{\link{plot_scree_adpc}}}{for plotting the model fits}
#' }
select_by_CHull <- function(model_fit, percentage_fit = 0.0001, ...) {
        if (substring(colnames(model_fit)[[1]], 1, 1) == "S") {
                fit_var <- "SSE"
        } else {
                fit_var <- "Unexplained_Variance"
        }
        if (ncol(model_fit) == 1) {
                multichull::CHull(cbind(strtoi(rownames(model_fit)), model_fit[, 1]), bound = "lower", PercentageFit = percentage_fit, ...)
        } else {
                clusters <- strtoi(rownames(model_fit))
                components <- c()
                model_fit_result <- c()
                for (i in 1:nrow(model_fit)) {
                        model_data <- cbind(readr::parse_number(colnames(model_fit)), model_fit[i,])
                        model_data <- model_data[rowSums(is.na(model_data)) != 1, , drop = FALSE]
                        if (nrow(model_data) == 1) {
                                components <- append(components, model_data[1, 1])
                                model_fit_result <- append(model_fit_result, model_data[1, 2])
                        } else if (nrow(model_data) == 2) {
                                components <- append(components, model_data[which.min(model_data[, 2]), 1])
                                model_fit_result <- append(model_fit_result, min(model_data[, 2]))
                        } else {
                                model <- multichull::CHull(data = model_data, bound = "lower", PercentageFit = percentage_fit, ...)
                                components <- append(components, model$Solution$complexity)
                                model_fit_result <- append(model_fit_result, model$Solution$fit)
                        }
                }
                results <- cbind(clusters, components, model_fit_result)
                rownames(results) <- paste("Cl:", clusters, "; Comp:", components, sep = "")
                colnames(results)[3] <- fit_var
                results
        }
}
