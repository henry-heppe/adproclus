

#' Title
#'
#' @param data
#' @param min_nclusters
#' @param max_nclusters
#' @param start_allocation
#' @param nrandomstart
#' @param nsemirandomstart
#' @param algorithm
#' @param save_all_starts
#' @param seed
#' @param return_models
#' @param unexplvar
#'
#' @return
#' @export
#'
#' @examples
mselect_adproclus <- function(data, min_nclusters, max_nclusters,
                              return_models = FALSE,
                              unexplvar = TRUE,
                              start_allocation = NULL,
                              nrandomstart = 1, nsemirandomstart = 1,
                              algorithm = "ALS1",
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

#' Title
#'
#' @param data
#' @param min_nclusters
#' @param max_nclusters
#' @param start_allocation
#' @param nrandomstart
#' @param nsemirandomstart
#' @param save_all_starts
#' @param seed
#' @param min_ncomponents
#' @param max_ncomponents
#' @param return_models
#' @param unexplvar
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#' only works when the unexplained variance/SSE is decreasing in the number of components
#' if for a certain number of clusters exactly 2 models where estimated, this procedure will always choose the one with the better model fit
#'
#' @param model_fit
#' @param PercentageFit
#'
#' @return
#' @export
#'
#' @examples
select_by_CHull <- function(model_fit, percentage_fit = 0.0001) {
        if (substring(colnames(model_fit)[[1]], 1, 1) == "S") {
                fit_var <- "SSE"
        } else {
                fit_var <- "Unexplained_Variance"
        }
        if (ncol(model_fit) == 1) {
                multichull::CHull(cbind(strtoi(rownames(model_fit)), model_fit[, 1]), bound = "lower", PercentageFit = percentage_fit)
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
                                model <- multichull::CHull(data = model_data, bound = "lower", PercentageFit = percentage_fit)
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
