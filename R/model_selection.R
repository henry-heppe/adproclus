

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
#' @param algorithm
#' @param save_all_starts
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
mselect_adproclus_low_dim <- function(data,
                                      min_nclusters, max_nclusters,
                                      min_ncomponents, max_ncomponents,
                                      return_models = FALSE,
                              start_allocation = NULL,
                              nrandomstart = 1, nsemirandomstart = 1,
                              save_all_starts = FALSE,
                              seed = NULL) {

        checkmate::assertCount(min_nclusters, positive = TRUE, coerce = TRUE)
        checkmate::assertCount(max_nclusters, positive = TRUE, coerce = TRUE)
        checkmate::assertCount(min_ncomponents, positive = TRUE, coerce = TRUE)
        checkmate::assertCount(max_ncomponents, positive = TRUE, coerce = TRUE)
        all_models <- list()
        for (i in min_nclusters:max_nclusters) {
                all_models_i <- list()
                for (j in min(min_ncomponents, i):min(max_ncomponents, i)) {
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

        }

}
