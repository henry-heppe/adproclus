

#' Title
#'
#' @param data
#' @param min_nclusters
#' @param max_nclusters
#' @param screeplot
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
                              screeplot = TRUE,
                              start_allocation = NULL,
                              nrandomstart = 1, nsemirandomstart = 1,
                              algorithm = "ALS1",
                              save_all_starts = FALSE,
                              seed = NULL) {

        checkmate::assertCount(min_nclusters, positive = TRUE, coerce = TRUE)
        checkmate::assertCount(max_nclusters, positive = TRUE, coerce = TRUE)
        checkmate::assertFlag(screeplot)
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
        all_models
}
