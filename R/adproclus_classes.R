#' An S4 class to represent an ADPROCLUS solution.
#'
#' @slot clusters A length-one numeric vector denoting the number of clusters in the solution
#' @slot Model A model matrix M of the same size as the input data matrix.
#' @slot Membs The n*k membership matrix of the clustering model
#' @slot Profs The k*j profile matrix of the clustering model
#' @slot sse A numeric denoting the residual sum of squares of the clustering model
#' @slot totvar A numeric denoting the total sum of squares of the input data matrix
#' @slot explvar A numeric denoting the proportion of variance in data that is accounted for by the clustering model.
#' @slot alg_iter A numeric denoing the number of iterations of the ALS algorithm
#' @slot timer A numeric denoting the elapsed time in seconds to obtain the clustering solution
#' @slot initialStart A list containing information about the stat used to initialize the algorithm

setClass("adproclus", representation(
        clusters = "numeric",
        Model = "matrix",
        Membs = "matrix",
        Profs = "matrix",
        sse = "numeric",
        totvar = "numeric",
        explvar = "numeric",
        alg_iter = "numeric",
        timer = "numeric",
        initialStart = "list",
        Runs = "list"))

setMethod ("show", "adproclus", function(object){
        cat ("ADPROCLUS overlapping clustering solution with", object@clusters, "clusters. \n")
        cat (" Type             :", class(object), '\n')
        cat ("This is a test and not fully implemented yet")
        cat ("\n")
})
