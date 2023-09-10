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
summary.adpc <- function(object, ...) {
        #print(model$sse)
        cat("hello")


}

print.summary.adpc <- function() {

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
#'
#' @export
#'
#' @examples
#' x <- ADPROCLUS::CGdata[1:100,]
#' model <- adproclus(x, 3)
#' print(model)
print.adpc <- function(x, title = "ADPROCLUS solution", ...) {
        if(!is.null(x$C)) {
                cat("Low Dimensional", title)
                print("Setup:")
                cat("number of clusters:", ncol(x$A))
                cat("number of dimensions:", ncol(x$C))
                cat("data format: ", nrow(x$model), "x", ncol(x$model))
                cat("number of runs:", length(x$runs))
                print("Results:")
                cat("explained variance:", x$explvar, "time:", x$timer, "iterations:", x$iterations)
                cat("First few rows of A:", x$A[1:5])
                cat("First few profiles in terms of components:", x$C[1:10])
                cat("Components by variables matrix B:", x$B)
        } else {
                print(title)
                print("Setup:")
                cat("number of clusters:", ncol(x$A))
                cat("data format: ", nrow(x$model), "x", ncol(x$model))
                cat("number of runs:", length(x$runs))
                print("Results:")
                cat("explained variance:", x$explvar, "time:", x$timer, "iterations:", x$iterations)
                cat("First few rows of A:", x$A[1:5])
                cat("First few profiles:", x$P[1:5])
        }
}

plot.adpc <- function() {

}
