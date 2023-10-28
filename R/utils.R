#' Repeat matrix in two dimensions
#'
#' @param a Matrix to repeat.
#' @param n Dimension 1.
#' @param m Dimension 2.
#'
#' @return Larger Matrix.
#'
#' @noRd
.repmat <- function(a, n, m) {
        kronecker(matrix(1, n, m), a)
}

#' Squared Loss Function
#'
#' @param w Input matrix.
#'
#' @return Sum of squares.
#'
#' @noRd
.lossL2 <- function(w) {
        sum(w^2)
}


