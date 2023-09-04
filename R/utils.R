#utility function
.repmat <- function(a, n, m) {
        kronecker(matrix(1, n, m), a)
}

.lossL2 <- function(w) {
        m <- sum(w^2)
        return(m)
}

