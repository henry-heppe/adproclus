#simulation study for low dimensional adproclus







all_data <- list()

gen_A <- function(data_shape, cluster_overlap, number_clusters) {
        J <- sqrt(4096/data_shape)
        I <- 4096/J
        k <- number_clusters
        possibA_rows = gtools::permutations(2, k, v = c(0, 1), repeats.allowed = TRUE)
        probsA <- rep(0, 2^k)
        probsA[1] <- 0.05
        probsA[where(rowSums(possibA_rows) == 1)] <- (1-cluster_overlap)/k
        probsA[probsA == 0] <- cluster_overlap / (2^k - (k+1))

        A <- possibA_rows[rmultinom(n = I, size = 1, prob = probsA),]


}

gen_P <- function(data_shape) {

}

gen_errors <- function(data_shape, M, noise_level, noise_correlation) {
        #size of E: I*J = 4096, I/J = data_shape
        J <- sqrt(4096/data_shape)
        I <- 4096/J
        E <- matrix(rep(0,I*J),I,J)
        sigma_E <- matrix(rep(noise_correlation, J^2), J, J)
        diag(sigma_E) <- noise_level

        set.seed(1)
        for (i in 1:I) {
                E[i,] <- NULL
        }
        E <- MASS::mvrnorm(n = I, mu = rep(0,J), Sigma = sigma_E)
        return(list(E, noise_level = var(E)/var(E+M)))

}
