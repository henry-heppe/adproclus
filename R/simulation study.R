#simulation study for low dimensional adproclus



data_shape_possib <- c(4, 1, 0.25)
number_clusters_possib <- c(4, 5, 6)
cluster_overlap_possib <- c(0.25, 0.5, 0.75)
number_dimensions_possib <- c(1, 2, 4)
noise_level_possib <- c(0, 0.1, 0.2, 0.3)
noise_correlation_possib <- c(0, 0.3)


all_data <- list()
for (i in 1:length(data_shape_possib)) {
        for (j in 1:length(number_clusters_possib)) {
                for (l in 1:length(cluster_overlap_possib)) {
                        for (m in 1:length(number_dimensions_possib)) {
                                for (n in 1:length(noise_level_possib)) {
                                        for (q in 1:length(noise_correlation_possib)) {
                                                A <- gen_A(data_shape = i, number_clusters = j, cluster_overlap = l)
                                                P <- gen_P(data_shape = i, number_clusters = j, number_dimensions = m)
                                                M <- A %*% P
                                                E <- gen_errors(data_shape = i, M = M, noise_level = n, noise_correlation = q)

                                                all_data <- append(all_data, M + E)
                                        }
                                }
                        }
                }
        }
}


#generate A depending on parameters
gen_A <- function(data_shape, number_clusters, cluster_overlap) {
        J <- sqrt(4096/data_shape)
        I <- 4096/J
        k <- number_clusters
        #all possible 2^k permutations where the rows of A are chosen from
        possibA_rows = gtools::permutations(2, k, v = c(0, 1), repeats.allowed = TRUE)
        probsA <- rep(0, 2^k)
        probsA[1] <- 0.05 #probability for all zeros (belonging to no cluster at all)
        probsA[which(rowSums(possibA_rows) == 1)] <- (1-cluster_overlap-0.05)/k #probabilities for rows that belong to one cluster only
        probsA[probsA == 0] <- cluster_overlap / (2^k - (k+1))  #probs for all remaining possible rows

        multinom_results <- rmultinom(n = I, size = 1, prob = probsA) #draw from multinomial distribution
        indeces_A <- which(multinom_results == 1, arr.ind = TRUE)[,1] #recover which rows to take from all possible rows
        A <- possibA_rows[indeces_A,]

}

#generate P depending on parameters
gen_P <- function(data_shape, number_clusters, number_dimensions) {
        J <- sqrt(4096/data_shape)
        I <- 4096/J
        k <- number_clusters
        s <- number_dimensions

        B_temp <- matrix(rnorm(J*s, mean = 0, sd = 1),J,s)
        B <- svd(B_temp)$u #issue: columns (should be rows?) of B are left singular vectors of B_temp
        repeat {
                C <- matrix(runif(k*s, -20, 20), k, s) #generate C until P is of rank s
                P <- C %*% t(B)
                P_decomp <- svd(P)$d
                P_decomp[abs(P_decomp) < 1e-10] <- 0
                rank_P <- length(which(P_decomp != 0))
                if (rank_P == s) {
                        break
                }
        }
        return(P)
}

gen_errors <- function(data_shape, M, noise_level, noise_correlation) {
        #size of E: I*J = 4096, I/J = data_shape
        J <- sqrt(4096/data_shape)
        I <- 4096/J
        E <- matrix(rep(0,I*J),I,J)
        sigma_E <- matrix(rep(noise_correlation, J^2), J, J) #covariance matrix of errors
        diag(sigma_E) <- noise_level #issue: what are the correct values for the diagonal?

        set.seed(1)
        E <- MASS::mvrnorm(n = I, mu = rep(0,J), Sigma = sigma_E)
        return(list(errors = E, noise_level = norm(E - mean(E), "f")^2/norm(M+E - mean(M+E), "f")^2))

}
