#' @include adproclus_classes.R
NULL

.repmat <- function(a, n, m) {
        kronecker(matrix(1, n, m), a)
}

.lossL2 <- function(w) {
        m <- sum(w^2)
        return(m)
}

.updateA_lf1 <- function(A, PossibA, data) {

        B <- A
        npos <- nrow(PossibA)
        n <- nrow(A)
        loss <- matrix(0, 1, npos)

        for (row in 1:n) {

                for (r in 1:npos) {

                        Apos <- B
                        Apos[row, ] <- PossibA[r, ]
                        if (any(colSums(Apos) == 0) == FALSE) {
                                Gpos <- NMFN::mpinv(Apos) %*% data
                                respos <- data - (Apos %*% Gpos)
                                loss[r] <- .lossL2(respos)
                        } else {
                                loss[r] <- NA
                        }

                }


                minloss <- min(loss[, !is.na(loss)])
                ind <- which.min(loss)
                B[row, ] <- PossibA[ind, ]
        }

        if (any(colSums(B) == 0) == TRUE) {
                zerocols <- 0
        } else {
                zerocols <- 1
        }
        result <- list(B, minloss, zerocols)
        return(result)
}

.updateA_lf2 <- function(n, P, replX, PossibA) {

        nA <- nrow(PossibA)

        posrec <- PossibA %*% P

        loss <- matrix(as.numeric(rowSums((replX - posrec)^2)),
                nrow = nA/n, ncol = n)

        index <- apply(loss, 2, which.min)
        A <- PossibA[index, ]

        return(A)

}

.adproclus_lf1 <- function(x, A) {

        t <- proc.time()

        data <- x
        n <- nrow(x)
        totvar <- norm(data - mean(data), "f")^2
        k <- ncol(A)

        PossibA <- gtools::permutations(2, k, v = c(0, 1), repeats.allowed = TRUE)
        PossibA <- apply(PossibA, 2, rev)
        if (k > 1) {
                PossibA <- t(apply(PossibA, 1, rev))
        }

        G <- NMFN::mpinv(A) %*% data
        res <- data - (A %*% G)
        f <- .lossL2(res)
        fold <- f + 1
        iter <- 1

        while (fold > f) {

                fold = f
                Aold = A

                update <- .updateA_lf1(Aold, PossibA, data)
                A <- update[[1]]
                f <- update[[2]]

                iter <- iter + 1
        }

        runs <- iter - 1
        Membs <- Aold
        Profs <- NMFN::mpinv(Aold) %*% data
        sse <- fold
        explvar <- 1 - (fold/totvar)

        timeruns <- (proc.time() - t)[1]

        model <- Membs %*% Profs
        result <- list(Model = model, Membs = Membs, Profs = Profs,
                sse = sum((model - x)^2), totvar = totvar,
                explvar = explvar, alg_iter = runs, timer = as.numeric(timeruns))

        return(result)
}

.adproclus_lf2 <- function(x, A, P) {

        t <- proc.time()

        data <- x
        n <- nrow(x)
        totvar <- norm(data - mean(data), "f")^2
        k <- ncol(A)
        npos <- 2^k

        G <- NMFN::mpinv(A) %*% data
        res <- data - (A %*% G)
        f <- .lossL2(res)
        fold <- f + 1
        iter <- 1

        PossibA <- gtools::permutations(2, k, v = c(0, 1), repeats.allowed = TRUE)
        PossibA <- apply(PossibA, 2, rev)
        if (k > 1) {
                PossibA <- t(apply(PossibA, 1, rev))
        }
        PossibA <- .repmat(PossibA, n, 1)

        replX <- data.frame()
        for (i in 1:n) {
                reps <- matrix(.repmat(data[i, ], npos, 1),
                        ncol = ncol(data), nrow = npos, byrow = TRUE)
                replX <- rbind(replX, reps)
        }
        replX <- as.matrix(replX)

        while (fold > f) {

                fold = f
                Aold = A
                Pold = P

                A <- .updateA_lf2(n, Pold, replX, PossibA)
                P <- NMFN::mpinv(A) %*% data
                f <- .lossL2(x - (A %*% P))

                iter <- iter + 1
        }

        runs <- iter - 1
        Membs <- Aold
        Profs <- Pold
        sse <- fold
        explvar <- 1 - (fold/totvar)

        timeruns <- (proc.time() - t)[1]

        model <- Membs %*% Profs
        result <- list(Model = model, Membs = Membs, Profs = Profs,
                sse = sum((model - x)^2), totvar = totvar,
                explvar = explvar, alg_iter = runs, timer = as.numeric(timeruns))

        return(result)

}


.ldadproclus <- function (x, A, s) {
        t <- proc.time()

        data <- x
        n <- nrow(x)
        totvar <- norm(data - mean(data), "f")^2
        k <- ncol(A)
        npos <- 2^k

        P <- .ldfindP(x, A, s)

        f <- .lossL2(x - (A %*% P))
        fold <- f + 1
        iter <- 1

        PossibA <- gtools::permutations(2, k, v = c(0, 1), repeats.allowed = TRUE)
        PossibA <- apply(PossibA, 2, rev)
        if (k > 1) {
                PossibA <- t(apply(PossibA, 1, rev))
        }
        PossibA <- .repmat(PossibA, n, 1)

        replX <- data.frame()
        for (i in 1:n) {
                reps <- matrix(.repmat(data[i, ], npos, 1),
                               ncol = ncol(data), nrow = npos, byrow = TRUE)
                replX <- rbind(replX, reps)
        }
        replX <- as.matrix(replX)

        while (fold > f) {

                fold = f
                Aold = A
                Pold = P

                A <- .updateA_lf2(n, Pold, replX, PossibA)
                P <- .ldfindP(data, A, s)
                f <- .lossL2(x - (A %*% P))

                iter <- iter + 1
        }

        runs <- iter - 1
        Membs <- A   #or Aold instead?
        Profs <- P   #or Pold instead?
        sse <- fold
        explvar <- 1 - (fold/totvar)

        decomp <- svd(P)

        ldBase <- decomp$v[,1:s] #B (correct to neglect columns and sv after s?)
        ldProfs <- decomp$u[,1:s] %*% diag(decomp$d[1:s]) #C


        timeruns <- (proc.time() - t)[1]

        model <- Membs %*% Profs
        result <- list(Model = model, Membs = Membs, Profs = Profs, low_dim_profs = ldProfs, low_dim_base = ldBase,
                       sse = sum((model - x)^2), totvar = totvar,
                       explvar = explvar, alg_iter = runs, timer = as.numeric(timeruns), svd = svd(Profs))

        return(result)




}

.ldfindP <- function (X, A, s) {
        Z <- A %*% matlib::inv(t(A) %*% A) %*% t(A)
        U <- Z %*% X %*% t(X) %*% Z
        Q <- eigen(U, symmetric = TRUE)[["vectors"]][,1:s]

        P <- matlib::inv((t(A) %*% A)) %*% t(A) %*% Q %*% t(Q) %*% X

        return(P)
}









.printoutput <- function (k, time, nstart) {

        print ("Clustering completed.")
        print (paste ("Additive Profile Clustering with",
                k, "overlapping clusters."))
        print (paste ("Total processing time:",
                round (time, digits = 1), "seconds."))
        print (paste ("Algorithm starts:", nstart[1], "random and",
                gtools::na.replace (nstart[2], 0), "rational."))

}

#' Generate initial random start
#'
#' Generate an initial random start for the Additive Profile Clustering
#' algorithm (see \code{\link{adproclus}}).
#'
#' \code{getRandom} generates an  random initial binary membership matrix
#' \strong{A} by drawing entries from a Bernoulli Distribution with \eqn{\pi =
#' 0.5}. A corresponding initial profile matrix \strong{P} is subsequently
#' estimated conditionally upon A (for details, see Depril et al., 2008, and
#' Wilderjans et al., 2010).
#'
#' For programming simplicity, this function provides the option to pass a
#' matrix of initial cluster centers to the \code{centers} argument. In this
#' case, the matrix is dismissed, the number of clusters \emph{k} is set to
#' \code{nrow(centers)} and a new profile matrix is returned, based on a
#' randomly generated membership matrix. For generating an initial start based
#' on a specific set of initial cluster centers, see \code{\link{getRational}}.
#'
#' \strong{Warning:} This function does \emph{not} obtain an ADPRCOLUS model. To
#' perform aditive profile clustering, see \code{\link{adproclus}}.
#'
#' @param data Object-by-variable data matrix of class \code{matrix} or \code{data.frame}.
#' @param centers either the number of clusters \emph{k}, or a matrix of initial
#'   (distinct) cluster centres.
#'
#' @return \code{getRandom} returns a list with the following components:
#'   \describe{ \item{\code{type}}{A character string denoting the type of start
#'   ('Random Start')} \item{\code{A}}{A randomly generated initial Membership
#'   matrix} \item{\code{P}}{The corresponding initial Profile matrix} }
#'
#' @export
#'
#' @references Wilderjans, T. F., Ceulemans, E., Van Mechelen, I., & Depril, D.
#' (2010). ADPROCLUS: a graphical user interface for fitting additive profile
#' clustering models to object by variable data matrices. \emph{Behavior
#' Research Methods, 43}(1), 56-65.
#'
#' Depril, D., Van Mechelen, I., & Mirkin, B.
#' (2008). Algorithms for additive clustering of rectangular data tables.
#' \emph{Computational Statistics and Data Analysis, 52,} 4923-4938.
#'
#' @examples
#' getRandom(ADPROCLUS::CGdata, 3)
#'
#' @seealso \code{\link{getRational}} for generating rational starts and
#'   \code{\link{adproclus}} for details about membership and profile matrices.
getRandom <- function(data, centers) {

        data <- as.matrix(data)

        if (length(centers) == 1L) {
                k <- centers
        } else {
                k <- nrow(centers)
        }
        n <- nrow(data)

        A <- (matrix(stats::runif(n * k), n, k) < 0.5) * 1
        while (any(colSums(A) == 0) || qr(A)$rank < k) {
                A <- (matrix(stats::runif(n * k), n, k) < 0.5) * 1
        }
        P <- NMFN::mpinv(A) %*% data

        return(list(type = "Random Start", A = A, P = P))
}

#' Generate initial rational start
#'
#' Generate an initial rational start for the Additive Profile Clustering
#' algorithm (see \code{\link{adproclus}}).
#'
#' If \code{centers} is a number of clusters \emph{k}, an initial profile matrix
#' \strong{P} is generated by drawing \emph{k} randomly chosen, distinct, rows
#' from \code{data}. Alternatively, you can pass a user selected matrix of size
#' \emph{k} * \code{ncol(data)} with initial cluster profiles to the
#' \code{centers} argument. In both cases, a corresponding membership matrix
#' \strong{A} is subsequently estimated conditionally upon \strong{P} (for
#' details, see Depril et al., 2008, and Wilderjans et al., 2010).
#'
#' \strong{Warning:} This function does \emph{not} obtain an ADPRCOLUS model. To
#' perform aditive profile clustering, see \code{\link{adproclus}}.
#'
#' @param data Object-by-variable data matrix of class \code{matrix} or
#'   \code{data.frame}.
#' @param centers either the number of clusters \emph{k}, or a matrix of initial
#'   (distinct) cluster centres.
#'
#' @return \code{getRational} returns a list with the following components:
#'   \describe{
#'   \item{\code{type}}{A character string denoting the type of start
#'   ('Rational Start')}
#'   \item{\code{A}}{An initial Membership matrix}
#'   \item{\code{P}}{The corresponding initial Profile matrix} }
#'
#' @export
#'
#' @references Wilderjans, T. F., Ceulemans, E., Van Mechelen, I., & Depril, D.
#'   (2010). ADPROCLUS: a graphical user interface for fitting additive profile
#'   clustering models to object by variable data matrices. \emph{Behavior
#'   Research Methods, 43}(1), 56-65.
#'
#'   Depril, D., Van Mechelen, I., & Mirkin, B. (2008). Algorithms for additive
#'   clustering of rectangular data tables. \emph{Computational Statistics and
#'   Data Analysis, 52,} 4923-4938.
#'
#' @examples
#' getRational(ADPROCLUS::CGdata, 3)
#'
#' @seealso \code{\link{getRandom}} for generating random starts and
#'   \code{\link{adproclus}} for details about membership and profile matrices.
getRational <- function(data, centers) {

        data <- as.matrix(data)
        n <- nrow(data)

        if (length(centers) == 1L) {
                k <- centers
                Permutation <- sample(n)
                P <- data[Permutation[1:k], ]
        } else {
                k <- nrow(centers)
                P <- centers
        }

        npos <- 2^k
        PossibA <- gtools::permutations(2, k, v = c(0, 1), repeats.allowed = TRUE)
        PossibA <- apply(PossibA, 2, rev)
        if (k > 1) {
                PossibA <- t(apply(PossibA, 1, rev))
        }
        PossibA <- .repmat(PossibA, n, 1)

        replX <- data.frame()
        for (i in 1:n) {
                reps <- matrix(.repmat(data[i, ], npos, 1),
                        ncol = ncol(data), nrow = npos, byrow = TRUE)
                replX <- rbind(replX, reps)
        }
        replX <- as.matrix(replX)

        A <- as.matrix(.updateA_lf2(n, P, replX, PossibA))
        result <- list(type = "Rational Start", A = A,
                P = P)
        return(result)
}

#' Additive profile clustering
#'
#' Perform ADditive PROfile CLUStering (ADPRCOLUS) on object by variable data.
#'
#' In this function, Mirkin's (1987, 1990) Aditive Profile Clustering
#' (ADPROCLUS) method is used to obtain an unrestricted overlapping clustering
#' model of the object by variable data provided by \code{data}.
#'
#' The ADPROCLUS model approximates an \emph{I} \eqn{x} \emph{J} object by
#' variable data matrix \strong{X} by an \emph{I} \eqn{x} \emph{J} model matrix
#' \strong{M} that can be decomposed into an \emph{I} \eqn{x} \emph{K} binary
#' cluster membership matrix \strong{A} and a \emph{K} \eqn{x} \emph{J}
#' real-valued cluster profile matrix \strong{P}, with \emph{K} indicating the
#' number of overlapping clusters. In particular, the aim of an ADPROCLUS
#' analysis is therefore, given a number of clusters \emph{k}, to estimate a
#' model matrix \strong{M} = \strong{AP} that reconstructs data matrix
#' \strong{X} as close as possible in a least squares sense (i.e. sum of squared
#' residuals). For a detailed illustration of the ADPROCLUS model and associated
#' loss function, see Wilderjans et al., 2011.
#'
#' The alternating least squares algorithms ("\code{ALS1}" and "\code{ALS2}")
#' that can be used for minimization of the loss function were proposed by
#' Depril et al. (2008). In "\code{ALS2}", starting from an initial random or
#' rational estimate of \strong{A} (see \code{\link{getRandom}} and
#' \code{\link{getRational}}), \strong{A} and \strong{P} are alternately
#' re-estimated conditionally upon each other until convergence. The
#' "\code{ALS1}" algorithm differs from the one previous one in that each row in
#' \strong{A} is updated independently and that the conditionally optimal
#' \strong{P} is recalculated after each row update, instead of the end of the
#' matrix. For a discussion and comparison of the different algorithms, see
#' Depril et al., 2008.
#'
#' \strong{Warning:} Computation time increases exponentially with increasing
#' number of clusters, \emph{k}! We recommend to determine the computation time
#' of a single start for each specific dataset and \emph{k} before employing a
#' multistart procedure.
#'
#' @param data object-by-variable data matrix of class \code{matrix} or
#'   \code{data.frame}.
#' @param centers either the number of clusters \emph{k}, or a matrix of initial
#'   (distinct) cluster centres. If a number \emph{k}, a random set of \emph{k}
#'   rows in \code{data} is chosen as initial centres.
#' @param nstart if \code{centers} is a number, a vector of length 1 or 2
#'   denoting the number of random and rational starts to be performed.
#' @param algorithm character string "\code{ALS1}" (default) or "\code{ALS2}",
#'   denoting the type of alternating least squares algorithm.
#' @param SaveAllStarts logical. If \code{TRUE}, the results of all algorithm
#'   starts are returned. By default, only the best solution is retained.
#'
#' @return By default, \code{adproclus} returns a list with the following
#'   components: (If \code{SaveAllStarts} is \code{TRUE}, a list is returned for
#'   each start of the algorithm) \describe{ \item{\code{Model}}{matrix. The
#'   obtained overlapping clustering model \strong{M} of the same size as
#'   \code{data}.} \item{\code{Membs}}{matrix. The membership matrix \strong{A}
#'   of the clustering model.} \item{\code{Profs}}{matrix. The profile matrix
#'   \strong{P} of the clustering model.} \item{\code{sse}}{numeric. The
#'   residual sum of sqares of the clustering model, which is minimised by the
#'   ALS algorithm.} \item{\code{totvar}}{numeric. The total sum of squares
#'   of\code{data}.} \item{\code{explvar}}{numeric. The proportion of variance
#'   in \code{data} that is accounted for by the clustering model.}
#'   \item{\code{alg_iter}}{numeric. The number of iterations of the algorithm.}
#'   \item{\code{timer}}{numeric. The amount of time (in seconds) the algorithm
#'   ran for.} \item{\code{initialStart}}{list. A list containing initial
#'   membership and profile matrices, as well as the type of start that was used
#'   to obtain the clustering solution. (as returned by \code{\link{getRandom}}
#'   or \code{\link{getRational}})}}
#'
#' @export
#'
#' @references Wilderjans, T. F., Ceulemans, E., Van Mechelen, I., & Depril, D.
#'   (2010). ADPROCLUS: a graphical user interface for fitting additive profile
#'   clustering models to object by variable data matrices. \emph{Behavior
#'   Research Methods, 43}(1), 56-65.
#'
#'   Depril, D., Van Mechelen, I., & Mirkin, B. (2008). Algorithms for additive
#'   clustering of rectangular data tables. \emph{Computational Statistics and
#'   Data Analysis, 52,} 4923-4938.
#'
#'   Mirkin, B. G. (1987). The method of principal clusters. \emph{Automation
#'   and Remote Control}, 10:131-143.
#'
#'   Mirkin, B. G. (1990). A sequential fitting procedure for linear data
#'   analysis models. \emph{Journal of Classification}, 7(2):167-195.
#'
#' @examples
#' # Loading a test dataset into the global environment
#' x <- ADPROCLUS::CGdata
#'
#' # Quick clustering with K = 3 clusters
#' clust <- adproclus(x, 3)
#'
#' # Clustering with K = 4 clusters,
#' # using the ALS2 algorithm,
#' # with 5 random and 5 rational starts
#' clust <- adproclus(x, 4, c(5,5), "ALS2")
#'
#' # Saving the results of all starts
#' clust <- adproclus(x, 3, c(2,2), SaveAllStarts = TRUE)
#'
#' # Clustering using a user-defined rational start
#' start <- getRational(x,3)
#' clust <- adproclus(x, start$P)
#'
#' @seealso \code{\link{getRandom}} and \code{\link{getRational}} for generating
#'   random and rational starts for ADORCLUS.
adproclus <- function(data, centers, nstart = 1L,
        algorithm = "ALS1", SaveAllStarts = FALSE) {

    t <- proc.time()
    results <- list()

    data <- as.matrix(data)
    n <- as.integer(nrow(data))
    p <- as.integer(ncol(data))
    if (is.na(n) || is.na(p))
        stop("'invalid data matrix")
    if (missing(centers))
        stop("'centers' must be a number or a matrix")

    alg <- switch(match.arg(algorithm, c("ALS1", "ALS2")),
        ALS1 = 1, ALS2 = 2)

    if (length(centers) != 1L) {
        if (ncol(data) != ncol(centers))
            stop("number of columns in 'centers' must match number of columns in 'data'")
        k <- nrow(centers)
        centers <- as.matrix(centers)
        nstart <- c(0,nstart)
        if (n < k)
            stop("number of clusters cannot exceed number of objects in 'data'")
        if (sum(nstart) != 1L) {
            nstart <- c(0,1)
            warning("'centers' is an initial profile matrix.
                    Number of starts has been set to one.")
        }
        initialStart <- getRational(data,
            centers)
        if (alg == 1) {
            results <- .adproclus_lf1(data, initialStart$A)
            results$initialStart <- initialStart

        }
        if (alg == 2) {
            results <- .adproclus_lf2(data, initialStart$A,
                initialStart$P)
            results$initialStart <- initialStart
        }
    } else {
        k <- centers
        if (n < k)
            stop("number of clusters cannot exceed number of objects in 'data'")

        if (length(nstart) > 2) {
            stop("'nstart' must be a vector of length 1 or 2")
        }
        if (sum(nstart) > 50) {
            warning("Number of starts is larger than 50. Computation might take a while")
        }

        BestSol <- list(sse = Inf)

        if (alg == 1) {
            for (i in 1:nstart[1]) {
                initialStart <- getRandom(data,
                  centers)
                res <- .adproclus_lf1(data, initialStart$A)
                res$initialStart <- initialStart
                if (res$sse < BestSol$sse) {
                  BestSol <- res
                }
                if (SaveAllStarts == TRUE) {
                  results <- append(results, list(res))
                }
            }
            remove(i)
            if (!is.na(nstart[2])) {
                for (i in 1:nstart[2]) {
                  initialStart <- getRational(data,
                    centers)
                  res <- .adproclus_lf1(data, initialStart$A)
                  res$initialStart <- initialStart
                  if (res$sse < BestSol$sse) {
                    BestSol <- res
                  }
                  if (SaveAllStarts == TRUE) {
                    results <- append(results, list(res))
                  }
                }
                remove(i)
            }
        }

        if (alg == 2) {
            for (i in 1:nstart[1]) {
                initialStart <- getRandom(data,
                  centers)
                res <- .adproclus_lf2(data, initialStart$A,
                  initialStart$P)
                res$initialStart <- initialStart
                if (res$sse < BestSol$sse) {
                  BestSol <- res
                }
                if (SaveAllStarts == TRUE) {
                  results <- append(results, list(res))
                }
                remove(i)
            }
            if (!is.na(nstart[2])) {
                for (i in 1:nstart[2]) {
                  initialStart <- getRational(data,
                    centers)
                  res <- .adproclus_lf2(data, initialStart$A,
                    initialStart$P)
                  res$initialStart <- initialStart
                  if (res$sse < BestSol$sse) {
                    BestSol <- res
                  }
                  if (SaveAllStarts == TRUE) {
                    results <- append(results, list(res))
                  }
                }
                remove(i)
            }
        }

        if (SaveAllStarts == TRUE) {
            results <- list(BestSol = BestSol, Runs = results)
            names(results$Runs) <- as.character(c(1:length(results$Runs)))
        } else {
            results <- BestSol
        }
    }
    time <- (proc.time() - t)[1]
    .printoutput(k,time, nstart)
    return(results)
}



#' Low dimensional ADPROCLUS
#'
#' @param data the data
#' @param nclusters the number of clusters
#' @param ncomponents the number of dimensions (the number of components)
#' @param start_allocation a possible start allocation
#' @param nrandomstart the number of random starts
#' @param randomstart type of random start
#' @param SaveAllStarts option to save all starts
#'
#' @return a list of relevant information about the resulting model
#' @export
#'
#' @examples some example
ldadproclus <- function(data, nclusters, ncomponents, start_allocation = NULL, nrandomstart = 1,
                      randomstart = c("random", "semi-random"), SaveAllStarts = FALSE) {

        t <- proc.time()
        results <- list()

        data <- as.matrix(data)
        n <- as.integer(nrow(data))
        p <- as.integer(ncol(data))
        if (is.na(n) || is.na(p))
                stop("invalid data matrix")
        if (!is.integer(as.integer(nclusters)))
                stop("'nclusters' must be a number")

        if (ncomponents >= min(n,nclusters)) {
                stop("'ncomponents' must be smaller than min(number of observations, number of clusters)")
        }
        if (!is.null(start_allocation)) {
                start_allocation <- as.matrix(start_allocation)
                m <- as.integer(nrow(start_allocation))
                q <- as.integer(ncol(start_allocation))
                if (is.na(m) || is.na(q) || !all(start_allocation %in% c(0,1)) || m != n || q != nclusters) {
                        stop("invalid start allocation matrix")
                }

        }
        randomstart <- match.arg(randomstart)












        #old code (for adproclus, not LDadrproclus)

        if (length(centers) != 1L) {
                if (ncol(data) != ncol(centers))
                        stop("number of columns in 'centers' must match number of columns in 'data'")
                k <- nclusters
                if (n < k)
                        stop("number of clusters cannot exceed number of objects in 'data'")
                if (sum(nstart) != 1L) {
                        nstart <- c(0,1)
                        warning("'centers' is an initial profile matrix.
                    Number of starts has been set to one.")
                }
                initialStart <- getRational(data,
                                            centers)
                if (alg == 1) {
                        results <- .adproclus_lf1(data, initialStart$A)
                        results$initialStart <- initialStart

                }
                if (alg == 2) {
                        results <- .adproclus_lf2(data, initialStart$A,
                                                  initialStart$P)
                        results$initialStart <- initialStart
                }
        } else {
                k <- centers
                if (n < k)
                        stop("number of clusters cannot exceed number of objects in 'data'")

                if (length(nstart) > 2) {
                        stop("'nstart' must be a vector of length 1 or 2")
                }
                if (sum(nstart) > 50) {
                        warning("Number of starts is larger than 50. Computation might take a while")
                }

                BestSol <- list(sse = Inf)

                if (alg == 1) {
                        for (i in 1:nstart[1]) {
                                initialStart <- getRandom(data,
                                                          centers)
                                res <- .adproclus_lf1(data, initialStart$A)
                                res$initialStart <- initialStart
                                if (res$sse < BestSol$sse) {
                                        BestSol <- res
                                }
                                if (SaveAllStarts == TRUE) {
                                        results <- append(results, list(res))
                                }
                        }
                        remove(i)
                        if (!is.na(nstart[2])) {
                                for (i in 1:nstart[2]) {
                                        initialStart <- getRational(data,
                                                                    centers)
                                        res <- .adproclus_lf1(data, initialStart$A)
                                        res$initialStart <- initialStart
                                        if (res$sse < BestSol$sse) {
                                                BestSol <- res
                                        }
                                        if (SaveAllStarts == TRUE) {
                                                results <- append(results, list(res))
                                        }
                                }
                                remove(i)
                        }
                }

                if (alg == 2) {
                        for (i in 1:nstart[1]) {
                                initialStart <- getRandom(data,
                                                          centers)
                                res <- .adproclus_lf2(data, initialStart$A,
                                                      initialStart$P)
                                res$initialStart <- initialStart
                                if (res$sse < BestSol$sse) {
                                        BestSol <- res
                                }
                                if (SaveAllStarts == TRUE) {
                                        results <- append(results, list(res))
                                }
                                remove(i)
                        }
                        if (!is.na(nstart[2])) {
                                for (i in 1:nstart[2]) {
                                        initialStart <- getRational(data,
                                                                    centers)
                                        res <- .adproclus_lf2(data, initialStart$A,
                                                              initialStart$P)
                                        res$initialStart <- initialStart
                                        if (res$sse < BestSol$sse) {
                                                BestSol <- res
                                        }
                                        if (SaveAllStarts == TRUE) {
                                                results <- append(results, list(res))
                                        }
                                }
                                remove(i)
                        }
                }

                if (SaveAllStarts == TRUE) {
                        results <- list(BestSol = BestSol, Runs = results)
                        names(results$Runs) <- as.character(c(1:length(results$Runs)))
                } else {
                        results <- BestSol
                }
        }
        time <- (proc.time() - t)[1]
        .printoutput(k,time, nstart)
        return(results)
}

basic_test <- function() {
        set.seed(1)
        k <- 6 #number of clusters
        s <- 3 #number of dimensions

        a <- matrix(rbinom(400*k, 1, 0.5), 400, k)
        a <- data.frame(a)
        a_ext <- cbind(a, rsum = rowSums(a))
        not_done = TRUE
        iterr <- 0
        while(not_done) {
                iterr <- iterr + 1
                a_filt <- a_ext[a_ext$rsum == 0 | a_ext$rsum == k,]
                a_ext[a_ext$rsum == 0 | a_ext$rsum == k,] <- data.frame(matrix(rbinom(k*nrow(a_filt),1,0.5), nrow(a_filt), k))
                if (nrow(a_ext[a_ext$rsum == 0 | a_ext$rsum == k,])) {
                     not_done <- FALSE
                }
                a <- a_ext[,1:k]
                a_ext <- cbind(a, rsum = rowSums(a))
        }
        print("A created")
        .ldadproclus(as.matrix(CGdata), as.matrix(a), s)
}

