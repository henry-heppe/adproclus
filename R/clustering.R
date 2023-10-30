#' @include adproclus_classes.R
NULL

#' Additive profile clustering
#'
#' Perform additive profile clustering (ADPROCLUS) on object by variable data.
#'
#' In this function, Mirkin's (1987, 1990) Additive Profile Clustering
#' (ADPROCLUS) method is used to obtain an unrestricted overlapping clustering
#' model of the object by variable data provided by \code{data}.
#'
#' The ADPROCLUS model approximates an \eqn{I \times J} object by
#' variable data matrix \eqn{\boldsymbol{X}} by an \eqn{I \times J} model matrix
#' \eqn{\boldsymbol{M}} that can be decomposed into an \eqn{I \times K} binary
#' cluster membership matrix \eqn{\boldsymbol{A}} and a \eqn{K \times J}
#' real-valued cluster profile matrix \eqn{\boldsymbol{P}}, with \eqn{K}
#' indicating the number of overlapping clusters.
#' In particular, the aim of an ADPROCLUS analysis is therefore,
#' given a number of clusters \eqn{K}, to estimate a
#' model matrix \deqn{M = AP} which reconstructs the data matrix
#' \eqn{\boldsymbol{X}} as close as possible in a least squares sense
#' (i.e. sum of squared residuals). For a detailed illustration of the
#' ADPROCLUS model and associated loss function, see Wilderjans et al., 2011.
#'
#' The alternating least squares algorithms ("\code{ALS1}" and "\code{ALS2}")
#' that can be used for minimization of the loss function were proposed by
#' Depril et al. (2008). In "\code{ALS2}", starting from an initial random or
#' rational estimate of \eqn{\boldsymbol{A}} (see \code{\link{get_random}} and
#' \code{\link{get_semirandom}}), \eqn{\boldsymbol{A}} and \eqn{\boldsymbol{P}}
#' are alternately re-estimated conditionally upon each other until convergence.
#' The "\code{ALS1}" algorithm differs from the one previous one in that each
#' row in \eqn{\boldsymbol{A}} is updated independently and that the
#' conditionally optimal \eqn{\boldsymbol{P}} is recalculated after each row
#' update, instead of the end of the matrix. For a discussion and comparison of
#' the different algorithms, see Depril et al., 2008.
#'
#' \strong{Warning:} Computation time increases exponentially with increasing
#' number of clusters, \eqn{K}. We recommend to determine the computation time
#' of a single start for each specific dataset and \eqn{K} before employing a
#' multistart procedure.
#'
#' @param data Object-by-variable data matrix of class \code{matrix} or
#'   \code{data.frame}.
#' @param nclusters Number of clusters to be used. Must be a positive integer.
#' @param start_allocation Optional matrix of binary values as starting
#' allocation for first run. Default is \code{NULL}.
#' @param nrandomstart Number of random starts (see \code{\link{get_random}}).
#' Can be zero. Increase for better results, though longer computation time.
#' Some research finds 500 starts to be a useful reference.
#' @param nsemirandomstart Number of semi-random starts
#' (see \code{\link{get_semirandom}})). Can be zero. Increase for better
#' results, though longer computation time.
#' Some research finds 500 starts to be a useful reference.
#' @param algorithm character string "\code{ALS1}" (default) or "\code{ALS2}",
#'   denoting the type of alternating least squares algorithm.
#' @param save_all_starts logical. If \code{TRUE}, the results of all algorithm
#'   starts are returned. By default, only the best solution is retained.
#' @param seed Integer. Seed for the random number generator.
#' Default: NULL, meaning no reproducibility
#'
#' @return \code{adproclus} returns a list with the following
#'   components, which describe the best model (from the multiple starts):
#'   \describe{
#'   \item{\code{model}}{matrix. The obtained overlapping clustering model
#'   \strong{M} of the same size as \code{data}.}
#'   \item{\code{A}}{matrix. The membership matrix \strong{A} of the clustering
#'   model. Clusters are sorted by size.}
#'   \item{\code{P}}{matrix. The profile matrix
#'   \strong{P} of the clustering model.}
#'   \item{\code{sse}}{numeric. The
#'   residual sum of squares of the clustering model, which is minimized by the
#'   ALS algorithm.}
#'   \item{\code{totvar}}{numeric. The total sum of squares
#'   of \code{data}.}
#'   \item{\code{explvar}}{numeric. The proportion of variance
#'   in \code{data} that is accounted for by the clustering model.}
#'   \item{\code{iterations}}{numeric. The number of iterations of the
#'   algorithm.}
#'   \item{\code{timer}}{numeric. The amount of time (in seconds) the complete
#'   algorithm ran for.}
#'   \item{\code{timer_one_run}}{numeric. The amount of time (in seconds) the
#'   relevant single start ran for.}
#'   \item{\code{initial_start}}{list. Containing the initial
#'   membership matrix, as well as the type of start that was used
#'   to obtain the clustering solution. (as returned by \code{\link{get_random}}
#'   or \code{\link{get_semirandom}})}
#'   \item{\code{runs}}{list. Each element represents one model obtained from
#'   one of the multiple starts.
#'   Each element contains all of the above information.}
#'   \item{\code{parameters}}{list. Containing the parameters used for the
#'   model.}}
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
#' x <- stackloss
#'
#' # Quick clustering with K = 2 clusters
#' clust <- adproclus(data = x, nclusters = 2)
#'
#' # Clustering with K = 3 clusters,
#' # using the ALS2 algorithm,
#' # with 2 random and 2 semi-random starts
#' clust <- adproclus(x, 3,
#'   nrandomstart = 2, nsemirandomstart = 2, algorithm = "ALS2"
#' )
#'
#' # Saving the results of all starts
#' clust <- adproclus(x, 3,
#'   nrandomstart = 2, nsemirandomstart = 2, save_all_starts = TRUE
#' )
#'
#' # Clustering using a user-defined rational start profile matrix
#' # (here the first 4 rows of the data)
#' start <- get_rational(x, x[1:4, ])$A
#' clust <- adproclus(x, 4, start_allocation = start)
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{adproclus_low_dim}}}{for low dimensional ADPROCLUS}
#'   \item{\code{\link{get_random}}}{for generating random starts}
#'   \item{\code{\link{get_semirandom}}}{for generating semi-random starts}
#'   \item{\code{\link{get_rational}}}{for generating rational starts}
#' }
adproclus <- function(data, nclusters,
                      start_allocation = NULL,
                      nrandomstart = 3, nsemirandomstart = 3,
                      algorithm = "ALS1",
                      save_all_starts = FALSE,
                      seed = NULL) {
  t <- proc.time()
  results <- list()

  data <- as.matrix(data)
  n <- as.integer(nrow(data))

  checkmate::assertCount(nclusters, positive = TRUE, coerce = TRUE)
  checkmate::assertCount(nrandomstart, coerce = TRUE)
  checkmate::assertCount(nsemirandomstart, coerce = TRUE)
  checkmate::assertInt(seed, null.ok = TRUE, coerce = TRUE)
  checkmate::assertFlag(save_all_starts)
  checkmate::assertMatrix(data, mode = "numeric", any.missing = FALSE)
  checkmate::assertChoice(algorithm, c("ALS1", "ALS2"))


  if (nrandomstart + nsemirandomstart > 50) {
    warning("Number of starts is larger than 50.
            Computation might take a while")
  }

  if (n < nclusters) {
    stop("Number of clusters must be less or equal the number of
         objects in 'data'.")
  }
  if (is.null(start_allocation) && nrandomstart == 0 && nsemirandomstart == 0) {
    stop("Must have either start allocation matrix or a non-zero number of
         (semi-) random starts")
  }
  if (nclusters > ncol(data)) {
    stop("Number of clusters must be less or equal the number of variables.")
  }

  alg <- switch(match.arg(algorithm, c("ALS1", "ALS2")),
    ALS1 = 1,
    ALS2 = 2
  )

  BestSol <- list(sse = Inf)
  results <- list()
  if (!is.null(start_allocation)) {
    start_allocation <- as.matrix(start_allocation)
    checkmate::assertMatrix(start_allocation)
    m <- as.integer(nrow(start_allocation))
    q <- as.integer(ncol(start_allocation))
    if (!all(start_allocation %in% c(0, 1)) || m != n || q != nclusters) {
      stop("invalid start allocation matrix: either non-binary matrix,
                         or number of rows or columns do not match data")
    }

    if (alg == 1) {
      res <- .adproclus_lf1(data, start_allocation)
    } else {
      res <- .adproclus_lf2(data, start_allocation)
    }

    res$initial_start <- list(
      type = "Rational Start",
      initialA = start_allocation
    )
    results <- append(results, list(res))
    BestSol <- res
  }
  if (alg == 1) {
    if (nrandomstart > 0) {
      seed_local <- seed
      for (i in 1:nrandomstart) {
        if (!is.null(seed)) {
          seed_local <- seed_local + 1
        }
        initial_start <- get_random(data,
          nclusters,
          seed = seed_local
        )
        initial_start$type <- paste("random_start_no_", i)
        res <- .adproclus_lf1(data, initial_start$A)
        res$initial_start <- initial_start
        if (res$sse < BestSol$sse) {
          BestSol <- res
        }
        results <- append(results, list(res))
      }
      remove(i)
    }
    if (nsemirandomstart > 0) {
      seed_local <- seed
      for (i in 1:nsemirandomstart) {
        if (!is.null(seed)) {
          seed_local <- seed_local + 1
        }
        initial_start <- get_semirandom(data,
          nclusters,
          seed = seed_local
        )
        initial_start$type <- paste("semi_random_start_no_", i)
        res <- .adproclus_lf1(data, initial_start$A)
        res$initial_start <- initial_start
        if (res$sse < BestSol$sse) {
          BestSol <- res
        }
        results <- append(results, list(res))
      }
      remove(i)
    }
  } else {
    if (nrandomstart > 0) {
      seed_local <- seed
      for (i in 1:nrandomstart) {
        if (!is.null(seed)) {
          seed_local <- seed_local + 1
        }
        initial_start <- get_random(data,
          nclusters,
          seed = seed_local
        )
        initial_start$type <- paste("random_start_no_", i)
        res <- .adproclus_lf2(data, initial_start$A)
        res$initial_start <- initial_start
        if (res$sse < BestSol$sse) {
          BestSol <- res
        }
        results <- append(results, list(res))
      }
      remove(i)
    }
    if (nsemirandomstart > 0) {
      seed_local <- seed
      for (i in 1:nsemirandomstart) {
        if (!is.null(seed)) {
          seed_local <- seed_local + 1
        }
        initial_start <- get_semirandom(data,
          nclusters,
          seed = seed_local
        )
        initial_start$type <- paste("semi_random_start_no_", i)
        res <- .adproclus_lf2(data, initial_start$A)
        res$initial_start <- initial_start
        if (res$sse < BestSol$sse) {
          BestSol <- res
        }
        results <- append(results, list(res))
      }
      remove(i)
    }
  }

  if (save_all_starts == TRUE) {
    BestSol$runs <- results
    results <- BestSol
    names(results$runs) <- as.character(seq_along(results$runs))
  } else {
    results <- BestSol
  }

  parameters <- list(
    nclusters = nclusters,
    nrandomstart = nrandomstart,
    nsemirandomstart = nsemirandomstart,
    start_allocation = start_allocation,
    seed = seed
  )
  results$parameters <- parameters

  time <- (proc.time() - t)[1]
  results$timer <- time
  results <- .adjust_row_col_names(results, data)
  results
}



#' Low dimensional ADPROCLUS
#'
#' Perform \strong{low dimensional} additive profile clustering (ADPROCLUS) on
#' object by variable data. Use case: data to cluster consists of a large set of
#' variables, where it can be useful to interpret the cluster profiles in terms
#' of a smaller set of components that represent the original variables well.
#'
#' In this function, an extension by Depril et al. (2012) of
#' Mirkins (1987, 1990) additive profile clustering method is used to obtain a
#' low dimensional overlapping clustering model of the object by variable data
#' provided by \code{data}.
#' More precisely, the low dimensional ADPROCLUS model approximates an
#' \eqn{I \times J} object by variable data matrix \eqn{\boldsymbol{X}} by an
#' \eqn{I \times J} model matrix \eqn{\boldsymbol{M}}. For \eqn{K} overlapping
#' clusters, \eqn{\boldsymbol{M}} can be decomposed into an \eqn{I \times K}
#' binary cluster membership matrix \eqn{\boldsymbol{A}} and a \eqn{K \times J}
#' real-valued cluster profile matrix \eqn{\boldsymbol{P}} s.t. \deqn{M = AP.}
#' With the simultaneous dimension reduction, \eqn{\boldsymbol{P}} is restricted
#' to be of reduced rank \eqn{S < \min(K,J)}, such that it can be decomposed
#' into \deqn{P = CB',} with \eqn{\boldsymbol{C}} a \eqn{K \times S} matrix and
#' \eqn{\boldsymbol{B}} a \eqn{J \times S} matrix. Now, a row in
#' \eqn{\boldsymbol{C}} represents the profile values associated with the
#' respective cluster in terms of the \eqn{S} components, while
#' the entries of \eqn{\boldsymbol{B}} can be used to interpret the components
#' in terms of the complete set of variables. In particular, the aim of an
#' ADPROCLUS analysis is therefore, given a number of clusters \eqn{K} and a
#' number of dimensions \eqn{S}, to estimate a model matrix \eqn{\boldsymbol{M}}
#' that reconstructs data matrix
#' \eqn{\boldsymbol{X}} as close as possible in a least squares sense and
#' simultaneously reduce the dimensions of the data.
#' For a detailed illustration of the low dimensional ADPROCLUS model and
#' associated loss function, see Depril et al. (2012).
#'
#' \strong{Warning:} Computation time increases exponentially with increasing
#' number of clusters, \eqn{K}. We recommend to determine the computation time
#' of a single start for each specific dataset and \eqn{K} before employing a
#' multistart procedure.
#'
#' @param data Object-by-variable data matrix of class \code{matrix} or
#'   \code{data.frame}.
#' @param nclusters Number of clusters to be used. Must be a positive integer.
#' @param ncomponents Number of components (dimensions) to which the profiles
#' should be restricted. Must be a positive integer.
#' @param start_allocation Optional matrix of binary values as starting
#' allocation for first run. Default is \code{NULL}.
#' @param nrandomstart Number of random starts (see \code{\link{get_random}}).
#' Can be zero. Increase for better results, though longer computation time.
#' Some research finds 500 starts to be a useful reference.
#' @param nsemirandomstart Number of semi-random starts
#' (see \code{\link{get_semirandom}})). Can be zero.
#' Increase for better results, though longer computation time.
#' Some research finds 500 starts to be a useful reference.
#' @param save_all_starts logical. If \code{TRUE}, the results of all algorithm
#'   starts are returned. By default, only the best solution is retained.
#' @param seed Integer. Seed for the random number generator.
#' Default: NULL, meaning no reproducibility
#'
#' @return \code{adproclus_low_dim} returns a list with the following
#'   components, which describe the best model (from the multiple starts):
#'   \describe{
#'   \item{\code{model}}{matrix. The obtained overlapping clustering model
#'   \eqn{\boldsymbol{M}} of the same size as \code{data}.}
#'   \item{\code{model_lowdim}}{matrix. The obtained low dimensional clustering
#'   model \eqn{\boldsymbol{AC}} of size \eqn{I \times S}}
#'   \item{\code{A}}{matrix. The membership matrix \eqn{\boldsymbol{A}} of the
#'   clustering model. Clusters are sorted by size.}
#'   \item{\code{P}}{matrix. The profile matrix \eqn{\boldsymbol{P}} of the
#'   clustering model.}
#'   \item{\code{c}}{matrix. The profile values in terms of the low dimensional
#'   components.}
#'   \item{\code{B}}{matrix. Base vectors connecting low dimensional components
#'   with original variables.
#'   Variables-by-components matrix. Warning: for computing
#'   \eqn{\boldsymbol{P}} use \eqn{\boldsymbol{B'}}}
#'   \item{\code{sse}}{numeric. The
#'   residual sum of squares of the clustering model, which is minimized by the
#'   ALS algorithm.}
#'   \item{\code{totvar}}{numeric. The total sum of squares
#'   of \code{data}.}
#'   \item{\code{explvar}}{numeric. The proportion of variance
#'   in \code{data} that is accounted for by the clustering model.}
#'   \item{\code{iterations}}{numeric. The number of iterations of the
#'   algorithm.}
#'   \item{\code{timer}}{numeric. The amount of time (in seconds) the complete
#'   algorithm ran for.}
#'   \item{\code{timer_one_run}}{numeric. The amount of time (in seconds) the
#'   relevant single start ran for.}
#'   \item{\code{initial_start}}{list. A list containing the initial
#'   membership matrix, as well as the type of start that was used
#'   to obtain the clustering solution. (as returned by \code{\link{get_random}}
#'   or \code{\link{get_semirandom}})}
#'   \item{\code{runs}}{list. Each element represents one model obtained
#'   from one of the multiple starts.
#'   Each element contains all of the above information.}
#'   \item{\code{parameters}}{list. Containing the parameters used for the
#'   model.}}
#'
#' @export
#'
#' @references Depril, D., Van Mechelen, I., & Wilderjans, T. F.
#'   (2012). Lowdimensional additive overlapping clustering.
#'   \emph{Journal of classification, 29,} 297-320.
#'
#' @examples
#' # Loading a test dataset into the global environment
#' x <- stackloss
#'
#' # Low dimensional clustering with K = 3 clusters
#' # where the resulting profiles can be characterized in S = 1 dimensions
#' clust <- adproclus_low_dim(x, 3, ncomponents = 1)
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{adproclus}}}{for full dimensional ADPROCLUS}
#'   \item{\code{\link{get_random}}}{for generating random starts}
#'   \item{\code{\link{get_semirandom}}}{for generating semi-random starts}
#'   \item{\code{\link{get_rational}}}{for generating rational starts}
#' }
adproclus_low_dim <- function(data, nclusters, ncomponents,
                              start_allocation = NULL,
                              nrandomstart = 3, nsemirandomstart = 3,
                              save_all_starts = FALSE,
                              seed = NULL) {
  t <- proc.time()
  results <- list()

  data <- as.matrix(data)
  n <- as.integer(nrow(data))

  checkmate::assertCount(nclusters, positive = TRUE, coerce = TRUE)
  checkmate::assertCount(ncomponents, positive = TRUE, coerce = TRUE)
  checkmate::assertCount(nrandomstart, coerce = TRUE)
  checkmate::assertCount(nsemirandomstart, coerce = TRUE)
  checkmate::assertInt(seed, null.ok = TRUE, coerce = TRUE)
  checkmate::assertFlag(save_all_starts)
  checkmate::assertMatrix(data, mode = "numeric", any.missing = FALSE)


  if (ncomponents > min(n, nclusters)) {
    stop("'ncomponents' must be smaller than min(number of observations,
         number of clusters).")
  }
  if (is.null(start_allocation) && nrandomstart == 0 && nsemirandomstart == 0) {
    stop("Must have either start allocation matrix or a non-zero number of
         (semi-) random starts")
  }
  if (n < nclusters) {
    stop("Number of clusters must be less or equal the number of objects
         in 'data'.")
  }
  if (nrandomstart + nsemirandomstart > 50) {
    warning("Number of starts is larger than 50.
            Computation might take a while")
  }

  best_sol <- list(sse = Inf)
  if (!is.null(start_allocation)) {
    start_allocation <- as.matrix(start_allocation)
    checkmate::assertMatrix(start_allocation)

    m <- as.integer(nrow(start_allocation))
    q <- as.integer(ncol(start_allocation))
    if (!all(start_allocation %in% c(0, 1)) || m != n || q != nclusters) {
      stop("invalid start allocation matrix: either non-binary matrix,
                             or number of rows or columns do not match data")
    }

    model_new <- .ldadproclus(data, start_allocation, ncomponents)
    model_new$initial_start <- list(
      type = "rational_start_model",
      A = start_allocation
    )
    results <- append(results, list(model_new))
    best_sol <- results[[1]]
  }
  if (nrandomstart > 0) {
    seed_local <- seed
    for (i in 1:nrandomstart) {
      if (!is.null(seed_local)) {
        seed_local <- seed_local + 1
      }
      random_start <- get_random(data, nclusters, seed = seed_local)$A
      model_new <- .ldadproclus(data, random_start, ncomponents)

      model_new$initial_start <- list(
        type = paste("random_start_no_", i),
        A = random_start
      )

      results <- append(results, list(model_new))
      if (model_new$sse < best_sol$sse) {
        best_sol <- model_new
      }
    }
  }

  if (nsemirandomstart > 0) {
    seed_local <- seed
    for (j in 1:nsemirandomstart) {
      if (!is.null(seed_local)) {
        seed_local <- seed_local + 1
      }
      semi_random_start <- get_semirandom(data, nclusters, seed = seed_local)$A
      model_new <- .ldadproclus(data, semi_random_start, ncomponents)

      model_new$initial_start <- list(
        type = paste("semi_random_start_no_", j),
        A = semi_random_start
      )

      results <- append(results, list(model_new))
      if (model_new$sse < best_sol$sse) {
        best_sol <- model_new
      }
    }
  }

  if (save_all_starts == TRUE) {
    best_sol$runs <- results
    results <- best_sol
    names(results$runs) <- as.character(seq_along(results$runs))
  } else {
    results <- best_sol
  }

  parameters <- list(
    nclusters = nclusters,
    ncomponents = ncomponents,
    nrandomstart = nrandomstart,
    nsemirandomstart = nsemirandomstart,
    start_allocation = start_allocation,
    seed = seed
  )
  results$parameters <- parameters
  time <- (proc.time() - t)[1]
  results$timer <- time
  results <- .adjust_row_col_names_LD(results, data)
  results
}

#' Update A for full dim ADPROCLUS ALS 1
#'
#' @param A Last A.
#' @param PossibA All possibilities for A.
#' @param data Input data.
#'
#' @return Updated A.
#'
#' @noRd
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
    B[row, ] <- PossibA[ind, , drop = FALSE]
  }

  if (any(colSums(B) == 0) == TRUE) {
    zerocols <- 0
  } else {
    zerocols <- 1
  }
  result <- list(B, minloss, zerocols)
  result
}

#' Update A for full dim ADPROCLUS ALS 1
#'
#' @param n Number of observations.
#' @param P Profile matrix.
#' @param replX Original data, where each row is multiplicated 2^k times,
#' for k clusters.
#' @param PossibA All possibilities for A.
#'
#' @return Updated A.
#'
#' @noRd
.updateA_lf2 <- function(n, P, replX, PossibA) {
  nA <- nrow(PossibA)

  posrec <- PossibA %*% P

  loss <- matrix(as.numeric(rowSums((replX - posrec)^2)),
    nrow = nA / n, ncol = n
  )

  index <- apply(loss, 2, which.min)
  A <- as.matrix(PossibA[index, , drop = FALSE])
  A
}

#' ADPROCLUS internal for ALS 1
#'
#' @param x Input Data.
#' @param A Starting cluster membership matrix.
#'
#' @return Object of class \code{"adpc"}.
#'
#' @noRd
.adproclus_lf1 <- function(x, A) {
  t <- proc.time()

  data <- x
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
    fold <- f
    Aold <- A

    update <- .updateA_lf1(Aold, PossibA, data)
    A <- update[[1]]
    f <- update[[2]]

    iter <- iter + 1
  }

  runs <- iter - 1
  Membs <- Aold
  Profs <- NMFN::mpinv(Aold) %*% data
  sse <- fold
  explvar <- 1 - (fold / totvar)

  # sort A columns and P rows by decreasing cluster size
  order <- order(colSums(Membs), decreasing = TRUE)
  Membs <- Membs[, order, drop = FALSE]
  Profs <- Profs[order, , drop = FALSE]

  timeruns <- (proc.time() - t)[1]

  model <- Membs %*% Profs
  result <- list(
    model = model, A = Membs, P = Profs,
    sse = sum((model - x)^2), totvar = totvar, explvar = explvar,
    iterations = runs, timer_one_run = as.numeric(timeruns),
    initial_start = NULL
  )
  class(result) <- "adpc"
  result
}

#' ADPROCLUS internal for ALS 2
#'
#' @param x Input Data.
#' @param A Starting cluster membership matrix.
#'
#' @return Object of class \code{"adpc"}.
#'
#' @noRd
.adproclus_lf2 <- function(x, A) {
  t <- proc.time()

  data <- x
  n <- nrow(x)
  totvar <- norm(data - mean(data), "f")^2
  k <- ncol(A)
  npos <- 2^k

  P <- NMFN::mpinv(A) %*% data
  res <- data - (A %*% P)
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
    reps <- .repmat(data[i, , drop = FALSE], npos, 1)
    replX <- rbind(replX, reps)
  }
  replX <- as.matrix(replX)

  while (fold > f) {
    fold <- f
    Pold <- P

    A <- .updateA_lf2(n, Pold, replX, PossibA)
    P <- NMFN::mpinv(A) %*% data
    f <- .lossL2(x - (A %*% P))

    iter <- iter + 1
  }

  runs <- iter - 1
  Membs <- A
  Profs <- P
  sse <- f
  explvar <- 1 - (f / totvar)

  timeruns <- (proc.time() - t)[1]

  # sort A columns and P rows by decreasing cluster size
  order <- order(colSums(Membs), decreasing = TRUE)
  Membs <- Membs[, order, drop = FALSE]
  Profs <- Profs[order, , drop = FALSE]


  model <- Membs %*% Profs
  result <- list(
    model = model, A = Membs, P = Profs,
    sse = sum((model - x)^2), totvar = totvar, explvar = explvar,
    iterations = runs, timer_one_run = as.numeric(timeruns),
    initial_start = NULL
  )
  class(result) <- "adpc"

  result
}


#' Low dimensional ADPROCLUS internal
#'
#' @param x Input Data.
#' @param A Starting cluster membership matrix.
#' @param s Number of dimensions.
#'
#' @return Object of class \code{"adpc"}.
#'
#' @noRd
.ldadproclus <- function(x, A, s) {
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

  # possible restriction against zero-rows
  # PossibA <- PossibA[which(rowSums(PossibA) > 0),]
  # npos <- npos - 1

  PossibA <- apply(PossibA, 2, rev)
  if (k > 1) {
    PossibA <- t(apply(PossibA, 1, rev))
  }

  PossibA <- .repmat(PossibA, n, 1)

  replX <- data.frame()
  for (i in 1:n) {
    reps <- .repmat(data[i, , drop = FALSE], npos, 1)
    replX <- rbind(replX, reps)
  }
  replX <- as.matrix(replX)

  while (fold > f) {
    fold <- f
    Pold <- P

    A <- .updateA_lf2(n, Pold, replX, PossibA)
    P <- .ldfindP(data, A, s)
    f <- .lossL2(x - (A %*% P))

    iter <- iter + 1
  }

  runs <- iter - 1
  Membs <- A
  Profs <- P
  sse <- f
  explvar <- 1 - (f / totvar)

  # sort A columns and P rows by decreasing cluster size
  order <- order(colSums(Membs), decreasing = TRUE)
  Membs <- Membs[, order, drop = FALSE]
  Profs <- Profs[order, , drop = FALSE]


  decomp <- svd(Profs)

  # B
  ldBase <- as.matrix(decomp$v[, 1:s, drop = FALSE])
  U <- as.matrix(decomp$u[, 1:s, drop = FALSE])
  if (s == 1) {
    D <- diag(as.matrix(decomp$d[1:s]))
  } else {
    D <- diag(decomp$d[1:s])
  }

  # C = U*D
  ldProfs <- as.matrix(U %*% D)


  timeruns <- (proc.time() - t)[1]

  model <- Membs %*% Profs
  model_lowdim <- Membs %*% ldProfs
  result <- list(
    model = model, model_lowdim = model_lowdim,
    A = Membs, P = Profs, C = ldProfs, B = ldBase,
    sse = sum((model - x)^2), totvar = totvar, explvar = explvar,
    iterations = runs, timer_one_run = as.numeric(timeruns),
    initial_start = NULL
  )
  class(result) <- "adpc"

  result
}

#' Update P for low dimensional ADPROCLUS
#'
#' @param X Input data.
#' @param A Last cluster membership matrix.
#' @param s Number of dimensions.
#'
#' @return New P.
#'
#' @noRd
.ldfindP <- function(X, A, s) {
  Z <- A %*% NMFN::mpinv(A)
  U <- Z %*% X %*% t(X) %*% Z
  Q <- eigen(U, symmetric = TRUE)[["vectors"]][, 1:s, drop = FALSE]

  P <- NMFN::mpinv(A) %*% Q %*% t(Q) %*% X
  P
}


#' Propagate input variable names to output
#'
#' @param object ADPROCLUS solution.
#' @param data Input data.
#'
#' @return ADPROCLUS solution with updated variable names.
#'
#' @noRd
.adjust_row_col_names <- function(object, data) {
  results <- object
  # row and column names
  colnames(results$model) <- colnames(data, do.NULL = FALSE, prefix = "V")
  colnames(results$A) <- colnames(results$A, do.NULL = FALSE, prefix = "Cl")
  rownames(results$P) <- rownames(results$P, do.NULL = FALSE, prefix = "Cl")

  if (!is.null(results$runs)) {
    for (i in seq_along(results$runs)) {
      colnames(results$runs[[i]]$model) <- colnames(data,
        do.NULL = FALSE,
        prefix = "V"
      )
      colnames(results$runs[[i]]$A) <- colnames(results$A,
        do.NULL = FALSE,
        prefix = "Cl"
      )
      rownames(results$runs[[i]]$P) <- rownames(results$P,
        do.NULL = FALSE,
        prefix = "Cl"
      )
    }
  }
  results
}

#' Propagate input variable names to output for low dimensional model
#'
#' @param object Low dimensional ADPROCLUS solution.
#' @param data Input data.
#'
#' @return Low dimensional ADPROCLUS solution with updated variable names.
#'
#' @noRd
.adjust_row_col_names_LD <- function(object, data) {
  results <- object

  # row and column names
  colnames(results$model) <- colnames(data, do.NULL = FALSE, prefix = "V")
  colnames(results$model_lowdim) <- colnames(results$model_lowdim,
    do.NULL = FALSE, prefix = "Comp"
  )
  colnames(results$A) <- colnames(results$A, do.NULL = FALSE, prefix = "Cl")
  rownames(results$P) <- rownames(results$P, do.NULL = FALSE, prefix = "Cl")
  colnames(results$C) <- colnames(results$C, do.NULL = FALSE, prefix = "Comp")
  rownames(results$C) <- rownames(results$C, do.NULL = FALSE, prefix = "Cl")
  colnames(results$B) <- colnames(results$B, do.NULL = FALSE, prefix = "Comp")
  rownames(results$B) <- rownames(results$B, do.NULL = FALSE, prefix = "V")

  if (!is.null(results$runs)) {
    for (i in seq_along(results$runs)) {
      colnames(results$runs[[i]]$model) <- colnames(data,
        do.NULL = FALSE,
        prefix = "V"
      )
      colnames(results$runs[[i]]$model_lowdim) <- colnames(results$model_lowdim,
        do.NULL = FALSE,
        prefix = "Comp"
      )
      colnames(results$runs[[i]]$A) <- colnames(results$A,
        do.NULL = FALSE,
        prefix = "Cl"
      )
      rownames(results$runs[[i]]$P) <- rownames(results$P,
        do.NULL = FALSE,
        prefix = "Cl"
      )
      colnames(results$runs[[i]]$C) <- colnames(results$C,
        do.NULL = FALSE,
        prefix = "Comp"
      )
      rownames(results$runs[[i]]$C) <- rownames(results$C,
        do.NULL = FALSE,
        prefix = "Cl"
      )
      colnames(results$runs[[i]]$B) <- colnames(results$B,
        do.NULL = FALSE,
        prefix = "Comp"
      )
      rownames(results$runs[[i]]$B) <- rownames(results$B,
        do.NULL = FALSE,
        prefix = "V"
      )
    }
  }

  results
}
