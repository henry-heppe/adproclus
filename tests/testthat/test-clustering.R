test_that("ADPROCLUS base case normal input", {
  x <- stackloss
  expect_no_error(adproclus(x, nclusters = 2, seed = 1))
  expect_no_error(adproclus(x,
    nclusters = 3,
    nrandomstart = 1, nsemirandomstart = 1, algorithm = "ALS2",
    seed = 1
  ))
  expect_no_error(adproclus(x,
    nclusters = 3,
    nrandomstart = 1, nsemirandomstart = 1, algorithm = "ALS1",
    seed = 1
  ))
  expect_no_error(adproclus(x,
    nclusters = 3,
    nrandomstart = 2, nsemirandomstart = 2, save_all_starts = TRUE,
    seed = 1
  ))
})

test_that("ADPROCLUS with start_allocation ", {
  x <- stackloss
  start <- get_rational(x, x[1:4, ])$A
  expect_no_error(adproclus(x, nclusters = 4, start_allocation = start, seed = 1))
  expect_no_error(adproclus(x,
    nclusters = 4,
    nrandomstart = 0, nsemirandomstart = 0,
    start_allocation = start,
    seed = 1
  ))
})

test_that("ADPROCLUS illegal inputs", {
  x <- stackloss

  # no random starts and no start_allocation
  expect_error(adproclus(x,
    nclusters = 2,
    nrandomstart = 0, nsemirandomstart = 0,
    seed = 1
  ))
  # A0 more rows than data
  start <- get_rational(x, x[1:4, ])$A
  expect_error(adproclus(x[seq_len(nrow(x) - 1), ],
    nclusters = 4,
    start_allocation = start,
    seed = 1
  ))

  # ncol(start_allocation) unequal nclusters
  start <- get_rational(x, x[1:4, ])$A
  expect_error(adproclus(x, nclusters = 2, start_allocation = start, seed = 1))
})

test_that("adproclus_low_dim base case normal input", {
  x <- stackloss
  expect_no_error(adproclus_low_dim(x, nclusters = 2, ncomponents = 1, seed = 1))
  expect_no_error(adproclus_low_dim(x,
    nclusters = 3, ncomponents = 2,
    nrandomstart = 1, nsemirandomstart = 1,
    seed = 1
  ))
  expect_no_error(adproclus_low_dim(x,
    nclusters = 3, ncomponents = 2,
    nrandomstart = 2, nsemirandomstart = 2, save_all_starts = TRUE,
    seed = 1
  ))
  expect_no_error(adproclus_low_dim(x, nclusters = 1, ncomponents = 1, seed = 1))
})

test_that("adproclus_low_dim with start_allocation ", {
  x <- stackloss
  start <- get_rational(x, x[1:3, ])$A
  expect_no_error(adproclus_low_dim(x, nclusters = 3, ncomponents = 1,
                                    start_allocation = start,
                                    seed = 1))
  expect_no_error(adproclus_low_dim(x,
    nclusters = 3, ncomponents = 1,
    nrandomstart = 0, nsemirandomstart = 0,
    start_allocation = start,
    seed = 1
  ))
})

test_that("adproclus_low_dim illegal inputs", {
  x <- stackloss

  # no random starts and no start_allocation
  expect_error(adproclus_low_dim(x,
    nclusters = 2, ncomponents = 1,
    nrandomstart = 0, nsemirandomstart = 0
  ))
  # A0 more rows than data
  start <- get_rational(x, x[1:4, ])$A
  expect_error(adproclus_low_dim(x[seq_len(nrow(x) - 1), ], nclusters = 4, ncomponents = 1, start_allocation = start))

  # ncol(start_allocation) unequal nclusters
  start <- get_rational(x, x[1:4, ])$A
  expect_error(adproclus_low_dim(x, nclusters = 2, ncomponents = 1, start_allocation = start))
})

test_that("reproducibility both functions", {
  x <- stackloss
  start <- get_rational(x, x[1:4, ])$A
  expect_equal(
    adproclus(x,
      nclusters = 4, nrandomstart = 1, nsemirandomstart = 1,
      start_allocation = start, save_all_starts = TRUE, seed = 10
    )$model,
    adproclus(x,
      nclusters = 4, nrandomstart = 1, nsemirandomstart = 1,
      start_allocation = start, save_all_starts = TRUE, seed = 10
    )$model
  )
  expect_equal(
    adproclus_low_dim(x,
      nclusters = 4, ncomponents = 1, nrandomstart = 1, nsemirandomstart = 1,
      start_allocation = start, save_all_starts = TRUE, seed = 10
    )$model,
    adproclus_low_dim(x,
      nclusters = 4, ncomponents = 1, nrandomstart = 1, nsemirandomstart = 1,
      start_allocation = start, save_all_starts = TRUE, seed = 10
    )$model
  )
})


test_that("order of clusters", {
  x <- stackloss
  model <- adproclus(x,
    nclusters = 4, nrandomstart = 1, nsemirandomstart = 1,
    save_all_starts = TRUE, seed = 10
  )
  expect_equal(
    unname(rank((-1) * colSums(model$A), ties.method = "first")),
    seq_len(ncol(model$A))
  )
  model2 <- adproclus(x,
    nclusters = 4, nrandomstart = 1, nsemirandomstart = 1,
    save_all_starts = TRUE, seed = 10, algorithm = "ALS2"
  )
  expect_equal(
    unname(rank((-1) * colSums(model2$A), ties.method = "first")),
    seq_len(ncol(model2$A))
  )
  modelLD <- adproclus_low_dim(x,
    nclusters = 4, ncomponents = 2, nrandomstart = 1, nsemirandomstart = 1,
    save_all_starts = TRUE, seed = 10
  )
  expect_equal(
    unname(rank((-1) * colSums(modelLD$A), ties.method = "first")),
    seq_len(ncol(modelLD$A))
  )
})
