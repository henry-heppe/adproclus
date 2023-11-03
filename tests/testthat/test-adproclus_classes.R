test_that("adpc constructor test", {
  model_complete <- adproclus(stackloss, 2)
  model_completeLD <- adproclus_low_dim(stackloss, 3, 2)
  expect_no_condition(adpc(model_complete$A, model_complete$P))
  expect_no_condition(adpc(
    A = model_completeLD$A, P = model_completeLD$P,
    C = model_completeLD$C, B = model_completeLD$B
  ))
  expect_no_condition(summary(adpc(model_complete$A, model_complete$P)))
  expect_no_condition(summary(adpc(
    A = model_completeLD$A, P = model_completeLD$P,
    C = model_completeLD$C, B = model_completeLD$B
  )))
})

test_that("print test", {
  x <- stackloss
  model <- adproclus(x, nclusters = 2)
  modelLD <- adproclus_low_dim(x, nclusters = 3, ncomponents = 2)
  expect_no_error(print.adpc(model))
  expect_no_error(print.adpc(modelLD))
  expect_equal(
    print(model, digits = 2, matrix_rows = 1, matrix_cols = 4),
    print.adpc(model, digits = 2, matrix_rows = 1, matrix_cols = 4)
  )
  expect_equal(
    print(modelLD, digits = 2, matrix_rows = 1, matrix_cols = 4),
    print.adpc(modelLD, digits = 2, matrix_rows = 1, matrix_cols = 4)
  )
})

test_that("summary.adpc test", {
  x <- stackloss
  model <- adproclus(x, nclusters = 2)
  modelLD <- adproclus_low_dim(x, nclusters = 3, ncomponents = 2)

  expect_no_error(summary.adpc(model))
  expect_no_error(summary.adpc(modelLD))
  expect_equal(summary(model), summary.adpc(model))
  expect_equal(summary(modelLD), summary.adpc(modelLD))
})

test_that("print.summary.adpc test", {
  x <- stackloss
  model <- adproclus(x, nclusters = 2)
  modelLD <- adproclus_low_dim(x, nclusters = 3, ncomponents = 2)
  sum_res <- summary(model, digits = 2, matrix_rows = 1, matrix_cols = 4)
  sum_resLD <- summary(modelLD, digits = 2, matrix_rows = 1, matrix_cols = 4)
  expect_error(print.summary.adpc(model))
  expect_no_error(print.summary.adpc(sum_res))
  expect_no_error(print.summary.adpc(sum_resLD))
  expect_equal(print(sum_res), print.summary.adpc(sum_res))
  expect_equal(print(sum_resLD), print.summary.adpc(sum_resLD))
})

test_that("plot test", {
  x <- stackloss
  model <- adproclus(x, nclusters = 2)
  modelLD <- adproclus_low_dim(x, nclusters = 3, ncomponents = 2)

  expect_no_condition(plot(model))
  expect_no_condition(plot(modelLD))
  expect_no_condition(plot(model, type = "Network"))
  expect_no_condition(plot(modelLD, type = "Network"))
  expect_no_condition(plot(model, type = "Profiles"))
  expect_no_condition(plot(modelLD, type = "Profiles"))
  expect_error(plot(model, type = "vars_by_comp"))
  expect_no_condition(plot(modelLD, type = "vars_by_comp"))
  expect_equal(plot(model), plot.adpc(model))
  expect_equal(plot(modelLD), plot.adpc(modelLD))
})
