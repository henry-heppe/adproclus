test_that("model selection, full dimensional", {
  expect_no_condition(mselect_adproclus(stackloss,
                                        min_nclusters = 1, max_nclusters = 4,
                                        return_models = TRUE,
                                        seed = 1))
        expect_no_condition(mselect_adproclus(stackloss,
                                              min_nclusters = 1, max_nclusters = 4,
                                              return_models = FALSE,
                                              seed = 1))
})

test_that("model selection, low dimensional", {
        expect_no_condition(mselect_adproclus_low_dim(stackloss,
                                                      min_nclusters = 1, max_nclusters = 4,
                                                      min_ncomponents = 1, max_ncomponents = 4,
                                                      return_models = TRUE,
                                                      seed = 1))

        expect_no_condition(mselect_adproclus_low_dim(stackloss,
                                                      min_nclusters = 1, max_nclusters = 4,
                                                      min_ncomponents = 1, max_ncomponents = 4,
                                                      return_models = FALSE,
                                                      seed = 1))
})
