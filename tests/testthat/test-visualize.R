test_that("network plot, basic functionality", {
        model <- adproclus(stackloss,
                           nclusters = 4,
                           nrandomstart = 1, nsemirandomstart = 1
        )
        modelLD <- adproclus_low_dim(stackloss, 3, 2)
        expect_no_condition(plot_cluster_network(model))
        expect_no_error(plot_cluster_network(model,
                                             title = "Test network",
                                             relative_overlap = FALSE
        ))
        expect_no_condition(plot_cluster_network(modelLD))
        expect_no_error(plot_cluster_network(modelLD,
                                             title = "Test network",
                                             relative_overlap = FALSE
        ))
        model_stackloss <- adproclus_low_dim(stackloss, 3, 1, seed = 1)
        expect_no_condition(plot(model_stackloss, type = "Network"))
})

test_that("profile plot, basic functionality", {
        model <- adproclus(stackloss,
                           nclusters = 4,
                           nrandomstart = 1, nsemirandomstart = 1
        )
        modelLD <- adproclus_low_dim(stackloss, 3, 2)
        expect_no_condition(plot_profiles(model))
        expect_no_condition(plot_profiles(modelLD))
})

test_that("VarsByComp plot, basic functionality", {
        modelLD <- adproclus_low_dim(stackloss, 3, 2)
        expect_no_condition(plot_vars_by_comp(modelLD))
})
