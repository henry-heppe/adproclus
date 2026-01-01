test_that("network plot, basic functionality", {
        model <- adproclus(stackloss,
                           nclusters = 4,
                           nrandomstart = 1, nsemirandomstart = 1,
                           seed = 1
        )


        modelLD <- adproclus_low_dim(stackloss, 3, 2, seed = 1)
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
                           nrandomstart = 1, nsemirandomstart = 1,
                           seed = 1
        )
        modelLD <- adproclus_low_dim(stackloss, 3, 2, seed = 1)
        expect_no_condition(plot_profiles(model))
        expect_no_condition(plot_profiles(modelLD))
})

test_that("VarsByComp plot, basic functionality", {
        modelLD <- adproclus_low_dim(stackloss, 3, 2, seed = 1)
        expect_no_condition(plot_vars_by_comp(modelLD))
})

test_that("Scree plot full dimensional", {
        model_selection <- mselect_adproclus(stackloss, min_nclusters = 1, max_nclusters = 4, seed = 1)
        expect_no_condition(plot_scree_adpc(model_selection))
})

test_that("Scree plots low dimensional", {

        model_selection <- mselect_adproclus_low_dim(stackloss,
                                                     min_nclusters = 1, max_nclusters = 4,
                                                     min_ncomponents = 1, max_ncomponents = 4,
                                                     seed = 1)
        expect_no_condition(plot_scree_adpc(model_selection))
        expect_no_error(plot_scree_adpc(model_selection, grid = TRUE))
})
