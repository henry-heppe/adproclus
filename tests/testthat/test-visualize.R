test_that("network plot, basic functionality", {
        model <- adproclus(adproclus::CGdata[1:100,], nclusters = 4, nrandomstart = 1, nsemirandomstart = 1)
        modelLD <- adproclusLD(adproclus::CGdata[1:100,], 3, 2)
        expect_no_condition(plot_cluster_network(model))
        expect_no_error(plot_cluster_network(model, title = "Test network", relative_overlap = FALSE, filetype = "pdf"))
        expect_no_condition(plot_cluster_network(modelLD))
        expect_no_error(plot_cluster_network(modelLD, title = "Test network", relative_overlap = FALSE, filetype = "pdf"))
        expect_no_condition(plot(adproclusLD(stackloss, 3, 1, seed = 1), type = "Network"))
})

test_that("profile plot, basic functionality", {
        model <- adproclus(adproclus::CGdata[1:100,], nclusters = 4, nrandomstart = 1, nsemirandomstart = 1)
        modelLD <- adproclusLD(adproclus::CGdata[1:100,], 3, 2)
        expect_no_condition(plot_profiles(model))
        expect_no_condition(plot_profiles(modelLD))
        })

test_that("VarsByComp plot, basic functionality", {
        modelLD <- adproclusLD(adproclus::CGdata[1:100,], 3, 2)
        expect_no_condition(plot_vars_by_comp(modelLD))
})
