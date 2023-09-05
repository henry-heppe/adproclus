test_that("network plot, basic functionality", {
        model <- adproclus(ADPROCLUS::CGdata[1:100,], nclusters = 4, nrandomstart = 1, nsemirandomstart = 1)
        #expect_no_error(plot_clusters_as_network(model))
})
