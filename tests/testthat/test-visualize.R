test_that("network plot, basic functionality", {
        model <- adproclus(ADPROCLUS::CGdata[1:100,], nclusters = 4, nrandomstart = 1, nsemirandomstart = 1)
        modelLD <- adproclusLD(ADPROCLUS::CGdata[1:100,], 3, 2)
        expect_no_condition(plotClusterNetwork(model))
        expect_no_error(plotClusterNetwork(model, title = "Test network", relative_overlap = FALSE, filetype = "pdf"))
        expect_no_condition(plotClusterNetwork(modelLD))
        expect_no_error(plotClusterNetwork(modelLD, title = "Test network", relative_overlap = FALSE, filetype = "pdf"))
        expect_no_condition(plot(adproclusLD(stackloss, 3, 1, seed = 1), type = "Network"))
})

test_that("profile plot, basic functionality", {
        model <- adproclus(ADPROCLUS::CGdata[1:100,], nclusters = 4, nrandomstart = 1, nsemirandomstart = 1)
        modelLD <- adproclusLD(ADPROCLUS::CGdata[1:100,], 3, 2)
        expect_no_condition(plotProfiles(model))
        expect_no_condition(plotProfiles(modelLD))
        })

test_that("VarsByComp plot, basic functionality", {
        modelLD <- adproclusLD(ADPROCLUS::CGdata[1:100,], 3, 2)
        expect_no_condition(plotVarsByComp(modelLD))
})
