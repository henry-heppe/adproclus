test_that("ADPROCLUS base case normal input", {
        x <- ADPROCLUS::CGdata[1:100,]
        expect_no_error(adproclus(data = x, nclusters= 2))
        expect_no_error(adproclus(data = x, nclusters = 3,
                                  nrandomstart = 1, nsemirandomstart = 1, algorithm = "ALS2"))
        expect_no_error(adproclus(data = x, nclusters = 3,
                                  nrandomstart = 1, nsemirandomstart = 1, algorithm = "ALS1"))
        expect_no_error(adproclus(data = x, nclusters= 3,
                                  nrandomstart = 2, nsemirandomstart = 2, saveAllStarts = TRUE))
})

test_that("ADPROCLUS with start_allocation ", {
        x <- ADPROCLUS::CGdata[1:100,]
        start <- getRational(x,x[1:4,])$A
        expect_no_error(adproclus(data = x, nclusters= 4, start_allocation = start))
        expect_no_error(adproclus(data = x, nclusters= 4,
                                  nrandomstart = 0, nsemirandomstart = 0,
                                  start_allocation = start))
})

test_that("ADPROCLUS illegal inputs", {
        x <- ADPROCLUS::CGdata[1:100,]

        #no random starts and no start_allocation
        expect_error(adproclus(data = x, nclusters= 2,
                                  nrandomstart = 0, nsemirandomstart = 0))
        #A0 more rows than data
        start <- getRational(x,x[1:4,])$A
        expect_error(adproclus(data = x[1:nrow(data)-1,], nclusters= 4, start_allocation = start))

        #ncol(start_allocation) unequal nclusters
        start <- getRational(x,x[1:4,])$A
        expect_error(adproclus(data = x, nclusters= 2, start_allocation = start))
})

test_that("adproclusLD base case normal input", {
        x <- ADPROCLUS::CGdata[1:100,]
        expect_no_error(adproclusLD(data = x, nclusters = 2, ncomponents = 1))
        expect_no_error(adproclusLD(data = x, nclusters = 3, ncomponents = 2,
                                  nrandomstart = 1, nsemirandomstart = 1))
        expect_no_error(adproclusLD(data = x, nclusters= 3, ncomponents = 2,
                                  nrandomstart = 2, nsemirandomstart = 2, saveAllStarts = TRUE))
})

test_that("adproclusLD with start_allocation ", {
        x <- ADPROCLUS::CGdata[1:100,]
        start <- getRational(x,x[1:3,])$A
        expect_no_error(adproclusLD(data = x, nclusters= 3, ncomponents = 1, start_allocation = start))
        expect_no_error(adproclusLD(data = x, nclusters= 3, ncomponents = 1,
                                  nrandomstart = 0, nsemirandomstart = 0,
                                  start_allocation = start))
})

test_that("adproclusLD illegal inputs", {
        x <- ADPROCLUS::CGdata[1:100,]

        #no random starts and no start_allocation
        expect_error(adproclusLD(data = x, nclusters= 2, ncomponents = 1,
                               nrandomstart = 0, nsemirandomstart = 0))
        #A0 more rows than data
        start <- getRational(x,x[1:4,])$A
        expect_error(adproclusLD(data = x[1:nrow(data)-1,], nclusters= 4, ncomponents = 1, start_allocation = start))

        #ncol(start_allocation) unequal nclusters
        start <- getRational(x,x[1:4,])$A
        expect_error(adproclusLD(data = x, nclusters= 2, ncomponents = 1, start_allocation = start))
})

test_that("reproducibility both functions", {
        x <- ADPROCLUS::CGdata[1:100,]
        start <- getRational(x,x[1:4,])$A
        expect_equal(adproclus(data = x, nclusters= 4, nrandomstart = 1, nsemirandomstart = 1,
                               start_allocation = start, saveAllStarts = TRUE, seed = 10)$model,
                     adproclus(data = x, nclusters= 4, nrandomstart = 1, nsemirandomstart = 1,
                               start_allocation = start, saveAllStarts = TRUE, seed = 10)$model)
        expect_equal(adproclusLD(data = x, nclusters= 4, ncomponents = 1, nrandomstart = 1, nsemirandomstart = 1,
                               start_allocation = start, saveAllStarts = TRUE, seed = 10)$model,
                     adproclusLD(data = x, nclusters= 4, ncomponents = 1, nrandomstart = 1, nsemirandomstart = 1,
                               start_allocation = start, saveAllStarts = TRUE, seed = 10)$model)
})

#issue: add test for (re-) naming of columns/rows

#other test ideas
test_that("ADPROCLUS S3 functionality", {

})

test_that("ADPROCLUS edge cases", {

})


