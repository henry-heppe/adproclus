test_that("test random start", {
        x <- ADPROCLUS::CGdata
        expect_no_error(getRandom(x, 6))
})

test_that("test semi-random start", {
        x <- ADPROCLUS::CGdata
        expect_no_error(getSemiRandom(x, 6))
        })

test_that("test rational start", {
        x <- ADPROCLUS::CGdata
        expect_no_error(getRational(x, x[1:6,]))
        })

test_that("reproducibility check for all three", {
        withr::local_seed(1)
        #set.seed(1)
        x <- ADPROCLUS::CGdata
        expect_equal(getRandom(x, 6),getRandom(x, 6))
        expect_equal(getSemiRandom(x, 6),getSemiRandom(x, 6))
        expect_equal(getRational(x, x[1:6,]),getRational(x, x[1:6,]))

})
