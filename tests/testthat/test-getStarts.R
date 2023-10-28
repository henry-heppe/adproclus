test_that("test random start", {
        x <- adproclus::CGdata
        expect_no_error(getRandom(x, 6))
})

test_that("test semi-random start", {
        x <- adproclus::CGdata
        expect_no_error(getSemiRandom(x, 6))
        })

test_that("test rational start", {
        x <- adproclus::CGdata
        expect_no_error(getRational(x, x[1:6,]))
        })

test_that("reproducibility check for all three", {
        x <- adproclus::CGdata
        expect_equal(getRandom(x, 6, seed = 1),getRandom(x, 6, seed = 1))
        expect_equal(getSemiRandom(x, 6, seed = 1),getSemiRandom(x, 6, seed = 1))
        expect_equal(getRational(x, x[1:6,]),getRational(x, x[1:6,]))

})
