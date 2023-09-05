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
