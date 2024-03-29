test_that("test random start", {
  x <- adproclus::CGdata
  expect_no_error(get_random(x, 2))
})

test_that("test semi-random start", {
  x <- adproclus::CGdata
  expect_no_error(get_semirandom(x, 2))
})

test_that("test rational start", {
  x <- adproclus::CGdata
  expect_no_error(get_rational(x, x[1:2, ]))
})

test_that("reproducibility check for all three", {
  x <- adproclus::CGdata
  expect_equal(get_random(x, 2, seed = 1), get_random(x, 2, seed = 1))
  expect_equal(get_semirandom(x, 2, seed = 1), get_semirandom(x, 2, seed = 1))
  expect_equal(get_rational(x, x[1:2, ]), get_rational(x, x[1:2, ]))
})
