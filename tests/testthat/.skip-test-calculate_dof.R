test_that("knot calculation works", {
  k <- 2
  n <- 20
  x <- 1:n
  alpha <- 3*(x < 6) - 1*(x > 5) + 6*(x > 15)
  expect_equal(calc_degrees_of_freedom(alpha, k, 1e-8), 2 + k + 1)
})
