test_that("kftf cpp implementation matches KFAS", {
  y <- rnorm(100)
  a1 <- kftfr(y, 3, 10)
  a2 <- kftf_cpp_test(y, 3, 10)
  expect_equal(a1, a2)
  expect_equal(kftfr(y, 5, 100), kftf_cpp_test(y, 5, 100))
})
