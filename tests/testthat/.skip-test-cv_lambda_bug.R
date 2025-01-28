test_that("cv_trendfilter accepts `lambda`", {
  x <- 1:100
  y <- sin(x / 101 * 2 * pi) + rnorm(100, 0, .2)
  tf <- trendfilter(y, x)
  expect_s3_class(cv_trendfilter(y, x, nfolds = 3), "cv_trendfilter")
  expect_s3_class(cv_trendfilter(y, x, nfolds = 3, lambda = tf$lambda), "cv_trendfilter")
  expect_equal(cv_trendfilter(2 * y, x, nfolds = 3, lambda = tf$lambda)$lambda, tf$lambda)
})
