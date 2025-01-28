test_that("admm control lists check args", {
  expect_snapshot(error = TRUE, admm_control_list(max_iter = 0L))
  expect_snapshot(error = TRUE, admm_control_list(max_iter = 1.1))
  expect_snapshot(error = TRUE, admm_control_list(max_iter = 100:200))

  expect_snapshot(error = TRUE, admm_control_list(rho_scale = 0))
  expect_snapshot(error = TRUE, admm_control_list(rho_scale = Inf))
  expect_snapshot(error = TRUE, admm_control_list(rho_scale = c(1e-6, 1e-5)))

  expect_snapshot(error = TRUE, admm_control_list(tolerance = 0))
  expect_snapshot(error = TRUE, admm_control_list(tolerance = Inf))
  expect_snapshot(error = TRUE, admm_control_list(tolerance = c(1e-6, 1e-5)))

  expect_snapshot(error = TRUE, admm_control_list(unknown = 1))
  expect_s3_class(admm_control_list(), "admm_control")
})

test_that("trendfilter control lists check args", {
  expect_snapshot(error = TRUE, trendfilter_control_list(obj_tol = 0))
  expect_snapshot(error = TRUE, trendfilter_control_list(obj_tol = Inf))
  expect_snapshot(error = TRUE, trendfilter_control_list(obj_tol = c(1e-6, 1e-5)))

  expect_snapshot(error = TRUE, trendfilter_control_list(x_cond = 0))
  expect_silent(trendfilter_control_list(x_cond = 1))
  expect_snapshot(error = TRUE, trendfilter_control_list(x_cond = Inf))
  expect_snapshot(error = TRUE, trendfilter_control_list(x_cond = c(5, 10)))

  expect_snapshot(error = TRUE, trendfilter_control_list(admm_control = list()))
  expect_snapshot(error = TRUE, trendfilter_control_list(unknown = 1))
  expect_s3_class(trendfilter_control_list(), "trendfilter_control")
})
