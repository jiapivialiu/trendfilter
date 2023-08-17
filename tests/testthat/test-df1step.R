test_that("df1steps coincide", {
  k <- 2
  lambda <- 2
  y <- rnorm(1)
  A <- rbind(-rev(dspline::d_mat(k, seq(k + 1)))[-1], diag(1, k, k)[-k, ])
  Z <- matrix(c(1, rep(0, k - 1)), nrow = 1)
  H <- matrix(1)
  R <- matrix(0, k, 1)
  R[1] <- 1
  Q <- matrix(1 / lambda, 1, 1)
  RQR <- R %*% Q %*% t(R)
  P1 <- matrix(rnorm(k^2), k, k)
  P1 <- crossprod(P1)
  P1inf <- diag(1, k, k)
  a1 <- as.matrix(double(k))
  rankp <- k

  Rcpp <- df1step_test(y, Z, H, A, RQR, a1, P1, P1inf, rankp)
  R <- df1step(y, Z, H, A, RQR, a1, P1, P1inf, rankp)
  nms <- names(Rcpp)
  for (i in seq_along(Rcpp)) {
    expect_equal(drop(Rcpp[[nms[i]]]), drop(R[[nms[i]]]))
  }
})
