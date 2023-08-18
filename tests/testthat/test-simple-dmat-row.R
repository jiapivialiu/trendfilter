test_that("simplistic dmat construction works", {

  f <- function(k) {
    # mimics C++ loop
    x = double(k)
    s = -1
    for (i in 1:k) {
      s = s * -1
      x[i] = choose(k, i) * s
    }
    return(x)
  }

  g <- function(k) -rev(dspline::d_mat(k, seq(k + 1)))[-1]

  for (j in 1:10) expect_equal(f(!!j), g(!!j))

})
