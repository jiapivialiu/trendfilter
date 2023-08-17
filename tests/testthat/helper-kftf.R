kftf <- function(y, k, lambda) {
  # correct diffuse initialization using {KFAS}
  n <- length(y)
  A <- rbind(-rev(dspline::d_mat(k, seq(k + 1)))[-1], diag(1, k, k)[-k, ])
  Z <- matrix(c(1, rep(0, k - 1)), nrow = 1)
  H <- matrix(1) # if weighted, then this is h_i = 1 / w_i
  R <- matrix(0, k, 1)
  R[1] <- 1
  Q <- matrix(1/lambda, 1, 1)
  a1 <- rep(0, k)
  P1 <- matrix(0, k, k)
  P1inf <- diag(1, k, k)

  o <- KFAS::SSModel(
    y ~ -1 + 
      KFAS::SSMcustom(Z, A, R, Q, a1, P1, P1inf, n = n), H = H
    )
  oo <- KFAS::KFS(o)
  as.numeric(drop(oo$alphahat[,1]))
}


kftfr <- function(y, k, lambda) {
  # correct diffuse initialization by hand
  n <- length(y)
  A <- rbind(-rev(dspline::d_mat(k, seq(k + 1)))[-1], diag(1, k, k)[-k, ])
  Z <- matrix(c(1, rep(0, k - 1)), nrow = 1)
  H <- matrix(1)
  R <- matrix(0, k, 1)
  R[1] <- 1
  Q <- matrix(1/lambda, 1, 1)
  at <- matrix(0, k, n + 1)
  Pt <- array(0, c(k, k, n + 1))
  vt <- double(n)
  Ft <- double(n)
  Kt <- matrix(0, k, n)
  Kinf <- matrix(0, k, n)
  RQR <- R %*% Q %*% t(R)
  P1 <- matrix(0, k, k)
  P1inf <- diag(1, k, k)
  a1 <- double(k)

  d <- 0
  Pinf <- array(0, c(k, k, n + 1))
  Finf <- double(n)

  rankp <- k
  Pt[,,1] <- P1
  Pinf[,,1] <- P1inf
  at[,1] <- a1

  while (rankp > 0 && d < n) {
    d <- d + 1
    filt <- df1step(y[d], Z, H, A, RQR, at[,d], Pt[,,d], Pinf[,,d], rankp)
    vt[d] <- filt$vt
    Ft[d] <- filt$Ft
    Finf[d] <- filt$Finf
    Kt[,d] <- A %*% filt$Mt
    Kinf[,d] <- A %*% filt$Minf
    at[,d + 1] <- filt$a
    Pt[,,d + 1] <- filt$P
    Pinf[,,d + 1] <- filt$Pinf
    rankp <- filt$rankp
  }

  for (i in (d + 1):n) {
    filt <- f1step(y[i], Z, H, A, RQR, at[,i], Pt[,,i])
    vt[i] <- filt$vt
    Ft[i] <- filt$Ft
    Kt[,i] <- filt$Kt
    at[,i + 1] <- filt$a
    Pt[,,i + 1] <- filt$P
  }

  r <- double(k)
  ahat <- matrix(0, k, n)
  for (i in n:(d + 1)) {
    L <- A - Kt[,i] %*% Z
    r <- t(Z) * vt[i] / Ft[i] + t(L) %*% r
    ahat[,i] <- at[,i] + Pt[,,i] %*% r
  }

  r1 <- double(k)
  ZZ <- crossprod(Z)
  for (i in d:1) {
    L0 <- A - A %*% Pinf[,,i] %*% ZZ / Finf[i]
    L1 <- -A %*% Pt[,,i] %*% ZZ / Finf[i] +
      A %*% Pinf[,,i] %*% ZZ * Ft[i] / Finf[i]^2
    r1 <- t(Z) %*% vt[i] / Finf[i] + t(L0) %*% r1 + t(L1) %*% r
    r <- t(L0) %*% r
    ahat[,i] <- at[,i] + Pt[,,i] %*% r + Pinf[,,i] %*% r1
  }
  ahat[1,]
}


f1step <- function(y, Z, H, A, RQR, a, P) {
  ZZ <- crossprod(Z)
  vt <- drop(y - Z %*% a)
  Ft <- drop(Z %*% P %*% t(Z) + H)
  Kt <- A %*% P %*% t(Z) / Ft
  a <- A %*% a + Kt * vt
  Ptt <- P - P %*% ZZ %*% P / Ft
  P <- A %*% Ptt %*% t(A) + RQR
  return(list(a = a, P = P, vt = vt, Ft = Ft, Kt = Kt))
}


df1step <- function(y, Z, H, A, RQR, a, P, Pinf, rankp) {
  tol <- sqrt(.Machine$double.eps)
  Mt <- P %*% t(Z)
  Ft <- drop(Z %*% P %*% t(Z) + H)
  Minf <- Pinf %*% t(Z)
  Finf <- drop(Z %*% Pinf %*% t(Z))
  vt <- drop(y - Z %*% a)

  if (Finf > tol) {
    # should always happen
    finv <- 1 / Finf
    a <- a + vt * finv * Minf
    P <- P + Ft * finv^2 * tcrossprod(Minf)
    P <- P - finv * tcrossprod(Mt, Minf) - finv * tcrossprod(Minf, Mt)
    Pinf <- Pinf - finv * tcrossprod(Minf)
    rankp <- rankp - 1
  } else {
    # should never happen
    Finf <- 0
    if (Ft > tol) {
      finv <- 1 / Ft
      a <- a + vt * finv * Mt
      P <- P - finv * tcrossprod(Mt)
    }
  }
  if (Ft < tol) Ft <- 0

  a <- A %*% a
  P <- A %*% P %*% t(A) + RQR
  Pinf <- A %*% Pinf %*% t(A)

  # Fix possible negative definiteness, should never happen
  pneg <- diag(P) < 0
  pinfneg <- diag(Pinf) < 0
  if (any(pneg)) {
    P[pneg,] <- 0
    P[,pneg] <- 0
  }
  if (any(pinfneg)) {
    Pinf[pinfneg,] <- 0
    Pinf[,pinfneg] <- 0
  }
  return(list(
    vt = vt, Ft = Ft, Finf = Finf, Mt = Mt, Minf = Minf,
    a = a, P = P, Pinf = Pinf, rankp = rankp
  ))
}
