gen_Deltahat_fun_g <- function(m) {
  Deltahat_fun_g <- function(W, theta) {

    mWtheta <- purrr::modify(W, m, theta=theta)
    mbar <- colMeans(mWtheta)
    Sigmahat <- tcrossprod(mWtheta - mbar) / n
    sigmahat <- sqrt(diag(Sigmahat))
    return(max(-mbar/sigmahat))
  }
  return(Deltahat_fun_g)
}



gen_Astarinf_obj_g <- function(W, m, n, kappa, rhatinf, iota_sd) {
  Astarinf_obj <- function(Wstar, theta) {

    #original data summary stats
    mWtheta <- purrr::modify(W, m, theta=theta)
    mbar <- colMeans(mWtheta)
    Sigmahat <- tcrossprod(mWtheta - mbar) / n
    sigmahat <- sqrt(diag(Sigmahat))
    mhat <- mbar/sigmahat

    #bootstrap sample:
    mWstar <- purrr::modify(Wstar, m, theta=theta)
    mbar_star <- colMeans(mWstar)
    Sigmastar <- tcrossprod(mWstar - mbar_star) / n
    sigmastar <- sqrt(diag(Sigma_star))
    mhat_star <- mbar_star/sigma_star



    nuhat_star <- sqrt(n) *
      ((mbar_star / sigmastar) - mhat)

    chistartheta <- ifelse(nuhat_star >= 0,
                           chi(nuhat_star, sqrt(n) * mhat - (kappa * sd2)),
                           chi(nuhat_star, sqrt(n) * mhat + (kappa * sd2)))


    # calculate sd2, sd3
    Vhatstar2 <- n * matrixStats::rowVars(mhatstarstars)
    sd2 <- sqrt(Vhatstar2)
    sd2[sd2 < iota_sd] <- iota_sd

    negmhatstars <- neg(mhatstarstars)
    maxmhatss <- matrixStats::colMaxs(negmhatstars)
    Vhatstar3 <- n * matrixStats::colVars(t(negmhatstars) - maxmhatss)
    sd3 <- sqrt(Vhatstar3)
    sd3[sd3 < iota_sd] <- iota_sd

    chistartheta <- ifelse(nuhat_star >= 0,
                           chi(nuhat_star, sqrt(n) * mhattheta - (kappa * sd2)),
                           chi(nuhat_star, sqrt(n) * mhattheta + (kappa * sd2)))

    bhattheta <- sqrt(n) * (neg(mhattheta) - rhatinf) - kappa * sd3
    xibtheta <- sqrt(n) * (neg(mhattheta) - rhatinf) / (kappa * sd3)
    phitheta <- phi(xibtheta)

    k <- length(mhattheta)
    temp <- rep(NA, k)
    for (j1 in 1:k) {
      temp[j1] <- max(chistartheta +
                        ifelse((1:k) == j1, phitheta, bhattheta))
    }

    jhat <- which(neg(mhattheta) >= max(neg(mhattheta)) - sd3 * kappa / sqrt(n))
    ret <- min(temp[jhat])

    if (ret == Inf) {ret <- 10^10}
    return(ret)
  }

  return(Astarinf_obj)
}
