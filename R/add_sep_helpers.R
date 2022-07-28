# Functions used in the case where m is additively separable


#' mhat calculation
#'
#' @param theta
#' @param gWbar Column Means of the gW matrix or the mean of g(W) over n
#' @param sigmahat Standard Deviation
#' @param f The function f(theta) that is a component of m(W,theta)
#' @param rho
#'
#' @return
#' @export
#'
#' @examples
mhat_fun <- function(theta, gWbar, sigmahat, f, rho = NULL) {
  ftheta <- f(theta)
  (gWbar + ftheta) / sigmahat
}

Deltahat_fun <- function(theta, gWbar, sigmahat, rho = NULL) {
  if (is.null(rho)) rho <- theta[9] #!
  max(-mhat_fun(theta, gWbar, sigmahat, rho = rho))
}

m_minus_mhat_fun <- function(gW, gWbar, sigmahat) {
  (gW - gWbar) / sigmahat
}

gen_Astarinf_obj <- function(gWbar, sigmahat, mhat_fun, Deltahatinf,
                             sqrt_ohplus_simZ, n, kappa) {
  Astarinf_obj <- function(theta, Wstarbar, sigmastar) {

    # Sample moments
    mhattheta <- mhat_fun(theta, Wbar, sigmahat)
    mhatstarstars <- mhat_fun(theta, sd_Wstarbars, sd_sigmastars)
    nuhat_star <- sqrt(n) * (mhat_fun(theta, Wstarbar, sigmastar) - mhattheta)

    # calculate sd2, sd3
    Vhatstar2 <- n * matrixStats::rowVars(mhatstarstars)
    sd2 <- sqrt(pmax1(Vhatstar2))

    negmhatstars <- neg(mhatstarstars)
    maxmhatss <- matrixStats::colMaxs(negmhatstars)
    Vhatstar3 <- n * matrixStats::colVars(t(negmhatstars) - maxmhatss)
    sd3 <- sqrt(pmax1(Vhatstar3))

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


gen_Astarinf_Deltas_obj <- function(gWbar, sigmahat, mhat_fun, Deltahatinf,
                                    sqrt_ohplus_simZ, n, kappa) {

  Astarinf_Deltas_obj <- function(ftheta, gWstarbar, sigmastar) {

    # Sample moments
    mhattheta <- mhat_fun(ftheta, gWbar, sigmahat)
    nuhat_star <- sqrt(n) * (mhat_fun(ftheta, gWstarbar, sigmastar) - mhattheta)

    # calculate sdhat
    Ghats <- sqrt_ohplus_simZ %*%
      rbind(diag(length(mhattheta)), diag(-.5 * mhattheta))
    maxGhats <- matrixStats::rowMaxs(Ghats)
    V_S <- matrixStats::colVars(Ghats - maxGhats)
    sdhat <- sqrt(pmax1(V_S))

    ehattheta <- sqrt(n) * ( - mhattheta - Deltahatinf) - sdhat * kappa
    xie <-  (sqrt(n)/(kappa * sdhat)) * ( - mhattheta - Deltahatinf)
    phitheta <- phi(xie)

    k <- length(mhattheta)
    temp <- rep(NA, k)
    for (j1 in 1:k) {
      temp[j1] <- max( - nuhat_star + ifelse((1:k)==j1, phitheta, ehattheta))
    }

    jhat <- which(-mhattheta >= max(-mhattheta) - sdhat * kappa / sqrt(n))
    ret <- min(temp[jhat])

    if (ret == Inf) {ret <- 10^10}
    return(ret)
  }

  return(Astarinf_Deltas_obj)
}


gen_Astarinf_DR_Deltas_obj <- function(gWbar, sigmahat, mhat_fun, Deltahatinf,
                                       sqrt_ohplus_simZ, n, kappa) {

  Astarinf_DR_Deltas_obj <- function(ftheta,  gWstarbar, sigmastar) {
    # Sample moments
    mhattheta <- mhat_fun(ftheta, gWbar, sigmahat)
    nuhat_star <- sqrt(n) * (mhat_fun(ftheta, gWstarbar, sigmastar) - mhattheta)

    # calculate sdhat
    Ghats <- sqrt_ohplus_simZ %*%
      rbind(diag(length(mhattheta)), diag(-.5 * mhattheta))
    maxGhats <- matrixStats::rowMaxs(Ghats)
    V_S <- matrixStats::colVars(Ghats - maxGhats)
    sdhat <- sqrt(pmax1(V_S))

    xie <-  (sqrt(n)/(kappa * sdhat)) * (- mhattheta - Deltahatinf)

    return(max(-nuhat_star - phi(-xie)))
  }

  return(Astarinf_DR_Deltas_obj)
}



gen_Astarinf_PR_Deltas_obj <- function(gWbar, sigmahat, mhat_fun, Deltahatinf,
                                       sqrt_ohplus_simZ, n, kappa) {

  Astarinf_PR_Deltas_obj <- function(ftheta, gWstarbar, sigmastar) {

    # Sample moments
    mhattheta <- mhat_fun(ftheta, gWbar, sigmahat)
    nuhat_star <- sqrt(n) * (mhat_fun(ftheta, gWstarbar, sigmastar) - mhattheta)

    # calculate sdhat
    Ghats <- sqrt_ohplus_simZ %*%
      rbind(diag(length(mhattheta)), diag(-.5 * mhattheta))
    maxGhats <- matrixStats::rowMaxs(Ghats)
    V_S <- matrixStats::colVars(Ghats - maxGhats)
    sdhat <- sqrt(pmax1(V_S))

    xiee <-  (sqrt(n)/(kappa * max(sdhat))) * (- mhattheta - Deltahatinf)

    return(max(-nuhat_star + xiee))
  }

  return(Astarinf_PR_Deltas_obj)
}
