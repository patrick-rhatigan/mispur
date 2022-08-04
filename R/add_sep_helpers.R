# Functions used in the case where m is additively separable


#' mhat calculation
#'
#' @param theta
#' @param gWbar Column Means of the gW matrix or the mean of g(W) over n
#' @param sigmahat Standard Deviation
#' @param f The function f(theta) that is a component of m(W,theta)
#'
#' @return
#' @export
#'
#' @examples
mhat_fun <- function(theta, gWbar, sigmahat, f) {
  ftheta <- f(theta)
  (gWbar + ftheta) / sigmahat
}

Deltahat_fun <- function(theta, gWbar, sigmahat) {
  max(-mhat_fun(theta, gWbar, sigmahat))
}

m_minus_mhat_fun <- function(gW, gWbar, sigmahat) {
  (gW - gWbar) / sigmahat
}

###
const_alt <- function(thetatil, gWbar, sigmahat, f) {
  theta <- thetatil[-1]
  gam <- thetatil[1]
  ftheta <- f(theta)
  gam + (gWbar + ftheta) / sigmahat
}
###


#' Calculate gWbar and sigmahat for use when bootstrapping
#' (Additively Seperable case)
#'
#' @param W Data Matrix
#' @param g the function g(W)
#'
#' @return
#' @export
#'
#' @examples
calc_summ_stats_as <- function(W, g) {
  n = nrow(W)
  k = ncol(W)
  gW <- matrix(NA,nrow = n, ncol = k)
  for(i in 1:n) {
    gW [i,] <- sapply(W[i,], g)
  }
  gWbar = colMeans(gW)

  Sigmahat <- tcrossprod(gW - gWbar) / n
  sigmahat <- sqrt(diag(Sigmahat))
  return(list(gWbar = gWbar, sigmahat = sigmahat))
}

#!!!

#' Generate Astarinf object.
#'
#' Parameters are identical for this whole class of functions,
#' small changes are made in the calculation for different particular
#' cases, I'm not sure how they differ.
#'
#' @param gWbar colmeans of gW
#' @param sigmahat
#' @param rhatinf
#' @param sd_gWstarbars
#' @param sd_sigmastars
#' @param n
#' @param kappa
#' @param f function f(theta)
#' @param iota_sd buffer term (default is 1e-06)
#'
#' @return
#' @export
#'
#' @examples
gen_Astarinf_obj <- function(gWbar, sigmahat, rhatinf, sd_gWstarbars,
                             sd_sigmastars, n, kappa, f, iota_sd) {
  Astarinf_obj <- function(theta, gWstarbar, sigmastar) {

    # Sample moments
    ftheta <- f(theta)
    mhattheta <- (gWbar + ftheta) / sigmahat
    mhatstarstars <- (sd_gWstarbars + ftheta) / sd_sigmastars
    nuhat_star <- sqrt(n) *
      ((gWstarbar + ftheta) / sigmastar - mhattheta)


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

gen_Astarinf_obj_comp <- function(gWbar, sigmahat, rhatinf, sd_gWstarbars,
                                  sd_sigmastars, n, kappa, f,
                                  sd2, sd3) {
  Astarinf_obj <- function(theta, gWstarbar, sigmastar) {

    # Sample moments
    mhattheta <- mhat_fun(theta, gWbar, sigmahat, f)
    nuhat_star <- sqrt(n) *
      (mhat_fun(theta, gWstarbar, sigmastar, rho, px) - mhattheta)


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


gen_Astarinf_Deltas_obj <- function(gWbar, sigmahat, Deltahatinf,
                                    sqrt_ohplus_simZ, n, kappa, f, iota_sd) {

  Astarinf_Deltas_obj <- function(theta, gWstarbar, sigmastar) {

    # Sample moments
    ftheta <- f(theta)
    mhattheta <- (gWbar + ftheta) / sigmahat
    nuhat_star <- sqrt(n) *
      ((gWstarbar + ftheta) / sigmastar - mhattheta)

    # calculate sdhat
    Ghats <- sqrt_ohplus_simZ %*%
      rbind(diag(length(mhattheta)), diag(-.5 * mhattheta))
    maxGhats <- matrixStats::rowMaxs(Ghats)
    V_S <- matrixStats::colVars(Ghats - maxGhats)
    sdhat <- sqrt(V_S)
    sdhat[sdhat < iota_sd] <- iota_sd

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


gen_Astarinf_DR_Deltas_obj <- function(gWbar, sigmahat, Deltahatinf,
                                       sqrt_ohplus_simZ, n, kappa, f, iota_sd) {

  Astarinf_DR_Deltas_obj <- function(theta,  gWstarbar, sigmastar) {
    # Sample moments
    ftheta <- f(theta)
    mhattheta <- (gWbar + ftheta) / sigmahat
    nuhat_star <- sqrt(n) *
      ((gWstarbar + ftheta) / sigmastar - mhattheta)

    # calculate sdhat
    Ghats <- sqrt_ohplus_simZ %*%
      rbind(diag(length(mhattheta)), diag(-.5 * mhattheta))
    maxGhats <- matrixStats::rowMaxs(Ghats)
    V_S <- matrixStats::colVars(Ghats - maxGhats)
    sdhat <- sqrt(V_S)
    sdhat[sdhat < iota_sd] <- iota_sd

    xie <-  (sqrt(n)/(kappa * sdhat)) * (- mhattheta - Deltahatinf)

    return(max(-nuhat_star - phi(-xie)))
  }

  return(Astarinf_DR_Deltas_obj)
}



## gen_Astarinf_PR_Deltas_obj <- function(gWbar, sigmahat, mhat_fun, Deltahatinf,
##                                        sqrt_ohplus_simZ, n, kappa) {

##   Astarinf_PR_Deltas_obj <- function(ftheta, gWstarbar, sigmastar) {

##     # Sample moments
##     mhattheta <- mhat_fun(ftheta, gWbar, sigmahat)
##     nuhat_star <- sqrt(n) * (mhat_fun(ftheta, gWstarbar, sigmastar) - mhattheta)

##     # calculate sdhat
##     Ghats <- sqrt_ohplus_simZ %*%
##       rbind(diag(length(mhattheta)), diag(-.5 * mhattheta))
##     maxGhats <- matrixStats::rowMaxs(Ghats)
##     V_S <- matrixStats::colVars(Ghats - maxGhats)
##     sdhat <- sqrt(pmax1(V_S))

##     xiee <-  (sqrt(n)/(kappa * max(sdhat))) * (- hattheta - Deltahatinf)

##     return(max(-nuhat_star + xiee))
##   }

##   return(Astarinf_PR_Deltas_obj)
## }

subvec <- function(theta, j, neg = FALSE) {
  ifelse(neg, -theta[j], theta[j])
}


subvec_extract <- function(thetas, p, sign) {
  if(!is.matrix(thetas)) thetas <- matrix(thetas, nrow = 1)
  sign * thetas[, p]
}
