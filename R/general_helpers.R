#' Generate the Deltahat function to be optimized
#'
#' @param m a function with inputs W and theta
#'  where theta is an individual entry in Theta
#'  (see spur general for implemenation of the function)
#' @return
#'
#'
#' @examples
gen_Deltahat_fun_g <- function(m) {
  Deltahat_fun_g <- function(W, theta) {
  #using purr here may creates challenges I (patrick)
  #can change this to not use modify if it causes issues
    mWtheta <- purrr::modify(W, m, theta=theta)
    mbar <- colMeans(mWtheta)
    Sigmahat <- tcrossprod(mWtheta - mbar) / n
    sigmahat <- sqrt(diag(Sigmahat))
    return(max(-mbar/sigmahat))
  }
  return(Deltahat_fun_g)
}



#' Generate Astar_inf objective function
#'
#' @param W original data matrix
#' @param m function m(W, theta)
#' @param n number of row in W
#' @param kappa parameter (will be passed in from spur function)
#' @param rhatinf minimal relaxation (will be passed in)
#' @param iota_sd buffer parameter (passed in)
#'
#' @return
#'
#'
#' @examples
gen_Astarinf_obj_g <- function(W, m, n, kappa, rhatinf, iota_sd) {
  #Wstar is the bootstrap sample data matrix
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


    # calculate sd2, sd3
    Vhatstar2 <- n * matrixStats::rowVars(mhat_star)
    sd2 <- sqrt(Vhatstar2)
    sd2[sd2 < iota_sd] <- iota_sd

    negmhatstars <- neg(mhat_star)
    maxmhatss <- matrixStats::colMaxs(negmhatstars)
    Vhatstar3 <- n * matrixStats::colVars(t(negmhatstars) - maxmhatss)
    sd3 <- sqrt(Vhatstar3)
    sd3[sd3 < iota_sd] <- iota_sd

    chistartheta <- ifelse(nuhat_star >= 0,
                           chi(nuhat_star, sqrt(n) * mhat_star - (kappa * sd2)),
                           chi(nuhat_star, sqrt(n) * mhat_star + (kappa * sd2)))

    bhattheta <- sqrt(n) * (neg(mhat_star) - rhatinf) - kappa * sd3
    xibtheta <- sqrt(n) * (neg(mhat_star) - rhatinf) / (kappa * sd3)
    phitheta <- phi(xibtheta)

    k <- length(mhat_star)
    temp <- rep(NA, k)
    for (j1 in 1:k) {
      temp[j1] <- max(chistartheta +
                        ifelse((1:k) == j1, phitheta, bhattheta))
    }

    jhat <- which(neg(mhat_star) >= max(neg(mhat_star)) - sd3 * kappa / sqrt(n))
    ret <- min(temp[jhat])

    if (ret == Inf) {ret <- 10^10}
    return(ret)
  }

  return(Astarinf_obj)
}
