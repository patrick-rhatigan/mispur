#' Calculate mbar, sigmahat, and mhat
#' (general case)
#'
#'
#' @param W Data matrix
#' @param m
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
gen_Deltahat_fun_g <- function(m) {
  Deltahat_fun_g <- function(W, theta) {
    mWtheta <- matrix(NA, nrow = n, ncol = k)
    for(i in 1:n) {
      mWtheta [i,] <- sapply(W[i,], m, theta=theta)
    }
    mbar <- colMeans(mWtheta)
    Sigmahat <- tcrossprod(mWtheta - mbar) / n
    sigmahat <- sqrt(diag(Sigmahat))
    max(-mbar/sigmahat)
  }
}


