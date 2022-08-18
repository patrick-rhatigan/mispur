gen_Deltahat_fun_g <- function(m) {
  Deltahat_fun_g <- function(W, theta) {

    mWtheta <- purrr::modify(W, m, theta=theta)
    mbar <- colMeans(mWtheta)
    Sigmahat <- tcrossprod(mWtheta - mbar) / n
    sigmahat <- sqrt(diag(Sigmahat))
    max(-mbar/sigmahat)
  }
}

gen_bootstrap_funs(m){
  mbar_star_fun <- function(W, theta) {
    mWtheta <- purrr::modify(W, m, theta=theta)
    return(colMeans(mWtheta))
  }
  sigma_star_fun <- function(W, theta) {
    mWtheta <- purrr::modify(W, m, theta=theta)
    mbar <- mbar_star_fun(W=W, theta=theta)
    Sigma_star <- tcrossprod(mWtheta - mbar) / n
    return(sqrt(diag(Sigmahat)))
  }
  mhat_star_fun <- function(W, theta) {
   return(mbar_star_fun(W, theta)/sigma_star_fun(W, theta))
  }
  return(list(mbar_star_fun = mbar_star_fun, sigma_star_fun = sigma_star_fun,
              mhat_star_fun = mhat_star_fun))
}
