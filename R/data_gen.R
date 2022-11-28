#' Data generation function for the upper lower bound model
#'
#' @param mu Mean vector of length k
#' @param n number of rows for W
#' @param k number of columns for W. Default = 2
#' @param Sigma The covariance matrix. The default is the identity matrix
#'
#' @return
#'
#'
#' @examples
gen_ulbd_data <- function(mu, n, k = 2, Sigma = diag(k)) {
  if(length(mu)!=k){
    stop("The mean vector must have a length equal to k")
  }

  if(nrow(Sigma)!=k || ncol(Sigma)!=k){
    stop("The covariance matrix must have dimensions k by k")
  }

  W <- mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma)

}
