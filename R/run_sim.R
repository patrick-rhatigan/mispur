#'Runs simulations for the generic upper bound/lower bound models
#'
#' Here, rather than providing r or length, the user provides the mean vector
#' and then r and length are calculated
#'
#'@param mu This gives the mean vector
#'@param k1 This gives the vector k1 which corresponds to
#' the moment "lower bounds"
#'@param k2 This gives the vector k2 which corresponds to
#' the moment "upper bounds"
#'@param ... These arguments are passed to \code{\link{momineq}}
#'@importFrom foreach %dopar%
#'

ulbd_sim <- function(Nsim = 1e03, n = 250, B = 1e03, sdS = 250,
                     k1, k2, mu, ncores) {

  k <- k1 + k2
  config <- paste0("mu=(",paste0(mu, collapse = ","),"),",
                   "_k1=", k1, "_k2=", k2)

  # Some sanity checks..sbatch -p pi_cowles ~/momineq/exclude/test_run_ws.sh

  if (length(mu) != k) {
    stop("mu must have length equal to ks[1]+ks[2]!!")
  }

  # Check misspecification
  mu1 <- mu[1:k1]
  mu2 <- mu[(k2 + 1):k]
  misspecified <- max(mu1) > min(mu2)
  Delta <- (max(mu1) - min(mu2))/2

  if (misspecified) {
    r <- (max(mu1) - min(mu2))/2
    length_id <- 0
  } else {
    r <- 0
    length_id <- min(mu2) - max(mu1)
  }

  # Set alternatives
  if (k > 2) {
    if ((r + length_id) >= 1) {
      theta0s <- seq(0, .6 * (r + length_id), by = 0.02)
    } else {
      theta0s <- seq(0, .6, by = 0.01)
    }
  } else if (k == 2) {
    if ((r + length_id) >= 1) {
      theta0s <- seq(0, .4 * (r + length_id), by = 0.01)
    } else {
      theta0s <- seq(0, .4, by = 0.01)
    }
  }

  # Setting seed for simulation
  set.seed(1001)
  doMC::registerDoMC(ncores)

  result <- foreach::foreach(j = 1:Nsim) %dopar% {

    # Generate data
    Sigma_sqrt <- diag(k)
    W <- Sigma_sqrt %*% matrix(rnorm(k*n), nrow = k) + mu

    result <-  momineq_sim(W = W, theta0s = theta0s, k1 = k1, k2 = k2,
                           sdS = sdS, B = B)
    rejects <- result$reject

    cat(paste0(paste0(result$CIs, collapse = " "), "\n"),
        file = paste0("output/", config, "_MI_CIs"),
        append = TRUE)

    return(rejects)
  }

  result <- Reduce("+", result) / Nsim
  attr(result, "mu") <- mu
  attr(result, "config") <- paste0("n=", n,
                                   "_Nsim=", Nsim,
                                   "_B=", B,
                                   "_", config,
                                   "_Delta=", Delta)

  return(result)
}
