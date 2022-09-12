
#' Spurious Precision Calculation
#'
#' @param W Matrix of the data (n by k)
#' @param m The Function of W and theta, if additively separable leave as NULL
#' @param additive Boolean indicating whether the function m
#'  is additively separable.
#'  i.e. it can be represented as m(W,theta) = g(W) + f(theta)
#' @param g The function g(W)
#' @param f The function f(theta)
#' @param Theta A matrix representing the parameter space of Theta.
#' The matrix should be 2 by k with the first row representing the lower bounds
#' of Theta and the second row the upper bounds.
#' @param theta0s
#' @param num_init_val Number of initial values to be used in each
#' bootstrap sample optimization. Choosing values above the default of one
#' may increase precision, but will also increase calculation time.
#' @param alpha Significance level of the SPUR2 test, default is 0.05
#' @param alpha2 Significance level used for the SPUR1 and GMS tests,
#' default is 0.045.
#'
#' alpha = alpha1 + alpha2 so alpha1 is computed by the function given the users
#' inputs of alpha and alpha2. alpha1 is the significance level used to construct
#' the CI for r_inf.
#' @param S
#' @param B Number of Bootstrap samples used. Default is 1000
#' @param sdS
#' @param c_tau Tau scaling parameter. The default value of Tau is sqrt(log(n)),
#' changing c_tau scales Tau up or down.
#' @param c_kappa Kappa scaling parameter. Identical to c_tau except for kappa.
#'  Kappa's default value is also sqrt(log(n))
#' @param iota_sd
#' @param iota_q

#'
#' @return
#' @export
#'
#' @examples
spur <- function(W, m = NULL,  additive = FALSE, g = NULL, f = NULL,
                 Theta, theta0s = NULL, num_init_val = 1,
                 alpha = .05, alpha2 = 0.045, S = S1, B = 1e03, sdS = 250,
                 c_tau = 1, c_kappa = 1, iota_sd = 1e-06, iota_q = 1e-06) {
  #Store important variables
  n <- nrow(W)
  k <- ncol(W)
  alpha1 <- alpha - alpha2

  # EGMS tuning parameters
  kappa <- c_kappa * sqrt(log(n))
  tau <- c_tau * sqrt(log(n))
  tau_DR <- 0

  #Sanity checks
  if((is.null(g)||is.null(f)) && additive){
    stop("If the moment function is additively separable you must provide
         its component functions, g=g(W) and f=f(theta)")
  }

  if(is.null(m) && !additive){
    stop("If the moment function is not additively separable you must provide
         the function m=m(W,theta)")
  }

  # Additively separable case and general case are delt with
  # in seperate files
  if(additive) spur_add_sep(W=W, g=g, f=f, Theta=Theta, theta0s=theta0s,
                            num_inti_val=num_inti_val, alpha=alpha,
                            alpha2=alpha2, S=S, B=B, sdS=sdS, iota_sd=iota_sd,
                            iota_q=iota_q)

  else spur_general(W=W, m=m, Theta=Theta, theta0s=theta0s,
                    num_inti_val=num_inti_val, alpha=alpha,
                    alpha2=alpha2, S=S, B=B, sdS=sdS, iota_sd=iota_sd,
                    iota_q=iota_q)
}


