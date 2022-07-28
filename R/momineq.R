#'Implementing SPUR1, SPUR2, and GMS
#'
#'This is the function that runs the inference method introduced by Andrews
#'and Kwon (2019) but the focus is to replicate the simulations given in the
#'paper. That is, this function is similar to \code{\link{momineq}} but always
#'runs SPUR1, SPUR2, and GMS (when r = 0). The following function is for the
#' lower/upper bound model, and is specialized for the simulation.
#'
#'@export

momineq_sim <- function(W, k1, k2, alpha = .05, alpha2 = 0.045, S = S1,
                        B = 1e03, sdS = 250,
                        Theta = c(-20,20), theta0s, r = NULL,
                        length_id = NULL) {

  k <- k1 + k2
  n <- ncol(W)
  alpha1 <- alpha - alpha2

  # Generate some commonly used functions
  mhat_fun <- gen_mhat_fun(k1, k2)
  m_minus_mhat_fun <- gen_m_minus_mhat_fun(k1, k2)
  Deltahat_fun <- gen_Deltahat_fun(k1, k2)

  # Tuning parameters
  kappa <- sqrt(log(n))
  tau <- sqrt(log(n))
  tau_DR <- 0

  # Some summary stats
  k <- k1 + k2
  Wbar <- rowMeans(W)
  Sigmahat <- tcrossprod(W - Wbar) / n
  sigmahat <- sqrt(diag(Sigmahat))

  # Deltahatinf, rhatinf
  Deltahatinf <- optimize(Deltahat_fun,
                          Theta, Wbar = Wbar, sigmahat = sigmahat)$obj
  rhatinf <- pos(Deltahatinf)

  # Omegahatplus
  mhat_n <- m_minus_mhat_fun(W, Wbar, sigmahat)
  mhat_n2 <- mhat_n^2 - 1
  mhat_mhat2 <- rbind(mhat_n, mhat_n2)
  omegahatplus <- tcrossprod(mhat_mhat2) / n
  sqrt_ohplus_simZ <- mvnfast::rmvn(n = sdS,
                                    mu = rep(0, 2 * k), sigma = omegahatplus)


  # Thetahat_min, Thetahat
  gen_Thetahat <- function(margin) {
    Thetahat <- matrix(NA, nrow = k, ncol = 2)
    for (i in 1:k) {
      bound <- sigmahat * margin
      if (i <= k1) {
        Thetahat[i,] <- c(-bound[i] + Wbar[i], Theta[2])
      } else {
        Thetahat[i,] <- c(Theta[1], Wbar[i] + bound[i])
      }
    }

    Thetahat <- c(max(Thetahat[,1]), min(Thetahat[,2]))

    return(Thetahat)
  }

  Thetahat_min <- gen_Thetahat(tau/sqrt(n) + Deltahatinf)
  Thetahat_min_DR <- gen_Thetahat(tau_DR/sqrt(n) + Deltahatinf)
  Thetahat <- gen_Thetahat(tau/sqrt(n) + rhatinf)

  # Bootstrap sample
  b_index <- matrix(sample.int(n, size = n*B, replace = TRUE), nrow = B)

  # Quantities irrelavant with the alternative
  Wstarbars <- matrix(NA, nrow = k, ncol = B)
  sigmastars <- matrix(NA, nrow = k, ncol = B)


  for (b in 1:B) {
    Wstar <- W[,b_index[b,]]
    Wstarbars[, b] <- Wstarbar <- rowMeans(Wstar)
    Sigmastar <- tcrossprod(Wstar - Wstarbar) / n
    sigmastars[,b] <- sqrt(diag(Sigmastar))
  }

  # Use only sdS samples for scale factors
  sd_Wstarbars <- Wstarbars[, 1:sdS]
  sd_sigmastars <- sigmastars[, 1:sdS]

  # Store A^\ast, \inf and A^\ast,\sup
  As_infs <- numeric(B)
  As_inf_Deltas <- numeric(B)

  As_DR_Deltas <- numeric(B)
  As_PR_Deltas <- numeric(B)
  As_MR_Deltas <- numeric(B)

  # Generate objective functions
  Astarinf_obj <-
    gen_Astarinf_obj(Wbar = Wbar, sigmahat = sigmahat, mhat_fun = mhat_fun,
                     rhatinf = rhatinf, sd_Wstarbars = sd_Wstarbars,
                     sd_sigmastars = sd_sigmastars, n = n, kappa = kappa)
  Astarinf_Deltas_obj <-
    gen_Astarinf_Deltas_obj(Wbar = Wbar, sigmahat = sigmahat,
                            mhat_fun = mhat_fun, Deltahatinf = Deltahatinf,
                            sqrt_ohplus_simZ = sqrt_ohplus_simZ,
                            n = n, kappa = kappa)
  Astarinf_DR_Deltas_obj <-
    gen_Astarinf_DR_Deltas_obj(Wbar = Wbar, sigmahat = sigmahat,
                               mhat_fun = mhat_fun, Deltahatinf = Deltahatinf,
                               sqrt_ohplus_simZ = sqrt_ohplus_simZ,
                               n = n, kappa = kappa)

  Astarinf_PR_Deltas_obj <-
    gen_Astarinf_PR_Deltas_obj(Wbar = Wbar, sigmahat = sigmahat,
                               mhat_fun = mhat_fun, Deltahatinf = Deltahatinf,
                               sqrt_ohplus_simZ = sqrt_ohplus_simZ,
                               n = n, kappa = kappa)

  for (b in 1:B) {
    As_infs[b] <- optimize(Astarinf_obj, Thetahat, Wstarbar = Wstarbars[, b],
                           sigmastar = sigmastars[, b])$obj

    As_inf_Deltas[b] <- optimize(Astarinf_Deltas_obj, Thetahat_min,
                                 Wstarbar = Wstarbars[, b],
                                 sigmastar = sigmastars[, b])$obj

    A_DR_Delta <- optimize(Astarinf_DR_Deltas_obj, Thetahat_min_DR,
                           Wstarbar = Wstarbars[, b],
                           sigmastar = sigmastars[, b])$obj

    A_PR_Delta <- optimize(Astarinf_PR_Deltas_obj, Theta,
                           Wstarbar = Wstarbars[, b],
                           sigmastar = sigmastars[, b])$obj

    As_DR_Deltas[b] <- A_DR_Delta
    As_PR_Deltas[b] <- A_PR_Delta
    As_MR_Deltas[b] <- min(A_DR_Delta, A_PR_Delta)
  }

  Deltahat_UB <- Deltahatinf + quantile(-As_inf_Deltas, 1 - alpha1)/sqrt(n)
  rn <- max(Deltahat_UB, 0)

  ## Deltahat_LB_DR <- Deltahatinf - quantile(As_DR_Deltas, 1-alpha)/sqrt(n)
  ## Deltahat_LB_PR <- Deltahatinf - quantile(As_PR_Deltas, 1-alpha)/sqrt(n)
  ## Deltahat_LB <- Deltahatinf - quantile(As_MR_Deltas, 1-alpha)/sqrt(n)

  CIs <- c(Deltahatinf - quantile(As_MR_Deltas, 1 - alpha/2) / sqrt(n),
           Deltahatinf - quantile(As_DR_Deltas, 1 - alpha/2) / sqrt(n),
           Deltahatinf - quantile(As_PR_Deltas, 1 - alpha/2) / sqrt(n),
           Deltahatinf + quantile(-As_inf_Deltas, 1 - alpha/2) / sqrt(n),
           Deltahatinf - quantile(As_MR_Deltas, 1 - alpha) / sqrt(n),
           Deltahatinf - quantile(As_DR_Deltas, 1 - alpha) / sqrt(n),
           Deltahatinf - quantile(As_PR_Deltas, 1 - alpha) / sqrt(n),
           Deltahatinf + quantile(-As_inf_Deltas, 1 - alpha)/sqrt(n),
           Deltahatinf)
  names(CIs) <- c("LB95", "LB95_DR", "LB95_PR", "UB95",
                  "LB90", "LB90_DR", "LB90_PR", "UB90",
                  "Deltahatinf")

  reject <- as.data.frame(matrix(NA, nrow = length(theta0s), ncol = 3))
  names(reject) <- c("SPUR1", "SPUR2", "GMS")

  row.names(reject) <- theta0s

  for (i in 1:length(theta0s)) {

    theta0 <- theta0s[i]

    # mhat and mbar
    mhattheta0 <- mhat_fun(theta0, Wbar, sigmahat)
    mhatstarstars <- mhat_fun(theta0, Wstarbars, sigmastars)
    nuhat_stars_theta0 <- sqrt(n) * (mhatstarstars - mhattheta0)

    # sd1
    negmhatstars <- neg(mhatstarstars[,1:sdS])
    maxmhatss <- matrixStats::colMaxs(negmhatstars)

    Vhatstar1 <- n * matrixStats::colVars(t(mhatstarstars[,1:sdS]) + maxmhatss)
    sd1 <- sqrt(pmax1(Vhatstar1))

    # GMS part for SPUR1
    # Add phitheta0 to nuhat_theta0
    rhattheta0 <- max(neg(mhat_fun(theta0, Wbar, sigmahat)))
    xitheta0 <- (sqrt(n)/kappa)*(mhattheta0 + rhattheta0)
    phitheta0_1EGMS <- phi((1/sd1) * xitheta0)
    Tstars <- nuhat_stars_theta0 + phitheta0_1EGMS

    # Std GMS ingredients
    nu_ms <- sqrt(n) * m_minus_mhat_fun(Wstarbars, Wbar, sigmahat)
    phitheta0_GMS <- phi((sqrt(n)/kappa) * mhattheta0)

    # Test stats
    Stheta0_SPUR1 <- S(sqrt(n) * (mhattheta0 + rhatinf))
    Stheta0_GMS <- S(sqrt(n) * mhattheta0)

    # Store bootstrapped test stats
    if (identical(S, S1)) {
      Ss_GMS <- colSums(neg(nu_ms + phitheta0_GMS)^2)
      Ss_SPUR1 <- rowSums(neg(t(Tstars) + As_infs)^2)
    } else if (identical(S, S4)) {
      Ss_GMS <- matrixStats::colMaxs(neg(nu_ms + phitheta0_GMS))
      Ss_SPUR1 <- matrixStats::rowMaxs(neg(t(Tstars) + As_infs))
    } else {
      Ss_GMS <-  sapply(1:B, function(b) S(nu_ms[,b] + phitheta0_GMS))
      Ss_SPUR1 <- sapply(1:B, function(b) S(Tstars[,b] + As_infs[b]))
    }

    # Rejections
    # SPUR1
    reject[i,1] <- 1 * c(Stheta0_SPUR1 > quantile(Ss_SPUR1, 1 - alpha))
    # SPUR2
    reject_SPUR_alpha2 <- 1 * (Stheta0_SPUR1 > quantile(Ss_SPUR1, 1 - alpha2))
    reject_GMS_alpha2 <- 1 * (Stheta0_GMS > quantile(Ss_GMS, 1 - alpha2))
    reject[i,2] <- ifelse(rn == 0,
                          reject_GMS_alpha2,
                          min(reject_GMS_alpha2, reject_SPUR_alpha2))
    # Standard GMS
    reject[i,3] <-  1 * c(Stheta0_GMS > quantile(Ss_GMS, 1 - alpha))
  }

  return(list(reject = reject, CIs = CIs))
}



#######

#' Spurious Precision Calculation
#'
#' @param W Matrix of the data (n*k)
#' @param m The Function of W and theta, if additively separable leave as NULL
#' @param additive Boolean indicating whether the function m
#'  is additively separable.
#'  i.e. it can be represented as m(W,theta) = g(W) + f(theta)
#' @param g The function g(W)
#' @param f The function f(theta)
#' @param alpha
#' @param alpha2
#' @param S
#' @param B
#' @param Theta
#' @param theta0s
#' @param rho
#'
#' @return
#' @export
#'
#' @examples
spur <- function(W, m = NULL,  additive = FALSE, g = NULL, f = NULL,
                 alpha = .05, alpha2 = 0.045, S = S1, B = 1e03,
                 Theta, theta0s, rho = NULL) {
  #Store important variables
  n <- nrow(W)
  k <- ncol(W)
  alpha1 <- alpha - alpha2


  # EGMS tuning parameter
  kappa <- sqrt(log(n))
  tau <- sqrt(log(n))
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
  if(additive) spur_add_sep(W=W, g=g, f=f, alpha=alpha, alpha2=alpha2,
                            S=S, B=B,Theta=Theta, theta0s=theta0s,
                            rho=rho)

  else spur_general(W=W, m=m, alpha=alpha, alpha2=alpha2,
                    S=S, B=B,Theta=Theta, theta0s=theta0s,
                    rho=rho)
}




