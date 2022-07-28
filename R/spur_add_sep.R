#This file contains the calculations for the additively seperable case

#' Spurious Precision Calculation for Additively Separable Case
#'
#' @param ... All params are passed down from the spur function call.
#' @return
#' @export
#'
#' @examples
spur_add_sep <- function(W, g, f, alpha, alpha2,
                         S, B ,Theta, theta0s, rho){
  #summary stats
  gW <- matrix(NA,nrow = n, ncol = k)
  for(i in 1:n) {
    gW [i,] <- sapply(W[i,], g)
  }
  gWbar = colMeans(gW)

  Sigmahat <- tcrossprod(gW - gWbar) / n
  sigmahat <- sqrt(diag(Sigmahat))


  init_val <- rowMeans(Theta)

  Deltahatinf_opt <- optim(init_val, Deltahat_fun,
                           lower = Theta[, 1],
                           upper = Theta[, 2],
                           method = "L-BFGS-B",
                           gWbar = gWbar, sigmahat = sigmahat, rho = rho)

  Deltahatinf <- Deltahatinf_opt$val

  rhatinf <- pos(Deltahatinf)

  # Omegahatplus
  mhat_n <- m_minus_mhat_fun(gW, gWbar, sigmahat)
  mhat_n2 <- mhat_n^2 - 1
  mhat_mhat2 <- rbind(mhat_n, mhat_n2)
  omegahatplus <- tcrossprod(mhat_mhat2) / n
  sqrt_ohplus_simZ <- mvnfast::rmvn(n = sdS,
                                    mu = rep(0, 2 * k),
                                    sigma = omegahatplus)

  mhat <- (gW - gWbar) / sigmahat
  mhat2 <- mhat^2 - 1
  mhat_mhat2 <- rbind(mhat, mhat2)

  omegahatplus <- tcrossprod(mhat_mhat2) / n
  ## eig_omegahatplus <- eigen(omegahatplus)
  ## eigvec_ohp <- eig_omegahatplus$vectors
  ## eigval_ohp <- eig_omegahatplus$values
  ## sqrt_ohplus <- eigvec_ohp %*% diag(sqrt(eigval_ohp)) %*% t(eigvec_ohp)
  ## NsimZ <- B
  ## simZ <- matrix(rnorm(2 * k * NsimZ), ncol = 2 * k)
  ## sqrt_ohplus_simZ <-  simZ %*% sqrt_ohplus # NsimZ x 2k
  ## sqrt_ohplus_simZ <- mvtnorm::rmvnorm(n = B, sigma = omegahatplus) # B x 2k
  sqrt_ohplus_simZ <- mvnfast::rmvn(n = B, mu = rep(0, 2 * k), sigma = omegahatplus)


  # Thetahat_min
  Thetahat_min <- matrix(NA, nrow = k, ncol = 2)
  for (i in 1:k) {
    bound <- sigmahat * (tau/sqrt(n) + Deltahatinf)
    if (i <= k1) {
      Thetahat_min[i,] <- c(-bound[i] + Wbar[i], Theta[2])
    } else {
      Thetahat_min[i,] <- c(Theta[1], Wbar[i] + bound[i])
    }
  }

  Thetahat_min <- c(max(Thetahat_min[,1]), min(Thetahat_min[,2]))

  # Thetahat_min_DR
  Thetahat_min_DR <- matrix(NA, nrow = k, ncol = 2)
  for (i in 1:k) {
    bound <- sigmahat * (tau_DR/sqrt(n) + Deltahatinf)
    if (i <= k1) {
      Thetahat_min_DR[i,] <- c(-bound[i] + Wbar[i], Theta[2])
    } else {
      Thetahat_min_DR[i,] <- c(Theta[1], Wbar[i] + bound[i])
    }
  }

  Thetahat_min_DR <- c(max(Thetahat_min_DR[,1]), min(Thetahat_min_DR[,2]))

  # Thetahat
  Thetahat <- matrix(NA, nrow = k, ncol = 2)
  for (i in 1:k) {
    bound <- sigmahat * (tau/sqrt(n) + rhatinf)
    if (i <= k1) {
      Thetahat[i,] <- c(-bound[i] + Wbar[i], Theta[2])
    } else {
      Thetahat[i,] <- c(Theta[1], Wbar[i] + bound[i])
    }
  }

  Thetahat <- c(max(Thetahat[,1]), min(Thetahat[,2]))

  b_index <- matrix(sample.int(n, size = n*B, replace = TRUE), nrow = B)


  # Quantities irrelavant with the alternative
  Wstarbars <- matrix(NA, nrow = k, ncol = B)
  sigmastars <- matrix(NA, nrow = k, ncol = B)
  Deltahatstars <- numeric(B)


  for (b in 1:B) {
    Wstar <- W[,b_index[b,]]
    Wstarbars[, b] <- Wstarbar <- rowMeans(Wstar)
    Sigmastar <- tcrossprod(Wstar - Wstarbar) / n
    sigmastars[,b] <- sqrt(diag(Sigmastar))
  }

  nu_ms <- sqrt(n) * (Wbar - Wstarbars) * c(rep(1,k1),rep(-1,k2)) /
    sigmahat
  nu_ms <- sqrt(n) * (Wbar - Wstarbars) * c(rep(1,k1),rep(-1,k2)) /
    sigmahat
  # Store A^\ast, \inf and A^\ast,\sup
  As_infs <- numeric(B)
  As_inf_Deltas <- numeric(B)

  As_DR_Deltas <- numeric(B)
  As_PR_Deltas <- numeric(B)
  As_MR_Deltas <- numeric(B)

  for (b in 1:B) {
    As_infs[b] <- optimize(Astarinf_obj, Thetahat, Wbar = Wbar,
                           Sigmahat = Sigmahat, rhatinf = rhatinf,
                           sigmastars = sigmastars, Wstarbars = Wstarbars,
                           n = n, b = b, B = B, kappa = kappa,
                           k1 = k1, k2 = k2)$obj

    As_inf_Deltas[b] <- optimize(Astarinf_Deltas_obj, Thetahat_min, Wbar = Wbar,
                                 Sigmahat = Sigmahat, Deltahatinf = Deltahatinf,
                                 sigmastars = sigmastars, Wstarbars = Wstarbars,
                                 n = n, b = b, B = B, kappa = kappa,
                                 sqrt_ohplus_simZ = sqrt_ohplus_simZ,
                                 k1 = k1, k2 = k2)$obj

    A_DR_Delta <- optimize(Astarinf_DR_Deltas_obj, Thetahat_min_DR, Wbar = Wbar,
                           Sigmahat = Sigmahat, Deltahatinf = Deltahatinf,
                           sigmastars = sigmastars, Wstarbars = Wstarbars,
                           n = n, b = b, B = B, kappa = kappa,
                           sqrt_ohplus_simZ = sqrt_ohplus_simZ,
                           k1 = k1, k2 = k2)$obj

    A_PR_Delta <- optimize(Astarinf_PR_Deltas_obj, Theta, Wbar = Wbar,
                           Sigmahat = Sigmahat, Deltahatinf = Deltahatinf,
                           sigmastars = sigmastars, Wstarbars = Wstarbars,
                           n = n, b = b, B = B, kappa = kappa,
                           sqrt_ohplus_simZ = sqrt_ohplus_simZ,
                           k1 = k1, k2 = k2)$obj

    As_DR_Deltas[b] <- A_DR_Delta
    As_PR_Deltas[b] <- A_PR_Delta
    As_MR_Deltas[b] <- min(A_DR_Delta, A_PR_Delta)
  }

  Deltahat_UB <- Deltahatinf + quantile(-As_inf_Deltas, 1-alpha1)/sqrt(n)
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
           Deltahatinf + quantile(-As_inf_Deltas, 1 - alpha)/sqrt(n))
  names(CIs) <- c("LB95", "LB95_DR", "LB95_PR", "UB95",
                  "LB90", "LB90_DR", "LB90_PR", "UB90")

  reject <- as.data.frame(matrix(NA, nrow = length(theta0s), ncol = 3))
  names(reject) <- c("SPUR1", "SPUR2", "GMS")

  row.names(reject) <- theta0s

  for (i in 1:length(theta0s)) {

    theta0 <- theta0s[i]

    # mhat and mbar
    mhattheta0 <- mhat(theta0, Wbar, Sigmahat, k1 = k1, k2 = k2)
    mbartheta0 <- mhattheta0 * sigmahat

    # bootstrap counterparts
    mbarstars <- (theta0 - Wstarbars) * c(rep(1,k1),rep(-1,k2))
    mhatstarstars <- mhat(theta0, Wstarbars, sigmas = sigmastars,
                          k1 = k1, k2 = k2)
    nuhat_star_theta0 <- sqrt(n) * (mhatstarstars - mhattheta0)

    # sd1
    negmhatstars <- neg(mhatstarstars)
    maxmhatss <- apply(negmhatstars, 2, function(mhat) max(mhat))

    Vstar1summand <- t(t(mhatstarstars) + maxmhatss)
    Vhatstar1 <- n * ((B - 1)/B) *
      sapply(1:k, function(l) var(Vstar1summand[l,]))
    sd1 <- sqrt(pmax(Vhatstar1, rep(1,k)))

    # GMS part for SPUR1
    # Add phitheta0 to nuhat_theta0
    rhattheta0 <- rhat(theta0, Wbar, Sigmahat, k1 = k1, k2 = k2)
    xitheta0 <- (sqrt(n)/kappa)*(mhattheta0 + rhattheta0)
    phitheta0_1EGMS <- phi((1/sd1) * xitheta0)
    Tstars <- nuhat_star_theta0 + phitheta0_1EGMS

    # GMS part for std GMS
    phitheta0_GMS <- phi((sqrt(n)/kappa)*mhattheta0)

    # Test stats
    Stheta0_1step <- S(sqrt(n)*(mhattheta0 + rhatinf))
    Stheta0_GMS <- S(sqrt(n) * mhattheta0)

    # Store bootstrapped test stats
    Ss_GMS <-  sapply(1:B, function(b) S(nu_ms[,b] + phitheta0_GMS))
    Ss_SPUR1 <- sapply(1:B, function(b) S(Tstars[,b] + As_infs[b]))

    # Rejections
    # SPUR1
    reject[i,1] <- 1 * c(Stheta0_1step > quantile(Ss_SPUR1, 1 - alpha))
    # SPUR2
    reject_SPUR_alpha2 <- 1 * c(Stheta0_1step > quantile(Ss_SPUR1, 1 - alpha2))
    reject_GMS_alpha2 <- 1 * (Stheta0_GMS > quantile(Ss_GMS, 1-alpha2))
    reject[i,2] <- ifelse(rn == 0,
                          reject_GMS_alpha2,
                          min(reject_GMS_alpha2, reject_SPUR_alpha2))
    # Standard GMS
    reject[i,3] <-  1 * c(Stheta0_GMS > quantile(Ss_GMS, 1 - alpha))
  }

  return(list(reject = reject, CIs = CIs))
}
