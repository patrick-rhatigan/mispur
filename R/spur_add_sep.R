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

  # Calculate Deltahatinf
  init_vals <- t(Theta[1, ] + (Theta[2, ] - Theta[1, ]) *
                   t(randtoolbox::sobol(num_init_vals, dim = d)))

  print("Calculating Deltahat_inf...")

  start_delinf <- Sys.time()
  all_opt_results <- foreach(j = 1:nrow(init_vals)) %dopar% {
    print(paste0("Working on initial value #", j, "..."))
    init_val <- init_vals[j, ]

    const <- function(thetatil) {
      const_alt(thetatil, gWbar=gWbar, sigmahat=sigmahat, f=f)
    }

    gam_init <- max(-mhat_fun(init_val, gWbar=gWbar, sigmahat=sigmahat, f=f))

    start_time <- Sys.time()

    const_jac <- function(thetatil) {
      cbind(1, f(thetatil[-1])/sigmahat)
    }

    Deltahatinf_opt <-
      nloptr::slsqp(c(0, init_val), function(thetatil) thetatil[1],
                    hin = const,
                    hinjac = const_jac,
                    lower = c(-Inf, Theta[1, ]),
                    upper = c(Inf, Theta[2, ]),
                    control = list(xtol_rel = 1e-16, maxeval = 2000))

    end_time <- Sys.time()
    delinf_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    print(paste0("Done with initial value #", j, "..."))
    print(paste0("Opt val: ", Deltahatinf_opt$val, ", time:", delinf_time))
    return(Deltahatinf_opt)
   }

  end_delinf <- Sys.time()

  print(paste0("Time taken to compute Deltahatinf with ", num_init_vals,
               " initial values"))
  print(as.numeric(difftime(end_delinf, start_delinf, units = "secs")))

  mins <- sapply(1:length(all_opt_results),
                 function(j) all_opt_results[[j]]$val,
                 simplify = TRUE)

  min_ind <- which(min(mins) == mins)[1]
  Deltahatinf <- all_opt_results[[min_ind]]$val
  Deltahatinf_par <- all_opt_results[[min_ind]]$par[-1]
  rhatinf <- pos(Deltahatinf)

  # Omegahatplus
  mhat_n <- m_minus_mhat_fun_as(gW, gWbar, sigmahat)
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
  sqrt_ohplus_simZ <- t(mhat_mhat2 %*% matrix(rnorm(sdS * n), n, sdS) / sqrt(n))
  # Create Thetahat, Thetahat_min, Thetahat_min_DR functions
  # Thetahat_min and Thetahat_min_DR can be constructed using the same function

  ## gen_Thetahat <-
  ##   function(theta, rhatinf, gWbar, sigmahat, rho, tau, n) {
  Thetahat <- function(theta) {
    tau / sqrt(n) - max(neg(mhat_fun(theta, gWbar, sigmahat, f) + rhatinf))
  }
  ## }

  Thetahat_noneg <- function(theta) {
    tau / sqrt(n) + mhat_fun(theta, gWbar, sigmahat, f) + rhatinf
  }

  Thetahat_vec <- function(theta) {
    tau / sqrt(n) - neg(mhat_fun(theta, gWbar, sigmahat, f) + rhatinf)
  }

  ## gen_Thetahat_min <-
  ##   function(theta, Deltahatinf, gWbar, sigmahat, rho, tau, n) {
  Thetahat_min <- function(theta) {
    Deltahatinf + tau / sqrt(n) - Deltahat_fun_as(theta, gWbar, sigmahat, f)
  }

  Thetahat_min_noneg <- function(theta) {
    Deltahatinf + tau / sqrt(n) + mhat_fun_as(theta, gWbar, sigmahat, f)
  }
  ## }


  Thetahat_min_DR <- function(theta) {
    Deltahatinf + tau_DR / sqrt(n) - Deltahat_fun_as(theta, gWbar, sigmahat, f)
  }

  Thetahat_min_DR_noneg <- function(theta) {
    Deltahatinf + tau_DR / sqrt(n) + mhat_fun_as(theta, gWbar, sigmahat, f)
  }

  const_jac <- function(theta) {
    f(theta)/sigmahat
  }

  # Bootstrap sample
  set.seed(1001)
  b_index <- matrix(sample.int(n, size = n*B, replace = TRUE), nrow = B)

  # Quantities irrelavant with the alternative

  gWstarbars <- matrix(NA, nrow = k, ncol = B)
  sigmastars <- matrix(NA, nrow = k, ncol = B)


  for (b in 1:B) {
    gWstarbars[, b] <- calc_summ_stats_as(W[b_index[b,],], g)$gWbar
    sigmastars[,b] <- calc_summ_stats_as(W[b_index[b,],], g)$sigmahat
  }

  # Use only sdS samples for scale factors
  sd_gWstarbars <- gWstarbars[, 1:sdS]
  sd_sigmastars <- sigmastars[, 1:sdS]

  # Generate objective functions
  Astarinf_obj <-
    gen_Astarinf_obj_as(gWbar = gWbar, sigmahat = sigmahat,
                     rhatinf = rhatinf,
                     sd_gWstarbars = sd_gWstarbars,
                     sd_sigmastars = sd_sigmastars,
                     n = n, kappa = kappa, f = f,
                     iota_sd = iota_sd)

  Astarinf_Deltas_obj <-
    gen_Astarinf_Deltas_obj_as(gWbar = gWbar, sigmahat = sigmahat,
                            Deltahatinf = Deltahatinf,
                            sqrt_ohplus_simZ = sqrt_ohplus_simZ,
                            n = n, kappa = kappa, f = f,
                            iota_sd = iota_sd)

  Astarinf_DR_Deltas_obj <-
    gen_Astarinf_DR_Deltas_obj_as(gWbar = gWbar, sigmahat = sigmahat,
                               Deltahatinf = Deltahatinf,
                               sqrt_ohplus_simZ = sqrt_ohplus_simZ,
                               n = n, kappa = kappa, f = f,
                               iota_sd = iota_sd)

  init_vals <- rbind(Deltahatinf_par, init_vals)
  start_As <- Sys.time()
  all_As <- foreach (b = 1:B, .combine = rbind) %dopar% {
    ## print(paste0("b = ", b, ": ", paste0(b_index[b,], collapse = ",")))
    print(paste0("b = ", b))
    As_inf <- NULL
    As_inf_Delta <- NULL
    As_DR_Delta <- NULL

    for (j in 1) {

      ## for(j in 1:nrow(init_vals)) {
      init_val <- init_vals[j, ]

      if (Thetahat(init_val) >= 0) {
        ## print("As_inf")
        ## print(j)

        Astarinf_obj_alt <- function(theta) {
          Astarinf_obj_as(theta, gWstarbars[, b], sigmastars[, b])
        }

        time_Astarinf <- system.time(As_inf_opt <-
                                       nloptr::slsqp(init_val, Astarinf_obj_alt,
                                                     lower = Theta[1, ],
                                                     upper = Theta[2, ],
                                                     hin = Thetahat_noneg,
                                                     hinjac = const_jac)
        )

        print("As_inf results:")
        print(paste(As_inf_opt$val, time_Astarinf[3]))
        #[3] used here (and below with) system.time to get time elapsed

        As_inf <- min(As_inf, As_inf_opt$val)
      }

      if (Thetahat_min(init_val) >= 0) {
        ## print("As_inf_Delta")
        ## print(j)

        Astarinf_Deltas_obj_alt <- function(theta) {
          Astarinf_Deltas_obj_as(theta, gWstarbar = gWstarbars[, b],
                              sigmastar = sigmastars[, b])
        }

        time_Delta <- system.time(As_inf_Delta_opt <-
                                    nloptr::slsqp(init_val, Astarinf_Deltas_obj_alt,
                                                  lower = Theta[1, ],
                                                  upper = Theta[2, ],
                                                  hin = Thetahat_min_noneg,
                                                  hinjac = const_jac)
        )
        print("As_inf_Delta results:")
        print(paste(As_inf_Delta_opt$val, time_Delta[3]))

        As_inf_Delta <- min(As_inf_Delta, As_inf_Delta_opt$val)
      }

      if (Thetahat_min_DR(init_val) >= -1e-03) {
        ## print("As_DR_Delta")
        ## print(j)

        Astarinf_DR_Deltas_obj_alt <- function(theta) {
          Astarinf_DR_Deltas_obj_as(theta, gWstarbar = gWstarbars[, b],
                                 sigmastar = sigmastars[, b])
        }

        time_DR <- system.time(As_DR_Delta_opt <-
                                 nloptr::slsqp(init_val, Astarinf_DR_Deltas_obj_alt,
                                               lower = Theta[1, ],
                                               upper = Theta[2, ],
                                               hin = Thetahat_min_DR_noneg,
                                               hinjac = const_jac))


        print("As_DR_Delta results:")
        print(paste(As_DR_Delta_opt$val, time_DR[3]))

        As_DR_Delta <- min(As_DR_Delta, As_DR_Delta_opt$val)
      }
    }


    return(c(As_inf, As_inf_Delta, As_DR_Delta, time_Astarinf[3], time_Delta[3], time_DR[3]))
  }

  end_As <- Sys.time()
  print("Time taken to compute As:")
  print(end_As - start_As)
  print(colSums(all_As[, 4:6]) / 60 / 36)
  As_infs <- all_As[, 1]
  As_inf_Deltas <- all_As[, 2]
  As_DR_Deltas <- all_As[, 3]

  Deltahat_UB <- Deltahatinf + quantile(-As_inf_Deltas, 1 - alpha1) / sqrt(n) + iota_q
  rn <- max(Deltahat_UB, 0)


  CIs <- c(Deltahatinf - (quantile(As_DR_Deltas, 1 - alpha/2) + iota_q) / sqrt(n),
           Deltahatinf + (quantile(-As_inf_Deltas, 1 - alpha/2) + iota_q) / sqrt(n),
           Deltahatinf - (quantile(As_DR_Deltas, 1 - alpha) + iota_q) / sqrt(n),
           Deltahatinf + (quantile(-As_inf_Deltas, 1 - alpha) + iota_q)/sqrt(n),
           Deltahatinf)

  names(CIs) <- c("LB95", "UB95", "LB90", "UB90", "Deltahatinf")

  return(list(Deltahatinf = Deltahatinf,
              Deltahatinf_par = Deltahatinf_par,
              CIs = CIs, all_As = all_As[,1:3], b_index = b_index))

}

