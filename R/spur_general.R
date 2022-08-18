#This file contains the calculations for the general case

spur_general <- function(W, m, Theta, theta0s, num_inti_val, alpha,
             alpha2, S, B, sdS, iota_sd, iota_q) {

  #Summary stats:
  summ_stats_obj <- calc_summ_stats_g(W=W, m=m, theta=Theta)
  mbar <- summ_stats_obj$mbar
  sigmahat <- summ_stats_obj$sigmahat
  mhat <- summ_stats_obj$mhat

  #Deltahatinf calculation:
  Deltahat_fun <- gen_Deltahat_fun_g(m=m)
  Deltahatinf <- rep(NA, ncol(Theta))
  for(i in 1:ncol(Theta)) {
    Deltahatinf[i] = optimize(Deltahat_fun, W=W,
                              lower = Theta[1,i],
                              upper = Theta[2,i])$objective
  }
  rhatinf <- pos(Deltahatinf)

#Omegahatplus?

  #Bootstrapping:
  set.seed(1001)
  b_index <- matrix(sample.int(n, size = n*B, replace = TRUE), nrow = B)
  for(b in 1:B) {
    Wstar <- W[b_index[b,],] #bootstrap sample data
    bootstrap_funs <- gen_bootstrap_funs(m=m)
    mbar_star_fun <- bootstrap_funs$mbar_star_fun
    sigma_star_fun <- bootstrap_funs$sigma_star_fun
    mhat_star_fun <- bootstrap_funs$mhat_star_fun

  }
}
