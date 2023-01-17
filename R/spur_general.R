#This file contains the calculations for the general case

spur_general <- function(W, m, Theta, theta0s, num_inti_val, alpha,
             alpha2, S, B, sdS, iota_sd, iota_q) {
#r and n defined as nrows and ncols in parent function spur

  #Deltahatinf calculation:
  Deltahat_fun <- gen_Deltahat_fun_g(m=m)
  Deltahat <- rep(NA, ncol(Theta))
  for(i in 1:ncol(Theta)) {
    Deltahatinf[i] = optimize(Deltahat_fun, W=W,
                              lower = Theta[1,i],
                              upper = Theta[2,i])
  }
  Deltahat_values <- sapply(1:length(Deltahat),
                            function(j) Deltahat[[j]]$obj,
                            simplify = TRUE)
  min_ind <- which(min(Deltahat_values) == Deltahat_values)[1]
  Deltahatinf <- Deltahat_values[[min_ind]]
  Deltahatinf_par <- Deltahat[[min_ind]]$min[-1]#?
  rhatinf <- pos(Deltahatinf)

#Omegahatplus?

  #Bootstrapping:
  set.seed(1001)
  b_index <- matrix(sample.int(n, size = n*B, replace = TRUE), nrow = B)

  #defining the Astarinf objective function:
  Astarinf_obj <- gen_Astarinf_obj_g(W, m, n, kappa, rhatinf, iota_sd)

  for(b in 1:B) {
    Wstar <- W[b_index[b,],] #bootstrap sample data


  }
}
