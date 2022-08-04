#This file contains the calculations for the general case
spur_general <- function(W, m, Theta, theta0s, num_inti_val, alpha,
             alpha2, S, B, sdS, iota_sd, iota_q) {

  #Deltahatinf calculation:
  Deltahat_fun <- gen_Deltahat_fun_g(m=m)
  Deltahatinf <- rep(NA, ncol(Theta))
  for(i in 1:ncol(Theta)) {
    Deltahatinf[i] = optimize(Deltahat_fun, W=W,
                              lower = Theta[1,i],
                              upper = Theta[2,i])$objective
  }
  rhatinf <- pos(Deltahatinf)

}