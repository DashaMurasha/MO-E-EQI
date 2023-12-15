tau_new_func <- function(n_rep, tau,l) {
  tau1 <- sqrt(tau[1]^2/(n_rep))
  tau2 <- sqrt(tau[2]^2/(n_rep))
  return(data.frame(tau1=rep(tau1,l),tau2=rep(tau2,l)))
}
