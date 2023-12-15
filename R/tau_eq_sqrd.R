tau_eq_sqrd <- function(tau1sqrd,tau2sqrd){
  tau_eq_sqrd <- tau1sqrd*tau2sqrd/(tau1sqrd+tau2sqrd)
  return(pmax(0,tau_eq_sqrd, na.rm=T))
}
