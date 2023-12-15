pred_Q <- function(x, model, beta, tau_new) {
  pred_f <- predict.km(model, x, type="UK")
  m_n <- pred_f$mean
  s <- pred_f$sd
  #m_Q <- rep(NA,length(s))
  m_Q <- m_n+qnorm(beta)*sqrt((tau_new^2*s^2)/(s^2+tau_new^2))
  #m_Q[!(tau_new!=0)] <- m_n+qnorm(beta)*s^2
  s_Q <- sqrt((s^2)^2/(s^2+tau_new^2))
  return(list(m_Q=m_Q, s_Q=s_Q))
}
