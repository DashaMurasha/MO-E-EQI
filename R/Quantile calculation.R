#quantile calculator

#quantile calculator
q_mean <- function(x, model, beta){
  pred_f <- predict.km(model, x, type="UK")
  m_n <- pred_f$mean
  s <- pred_f$sd
  q <- m_n+qnorm(beta)*s
  return(q)
}
