mean_obs <- function(y1, y2, tau1sqrd, tau2sqrd) {
  if (tau1sqrd == 0 & tau2sqrd == 0)
    y_eq <- (y1 + y2) / 2
  else
    y_eq <- (y1 / tau1sqrd + y2 / tau2sqrd) / (1 / tau1sqrd + 1 / tau2sqrd)
  return(y_eq)
}
