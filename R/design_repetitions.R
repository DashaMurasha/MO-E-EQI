design_repetitions <- function(x, design) {
  design_rep <- NULL
  for (i in 1:length(design[, 1])) {
    design_rep <- rbind(design_rep,
                        c(i, prodlim::row.match(design[i, ], x)))
  }
  design_rep <- na.exclude(design_rep)
  return(design_rep)
}
