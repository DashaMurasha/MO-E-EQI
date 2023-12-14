pareto_front <- function(y1_obs, y2_obs, design){
  k=length(design[1,])
  Xtemp=design
  #find Pareto set
  b=order(y1_obs)
  PX <- Xtemp[b[1],1:k]
  Y1 <- y1_obs[b[1]]
  Y2 <- y2_obs[b[1]]
  Pnum=1
  for (i in 2:length(y1_obs)){
    if (y2_obs[b[i]]<=Y2[Pnum]){
      Pnum=Pnum+1
      PX <- rbind(PX, Xtemp[b[i],1:k])
      Y1 <- c(Y1,y1_obs[b[i]])
      Y2 <- c(Y2,y2_obs[b[i]])
    }
  }
  return(list(y1=Y1, y2=Y2, X=PX))
}
