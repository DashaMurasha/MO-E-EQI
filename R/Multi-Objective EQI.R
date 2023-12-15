#  Calculates the expectation of Objective 1 and Objective 2 (functions of X) improving on the Pareto front
mult_EQI = function(newdata,design_X, model_f1, model_f2, beta, tau_new, ConstraintInfo=NULL){
  
  # model_f1=model_dead; model_f2=model_cost
  
  k=length(design_X[1,])
  
  #Option=ModelInfo.Option
  
  #find points which satisfy constraint (if present)
  q1temp=q_mean(design_X, model_f1, beta)
  q2temp=q_mean(design_X, model_f2, beta)
  Xtemp=design_X
  
  if (!is.null(ConstraintInfo)) {
    for (i in 1:length(q1temp)){
      for (j in 1:length(ConstraintInfo$ConstraintLimits)){
        #We check if all the sampling points satisfy given constraints
        if (ConstraintInfo$y[i,j]>ConstraintInfo$ConstraintLimit[j]) {
          q1temp[i]=NaN
          q2temp[i]=NaN
        }
      }
    }
    Xtemp=Xtemp[!is.nan(q2temp),]
    q1temp=q1temp[!is.nan(q1temp)]
    q2temp=q2temp[!is.nan(q2temp)]
    
  }
  if (length(q1temp)==0){ 
    stop("No points satisfy constraints. Consider changing constraints.")
  }

  #find Pareto set
  find_pareto <- pareto_front(q1temp,q2temp,Xtemp)
  PX <- find_pareto$X
  Pq1 <- find_pareto$y1
  Pq2 <- find_pareto$y2
  Pnum <- length(find_pareto$y2)
  # b=order(q1temp)
  # PX <- Xtemp[b[1],1:k]
  # Pq1 <- q1temp[b[1]]
  # Pq2 <- q2temp[b[1]]
  # Pnum=1
  # for (i in 2:length(q1temp)){
  #   if (q2temp[b[i]]<=Pq2[Pnum]){
  #     Pnum=Pnum+1
  #     PX <- rbind(PX, Xtemp[b[i],1:k])
  #     Pq1 <- c(Pq1,q1temp[b[i]])
  #     Pq2 <- c(Pq2,q2temp[b[i]])
  #   }
  # }
  
  # prediction of each objective's quantile mean and sd value at x
  pred1=pred_Q(newdata, model_f1, beta, tau_new$tau1)
  pred2=pred_Q(newdata, model_f2, beta, tau_new$tau2)
  
  m_Q1 <- pred1$m_Q
  s_Q1 <- pred1$s_Q
  
  m_Q2<- pred2$m_Q
  s_Q2 <- pred2$s_Q
  
  # probability of improvement calculation
  PITerm1=pnorm((Pq1[1]-m_Q1)/s_Q1)
  
  PITerm3=(1-pnorm((Pq1[Pnum]-m_Q1)/s_Q1))*pnorm((Pq2[Pnum]-m_Q2)/s_Q2)
  
  if (Pnum>1) {
    PITerm2calc <- NULL
    for (I in 1:(length(Pq1)-1)){
      PITerm2calc <- rbind(PITerm2calc,(
        pnorm((Pq1[I+1]-m_Q1)/s_Q1)-
          pnorm((Pq1[I]-m_Q1)/s_Q1))*
          pnorm((Pq2[I+1]-m_Q2)/s_Q2)
      )
    }
    PITerm2=colSums(PITerm2calc)
    PI=(PITerm1+PITerm2+PITerm3)
  }
  else {
    PI=(PITerm1+PITerm3)
  }
  
  if (Option=='NegLogEQI'){
    # Qbar1 calculation
    Qbar1Term1=m_Q1*pnorm((Pq1[1]-m_Q1)/s_Q1)-
      s_Q1*dnorm((Pq1[1]-m_Q1)/s_Q1)
    
    Qbar1Term3=(m_Q1*pnorm((m_Q1-Pq1[Pnum])/s_Q1)+
                  s_Q1*dnorm((m_Q1-Pq1[Pnum])/s_Q1))*
      pnorm((Pq2[Pnum]-m_Q2)/s_Q2)
    
    if (Pnum>1){
      Qbar1Term2calc <- NULL
      for (I in 1:(length(Pq1)-1)){
        Qbar1Term2calc <- rbind(Qbar1Term2calc,
                                ((m_Q1*pnorm((Pq1[I+1]-m_Q1)/s_Q1)-
                                    s_Q1*dnorm((Pq1[I+1]-m_Q1)/s_Q1))-
                                   
                                   (m_Q1*pnorm((Pq1[I]-m_Q1)/s_Q1)-
                                      s_Q1*dnorm((Pq1[I]-m_Q1)/s_Q1))
                                )*
                                  pnorm((Pq2[I+1]-m_Q2)/s_Q2))
      }
      Qbar1Term2=colSums(Qbar1Term2calc)
      Qbar1=(Qbar1Term1+Qbar1Term2+Qbar1Term3)/PI
    }
    else{
      Qbar1=(Qbar1Term1+Qbar1Term3)/PI;
    }
    
    Qbar2Term1=m_Q2*pnorm((Pq2[Pnum]-m_Q2)/s_Q2)-
      s_Q2*dnorm((Pq2[Pnum]-m_Q2)/s_Q2)
    
    Qbar2Term3=(
      m_Q2*pnorm((m_Q2-Pq2[1])/s_Q2)+
        s_Q2*dnorm((m_Q2-Pq2[1])/s_Q2))*pnorm((Pq1[1]-m_Q1)/s_Q1)
    
    if (Pnum>1){
      Qbar2Term2calc <- NULL
      for (I in length(Pq2):2){
        Qbar2Term2calc <- rbind(Qbar2Term2calc,(
          (m_Q2*pnorm((Pq2[I-1]-m_Q2)/s_Q2)-
             s_Q2*dnorm((Pq2[I-1]-m_Q2)/s_Q2))-
            (m_Q2*pnorm((Pq2[I]-m_Q2)/s_Q2)-
               s_Q2*dnorm((Pq2[I]-m_Q2)/s_Q2)))*
            pnorm((Pq1[I-1]-m_Q1)/s_Q1)
        )
      }
      Qbar2Term2=colSums(Qbar2Term2calc);
      Qbar2=(Qbar2Term1+Qbar2Term2+Qbar2Term3)/PI
    }
    else{
      Qbar2=(Qbar2Term1+Qbar2Term3)/PI
    }
    # find closest point on front
    dist <- NULL
    for (i in 1:Pnum){
      dist=rbind(dist,sqrt((Qbar1-Pq1[i])^2+(Qbar2-Pq2[i])^2))
    }
    b <- apply(dist,2, min)
    # expected improvement calculation
    EQI=PI*b
    
    EQI[which(PI==0)]=0
    
    # low_sd_points <- set.intersection(which(predict.km(model_f1, newdata, type="UK")$sd<=epsilon_noise),
    #                                    which(predict.km(model_f2, newdata, type="UK")$sd<=epsilon_noise))
    # if (length(low_sd_points)!=0) EQI[low_sd_points]=0
  }
  if (Option=='NegLogEQI'){
    
    metric=-log10(EQI)
    
  }
  else{
    metric=-PI
  }
  
  return(list(metric=metric,Pq1=Pq1,Pq2=Pq2,PX=PX, s_Q1=s_Q1, s_Q2=s_Q2))
}



