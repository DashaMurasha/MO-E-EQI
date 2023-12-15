---
title: "MO-E-EQI"
output: html_document
date: "2023-12-15"
---

Supplemental materials for (!!!!!!!!!!!!!!! ref !!!!!!!!!!!!!!!!!!)

Authors: D. Semochkina, A.I.J. Forester, D.C. Woods

Corresponding Author: D. Semochkina (d.semochkina@soton.ac.uk), School of Mathematical Sciences, University of Southampton, UK 

This README.md file is a copy of an Rmarkdown file 'MO-E-EQI.Rmd' available with the package.

This repository provides R package to caclulate multi-objective Euclidian expected quantile improvement (MO-E-EQI)  presented in the manuscript. 

First install our package MOEEQI from git. We need the standard package 'devtools' to add our package off git.

```{r}
install.packages("devtools")
```

This is the standard way to import an R package into the current session.

```{r}
library("devtools")
```

Now we need to build our package MOEEQI from git.
```{r}
install_github("DashaMurasha/MO-E-EQI")
```

Note that the above instructions should only need running once in order to install our package.  After which we can just run:
```{r}
library("MOEEQI")
```

Next, we move on to the example accompanying the paper (!!!!!!!!!!!!!!! ref !!!!!!!!!!!!!!!!!!)

We first set the level of noise and define out funstions 
```{r}
# Set a
a=0.25
# Test functions
f1 <- function(x,theta){
  x1 <- unlist(x[,1])
  x2 <- unlist(x[,2])
  theta1 <- unlist(theta[,1])
  theta2 <- unlist(theta[,2])
  1-sin(x1)+x2/10+a*cos(theta1)+theta2/10
}
f2 <- function(x,theta){
  x1 <- unlist(x[,1])
  x2 <- unlist(x[,2])
  theta1 <- unlist(theta[,1])
  theta2 <- unlist(theta[,2])
  1-cos(x1)+x2/3+a*sin(theta1)+theta2/3
  
}
```
Choose the number of repetitions of each model run
```{r}
# Number of repetitions of each observation (model run)
MC_sample_size <- 100
```
Choose the number of steps of the MO-E-EsQI sequential design loop
```{r}
# Computational budget (steps in the loop)
Nsteps <- 9
```
Select initial design points.
```{r}
# Input parameters ranges
x_c_1_range <- c(0, pi/2)
x_c_2_range <-  c(0, 1)

# Number of original design points
n_sample_points <- 5

# Generate original design (maximin Latin hypercube)
design_X <- MaxPro::MaxProLHD( n = n_sample_points, p = 2, itermax = 20 )
design_X <- MaxPro::MaxPro( InitialDesign = design_X$Design, iteration = 10 )$Design

design_X <- as.matrix(t(t(design_X)*c(diff(x_c_1_range), diff(x_c_2_range))+
                          c(x_c_1_range[1],x_c_2_range[1])))
# Make it a data.frame
orig_design_X <- data.frame(x=design_X)
```

Choose the metric option. Currently two options available,  -log(EQI) ('NegLogEQI') or -EQI ('NegEQI')
```{r}
Option <- 'NegLogEQI'
```

select quantile level (see EQI (!!!!!!!!!!!!!!! ref !!!!!!!!!!!!!!!!!!) for details)
```{r}
beta <- .8
```

Define some objects to store results
```{r}
y1_orig <- y2_orig <- epsilons_orig <- NULL
y1_all <- y2_all <- NULL
var1 <- var2 <- NULL
y_plot <- NULL
```
Define the noise level for the second environmental variable
```{r}
mu_2 <- 0
sigma_2 <- 0.5
```

Run the model at the original design points
```{r}
for (i in 1:n_sample_points) {
  # data_points <- expand.grid(theta=orig_design_X$theta[i],x=rnorm(MC_sample_size,0.4,1))
  data_points <- cbind(x <- matrix(rep(orig_design_X[i,],each=MC_sample_size),MC_sample_size,2),
                       theta.1=runif(MC_sample_size,-pi,pi),
                       theta.2=rnorm(MC_sample_size,mu_2,sigma_2))
  
  y1_orig <- c(y1_orig, mean(f1(data_points[,1:2],data_points[,3:4])))
  y2_orig <- c(y2_orig, mean(f2(data_points[,1:2],data_points[,3:4])))
  y1_all <- rbind(y1_all, unlist(c(f1(data_points[,1:2],data_points[,3:4]), orig_design_X[i,])))
  y2_all <- rbind(y2_all, unlist(c(f2(data_points[,1:2],data_points[,3:4]), orig_design_X[i,])))
  var1 <- c(var1,var(f1(data_points[,1:2],data_points[,3:4])))
  var2 <- c(var2,var(f2(data_points[,1:2],data_points[,3:4])))
}
```
Provide noise sd for both objectives.
```{r}
noise_sd <- sqrt(c(mean(var1),mean(var2)))
noise_orig_design_sd <- cbind(apply(y1_all[,1:MC_sample_size],1,sd),
                              apply(y2_all[,1:MC_sample_size],1,sd))

```
Write original design and the starting point of the final design.
```{r}
design_X <- orig_design_X
```
Calculate future noise. Note that tau is standard deviation (not variance).
```{r}
tau_new <- tau_new_func(MC_sample_size, noise_sd, 1)
tau_orig_design <- tau_new_func(MC_sample_size, noise_sd, nrow(orig_design_X))
# noise.var is the only way to have a stochastic emulator
noise.var <- list(tau1 = tau_orig_design$tau1^2,
                  tau2 = tau_orig_design$tau2^2)
```
Fit emulators
```{r}
model_f1 <- DiceKriging::km(formula=~1, design=orig_design_X, response=y1_orig, covtype="gauss", noise.var=noise.var$tau1)
model_f2 <- DiceKriging::km(formula=~1, design=orig_design_X, response=y2_orig, covtype="gauss", noise.var=noise.var$tau2)
```
y1_new and y2_new will record new observations, prompted by EQI
```{r}
y1_new <- y1_orig
y2_new <- y2_orig
y_plot <- cbind(y1=as.vector(y1_new),y2=as.vector(y2_new), x.1=orig_design_X[,1], x.2=orig_design_X[,2])
```

Now we move to the EQI loop to sequentially add desing points to alter the Pareto front.

First, we select new points to calculate EQI at. Covers all the points in the ranges.
```{r}
newdata <- expand.grid(x.1 = seq(from=x_c_1_range[1], to=x_c_1_range[2], length.out = 100),
                       x.2 = seq(from=x_c_2_range[1], to=x_c_2_range[2], length.out = 100))
n_sample <- length(newdata[,1])
```

The next line checks which of the current design points exists in the newdata. This is nessessary for tau_new function
```{r}
des_rep <- design_repetitions(newdata, design_X)
```
The next line calculates the default tau_new if there were no repetitions.
```{r}
tau_new <- tau_new_func(MC_sample_size, noise_sd, n_sample)
```
Update the design locations that were repeated
```{r}
if(sum(des_rep)!=0){
  tau_new[des_rep[,2],] <- cbind(tau1=sqrt(tau_eq_sqrd(noise.var$tau1[des_rep[,1]],noise.var$tau1[des_rep[,1]])),
                                 tau2=sqrt(tau_eq_sqrd(noise.var$tau2[des_rep[,1]],noise.var$tau2[des_rep[,1]])))
}
```
Add constraint info for objectives (currently set to no constraints).
```{r}
ConstraintInfo <- NULL
# ConstraintInfo$ConstraintLimits<-matrix(c(2, 2),1,2)
# #Current observations to be compared against ConstraintLimits
# ConstraintInfo$y <- cbind(y1_new, y2_new)
```
Start the EQI loop
```{r}
reps <- NULL
for (i in 1:Nsteps) {
  #calculate EQI metric. Note that other outputs are Pareto front, design and quantile sd
  EQI_newdata <- mult_EQI(newdata,design_X, model_f1, model_f2, beta, tau_new)
  
  #stopping criterion
  # If all expected improvements are 0 -- stop (i.e. -log(0)=Inf)
  if (sum(EQI_newdata$metric==Inf, na.rm = T)==n_sample) break
  # If all expected improvements are the same - select point at random
  else if (length(unique(EQI_newdata$metric)) == 1) {
    # Select a point to add to design
    select_point <- sample(1:n_sample,1)
    # Add selected point to design
    design_X <- rbind(newdata[select_point,],design_X)
  }
  # If not all EQI are zero and not all the same -- standard case
  else{#find the design point with the highest EQI metric (min(-log(EQI)))
    best_X <- which.min(EQI_newdata$metric)
    #find the values of the best design points
    impr_x <- newdata[best_X,]
    repetition <- row.match(impr_x, design_X)
    # Update the design_X
    design_X <- rbind(impr_x,design_X)
  }
  
  
  # Run the model at the new design point (MC over x)
  data_points <-cbind(x <- matrix(rep(design_X[1,],each=MC_sample_size),MC_sample_size,2),
                      theta.1=runif(MC_sample_size,-pi,pi),
                      theta.2=rnorm(MC_sample_size,mu_2,sigma_2))
  
  
  
  if (is.na(repetition)) {
    # Update observations
    y1_new <- c(mean(f1(data_points[,1:2],data_points[,3:4])), y1_new)
    y2_new <- c(mean(f2(data_points[,1:2],data_points[,3:4])), y2_new)
    y1_all <- rbind(unlist(c(f1(data_points[,1:2],data_points[,3:4]), design_X[1,])),y1_all)
    y2_all <- rbind(unlist(c(f2(data_points[,1:2],data_points[,3:4]), design_X[1,])),y2_all)
    # Update the tunable future noise
    tau_at_best_X <-  tau_new_func(MC_sample_size, c(sd(y1_all[1,1:MC_sample_size]),sd(y2_all[1,1:MC_sample_size])), 1)
    tau_new[best_X,] <- cbind(tau1=sqrt(tau_eq_sqrd(tau_at_best_X$tau1^2,tau_at_best_X$tau1^2)),
                              tau2=sqrt(tau_eq_sqrd(tau_at_best_X$tau2^2,tau_at_best_X$tau2^2)))
    # Update the observations noise
    noise.var <- data.frame(tau1 = c(tau_at_best_X$tau1^2,noise.var$tau1),
                            tau2 = c(tau_at_best_X$tau2^2,noise.var$tau2))
    y_plot <- rbind(c(y1=mean(f1(data_points[,1:2],data_points[,3:4])), y2=mean(f2(data_points[,1:2],data_points[,3:4])),
                      x.1=design_X[1,1], x.2=design_X[1,2]), y_plot)
  }else{
    # Update observations
    y1_all <- rbind(unlist(c(f1(data_points[,1:2],data_points[,3:4]), design_X[1,])),y1_all)
    y2_all <- rbind(unlist(c(f2(data_points[,1:2],data_points[,3:4]), design_X[1,])),y2_all)
    y_plot <- rbind(c(y1=mean(f1(data_points[,1:2],data_points[,3:4])), y2=mean(f2(data_points[,1:2],data_points[,3:4])),
                      x.1=design_X[1,1], x.2=design_X[1,2]), y_plot)
    design_X <- design_X[-1,]
    y1_new[repetition] <- mean_obs(mean(f1(data_points[,1:2],data_points[,3:4])),y1_new[repetition],noise_sd[1]^2/MC_sample_size,noise.var$tau1[repetition])
    y2_new[repetition] <- mean_obs(mean(f2(data_points[,1:2],data_points[,3:4])),y2_new[repetition],noise_sd[2]^2/MC_sample_size,noise.var$tau2[repetition])
    
    # Update the tunable future noise
    tau_at_best_X <- tau_new_func(MC_sample_size, c(sd(y1_all[1,1:MC_sample_size]),sd(y2_all[1,1:MC_sample_size])), 1)
    tau_new[best_X,] <- cbind(tau1=sqrt(tau_eq_sqrd((tau_new[best_X,]$tau1)^2,tau_at_best_X$tau1^2)),
                              tau2=sqrt(tau_eq_sqrd((tau_new[best_X,]$tau2)^2,tau_at_best_X$tau2^2)))
    # Update the observations noise
    noise.var$tau1[repetition] <- tau_eq_sqrd(noise.var$tau1[repetition], tau_at_best_X$tau1^2)
    noise.var$tau2[repetition] <- tau_eq_sqrd(noise.var$tau2[repetition], tau_at_best_X$tau2^2)
    
    reps <- c(reps, repetition)
  }
  model_f1 <- km(formula=~1, design=design_X, response=y1_new, covtype="gauss", noise.var=noise.var$tau1)
  model_f2 <- km(formula=~1, design=design_X, response=y2_new, covtype="gauss", noise.var=noise.var$tau2)
}
```

The following code reproduce the plots from the paper (!!!!!!!!!!!!!!! ref !!!!!!!!!!!!!!!!!!) . Note that the legend is conditional on whether there were repeated observations.
```{r}
###################################################################################################
#################################              Plot 2-D           #################################
###################################################################################################
grey <- '#e6e6e6' #230, 230, 230
darkgrey <- '#bdbdbd' #189, 189, 189
green <- '#9ad7cf' #154, 215, 207
purple <- '#b8b1db' #184, 177, 219
library(scales)


uncert <- alpha(purple, 0.2)
pareto <- green
observ <- darkgrey

#Get the EQI's outputs
EQI_newdata <- mult_EQI(newdata,design_X, model_f1, model_f2, beta, tau_new, ConstraintInfo)
p1 <- predict.km(model_f1, newdata, "UK")
p2 <- predict.km(model_f2, newdata, "UK")

# Plot sampled points
epsilon <- .1
par(mfrow=c(1,1))
# Create empty plot with labels and limits
plot(1, type="n", main=bquote(a==.(a)),
     xlab=expression(f[1](bold(x)[c],bold(x)[e])),
     ylab=expression(f[2](bold(x)[c],bold(x)[e])),
     ylim=c(min(y2_new)-2*epsilon,max(y2_new)+2*epsilon),
     xlim=c(min(y1_new)-4*epsilon,max(y1_new)+4*epsilon))

# Add uncertainty (MC samples)
Pareto_front_X <- EQI_newdata$PX
density_all <- list(NA)
for(i in 1:dim(Pareto_front_X)[1]){
  cont_index <- which(y1_all[,MC_sample_size+1]==Pareto_front_X[i,1] & y1_all[,MC_sample_size+2]==Pareto_front_X[i,2])
  y <- as.vector(unlist(y2_all[cont_index,1:MC_sample_size]))
  x <- as.vector(unlist(y1_all[cont_index,1:MC_sample_size]))
  density <- kde2d(x,y,n=20)
  density_all[[i]] <- density
  points(x,y,col=uncert, lwd=.1, pch=19)
}

# Calculate function values without environmental variables
f1_no_noise <- function(x){
  x1 <- unlist(x[,1])
  x2 <- unlist(x[,2])
  1-sin(x1)+x2/10
}
f2_no_noise <- function(x){
  x1 <- unlist(x[,1])
  x2 <- unlist(x[,2])
  1-cos(x1)+x2/3
}
y_1_for_pareto <- f1_no_noise(newdata)
y_2_for_pareto <- f2_no_noise(newdata)
# Find the actual Pareto front
find_pareto <- pareto_front(y_1_for_pareto,y_2_for_pareto,newdata)
# Add the actual pareto fron to the plot
points(x=find_pareto$y1, y=find_pareto$y2,col=purple, lwd=6, pch=19, type = "l")

# Add observations (y's)
points(y1_new,y2_new,col=observ,cex=2, pch=19, xlab='f1', ylab='f2',
       ylim=c(min(y2_new)-.3*epsilon,max(y2_new)+epsilon),
       xlim=c(min(y1_new)-epsilon,max(y1_new)+3*epsilon))

# Add Pareto front (qunatiles from EQI algorithm)
points(EQI_newdata$Pq1,EQI_newdata$Pq2,col=pareto, lwd=6, pch=19, type = "s")

# Find, add and label the order of new observations, added by the EQI algorithm
y_plot <- cbind(y_plot, i=dim(y_plot)[1]:1)
points(y_plot[1:Nsteps,1], y_plot[1:Nsteps,2], cex=2, pch=19, col=2)
text(y_plot[1:Nsteps,2]~y_plot[1:Nsteps,1], labels=y_plot[1:Nsteps,5]-n_sample_points, cex=1.5, font=2, pos=4)

# Find, add and connect the repeatedobservations to the final mean obseration
legend_ind <- NULL
for (i in 1:dim(design_X)[1]){
  cont_index <- which(y_plot[,3]==EQI_newdata$PX[i,1] & y_plot[,4]==EQI_newdata$PX[i,2])
  if(length(cont_index)>1){
    j <- which(design_X[,1]==EQI_newdata$PX[i,1] & design_X[,2]==EQI_newdata$PX[i,2])
    # connect repeated observations to the mean (final) observation
    segments(x0 = y_plot[cont_index,1], y0 = y_plot[cont_index,2],
             x1 = rep(y1_new[j],length(cont_index)), y1 = rep(y2_new[j],length(cont_index)),
             col = 2,
             lwd = 2,
             lty = 3)
    points(y1_new[j],y2_new[j], pch = 1, cex=2, col=1, lwd=2)
    legend_ind <- cont_index
  }
}
# Add legend
if(length(legend_ind)>1){
  legend(x= "topright",
         legend=c('observations','EQI pareto front',
                  'uncertainty', "new design points",
                  "real pareto front", "repeated observations",
                  "pooled observations"),
         col=c(observ, pareto, uncert, 2, purple, 1, 2),
         pch = c(19, 19, 19, 19, NA, 1, NA),
         lty=c(NA, 1 , NA, NA, 1,NA, 3),
         lwd=c(6,6,.1, 6, 6,2,2)
  )
}else{
  legend(x= "topright", 
         legend=c('observations','EQI pareto front',
                  'uncertainty', "new design points"),
         col=c(observ, pareto, uncert, 2),
         pch = c(19, 19, 19, 19),
         lty=c(NA, 1 , NA, NA),
         lwd=c(6,6,.1, 6)
  )
}
```
New design points plot.
```{r}
###################################################################################################
#################################     Sequential design  points   #################################
###################################################################################################

plot(y_plot[(Nsteps+1):nrow(y_plot),3],y_plot[(Nsteps+1):nrow(y_plot),4], xlab=expression(x[1]), ylab=expression(x[2]),
     xlim=(x_c_1_range+c(-.1,.1)), ylim=(x_c_2_range+c(-.1,.1)), main=paste("Original and added new design points") , pch=19, cex=3, col=observ)
for(i in Nsteps:1){
  text(y_plot[i,4]~y_plot[i,3], labels=nrow(y_plot)-nrow(orig_design_X)-i+1, cex=1.5, font=2, pos=4, offset=1)
  points(y_plot[i,3],y_plot[i,4], xlab=expression(x[1]), ylab=expression(x[2]), pch=19, cex=3, col=2)
}
```
Simple plot.
```{r}
###################################################################################################
################################              Simple plot           ###############################
###################################################################################################

# Plot sampled points
epsilon <- .1
par(mfrow=c(1,1))
# Create empty plot with labels and limits
plot(1, type="n", main=bquote(a==.(a)),
     xlab=expression(f[1](bold(x)[c],bold(x)[e])),
     ylab=expression(f[2](bold(x)[c],bold(x)[e])),
     ylim=c(min(y2_new)-4*epsilon,max(y2_new)+4*epsilon),
     xlim=c(min(y1_new)-4*epsilon,max(y1_new)+4*epsilon))

# Add uncertainty (MC samples)
Pareto_front_X <- EQI_newdata$PX
density_all <- list(NA)
for(i in 1:dim(Pareto_front_X)[1]){
  cont_index <- which(y1_all[,MC_sample_size+1]==Pareto_front_X[i,1])
  y <- as.vector(unlist(y2_all[cont_index,1:MC_sample_size]))
  x <- as.vector(unlist(y1_all[cont_index,1:MC_sample_size]))
  density <- kde2d(x,y,n=20)
  density_all[[i]] <- density
  points(x,y,col=uncert, lwd=.1, pch=19)
}

# Add the actual pareto front to the plot
points(x=find_pareto$y1, y=find_pareto$y2,col=purple, lwd=6, pch=19, type = "l")

# Add observations (y's)
points(y1_new,y2_new,col=observ,cex=2, pch=19, xlab='f1', ylab='f2',
       ylim=c(min(y2_new)-.3*epsilon,max(y2_new)+epsilon),
       xlim=c(min(y1_new)-epsilon,max(y1_new)+3*epsilon))

# Add Pareto front (qunatiles from EQI algorithm)
points(EQI_newdata$Pq1,EQI_newdata$Pq2,col=pareto, lwd=6, pch=19, type = "s")

# Add legend
legend(x= "topright",
       legend=c('observations','EQI pareto front','uncertainty', 'real pareto front'),
       col=c(observ, pareto, uncert, purple),
       pch = c(19, 19, 19, NA),
       lty=c(NA, 1 , NA, 1),
       lwd=c(6,6,.1, 6)
)
```
