MO-E-EQI
================
2023-12-15

Supplemental materials for (!!!!!!!!!!!!!!! ref !!!!!!!!!!!!!!!!!!)

Authors: D. Semochkina, A.I.J. Forester, D.C. Woods

Corresponding Author: D. Semochkina (<d.semochkina@soton.ac.uk>), School
of Mathematical Sciences, University of Southampton, UK

This repository provides R package to caclulate multi-objective
Euclidian expected quantile improvement (MO-E-EQI) presented in the
manuscript.

First install our package MOEEQI from git. We need the standard package
‘devtools’ to add our package off git.

``` r
#install.packages("devtools")
```

This is the standard way to import an R package into the current
session.

``` r
library("devtools")
```

    ## Loading required package: usethis

Now we need to build our package MOEEQI from git.

``` r
install_github("DashaMurasha/MO-E-EQI")
```

    ## Skipping install of 'MOEEQI' from a github remote, the SHA1 (9dccdd55) has not changed since last install.
    ##   Use `force = TRUE` to force installation

Note that the above instructions should only need running once in order
to install our package. After which we can just run:

``` r
library("MOEEQI")
```

    ## Loading required package: MASS

    ## Loading required package: DiceKriging

    ## Loading required package: prodlim

Next, we move on to the example accompanying the paper (!!!!!!!!!!!!!!!
ref !!!!!!!!!!!!!!!!!!)

We first set the level of noise and define out funstions

``` r
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

``` r
# Number of repetitions of each observation (model run)
MC_sample_size <- 100
```

Choose the number of steps of the MO-E-EsQI sequential design loop

``` r
# Computational budget (steps in the loop)
Nsteps <- 9
```

Select initial design points.

``` r
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

Choose the metric option. Currently two options available, -log(EQI)
(‘NegLogEQI’) or -EQI (‘NegEQI’)

``` r
Option <- 'NegLogEQI'
```

select quantile level (see EQI (!!!!!!!!!!!!!!! ref !!!!!!!!!!!!!!!!!!)
for details)

``` r
beta <- .8
```

Define some objects to store results

``` r
y1_orig <- y2_orig <- epsilons_orig <- NULL
y1_all <- y2_all <- NULL
var1 <- var2 <- NULL
y_plot <- NULL
```

Define the noise level for the second environmental variable

``` r
mu_2 <- 0
sigma_2 <- 0.5
```

Run the model at the original design points

``` r
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

``` r
noise_sd <- sqrt(c(mean(var1),mean(var2)))
noise_orig_design_sd <- cbind(apply(y1_all[,1:MC_sample_size],1,sd),
                              apply(y2_all[,1:MC_sample_size],1,sd))
```

Write original design and the starting point of the final design.

``` r
design_X <- orig_design_X
```

Calculate future noise. Note that tau is standard deviation (not
variance).

``` r
tau_new <- tau_new_func(MC_sample_size, noise_sd, 1)
tau_orig_design <- tau_new_func(MC_sample_size, noise_sd, nrow(orig_design_X))
# noise.var is the only way to have a stochastic emulator
noise.var <- list(tau1 = tau_orig_design$tau1^2,
                  tau2 = tau_orig_design$tau2^2)
```

Fit emulators

``` r
model_f1 <- DiceKriging::km(formula=~1, design=orig_design_X, response=y1_orig, covtype="gauss", noise.var=noise.var$tau1)
```

    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ## [1] 0.0003524746 0.0003524746 0.0003524746 0.0003524746 0.0003524746
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01904093 3.124538 
    ##   - best initial criterion value(s) :  -0.4675368 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      0.46754  |proj g|=       0.5852
    ## At iterate     1  f =       0.3414  |proj g|=       0.63238
    ## At iterate     2  f =      0.29331  |proj g|=       0.63234
    ## At iterate     3  f =      0.13326  |proj g|=        1.1388
    ## At iterate     4  f =      0.11265  |proj g|=       0.48233
    ## At iterate     5  f =      0.10144  |proj g|=       0.22503
    ## At iterate     6  f =     0.097903  |proj g|=      0.042449
    ## At iterate     7  f =     0.097875  |proj g|=      0.013125
    ## At iterate     8  f =     0.097872  |proj g|=    0.00015417
    ## At iterate     9  f =     0.097872  |proj g|=    5.8032e-07
    ## 
    ## iterations 9
    ## function evaluations 12
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 5.80318e-07
    ## final function value 0.0978724
    ## 
    ## F = 0.0978724
    ## final  value 0.097872 
    ## converged

``` r
model_f2 <- DiceKriging::km(formula=~1, design=orig_design_X, response=y2_orig, covtype="gauss", noise.var=noise.var$tau2)
```

    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ## [1] 0.0005667531 0.0005667531 0.0005667531 0.0005667531 0.0005667531
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01772042 2.436042 
    ##   - best initial criterion value(s) :  -1.183274 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=       1.1833  |proj g|=       2.2722
    ## At iterate     1  f =       0.8479  |proj g|=       0.67667
    ## At iterate     2  f =      0.80505  |proj g|=       0.46317
    ## At iterate     3  f =       0.7774  |proj g|=       0.89553
    ## At iterate     4  f =      0.76783  |proj g|=       0.33474
    ## At iterate     5  f =      0.75119  |proj g|=       0.25486
    ## At iterate     6  f =      0.69685  |proj g|=       0.31439
    ## At iterate     7  f =      0.67788  |proj g|=       0.35814
    ## At iterate     8  f =      0.67549  |proj g|=       0.16854
    ## At iterate     9  f =      0.67473  |proj g|=      0.019548
    ## At iterate    10  f =      0.67468  |proj g|=     0.0086753
    ## At iterate    11  f =      0.67467  |proj g|=     0.0020236
    ## At iterate    12  f =      0.67467  |proj g|=    0.00012716
    ## At iterate    13  f =      0.67467  |proj g|=    7.3404e-06
    ## 
    ## iterations 13
    ## function evaluations 17
    ## segments explored during Cauchy searches 15
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 7.34041e-06
    ## final function value 0.674673
    ## 
    ## F = 0.674673
    ## final  value 0.674673 
    ## converged

y1_new and y2_new will record new observations, prompted by EQI

``` r
y1_new <- y1_orig
y2_new <- y2_orig
y_plot <- cbind(y1=as.vector(y1_new),y2=as.vector(y2_new), x.1=orig_design_X[,1], x.2=orig_design_X[,2])
```

Now we move to the EQI loop to sequentially add desing points to alter
the Pareto front.

First, we select new points to calculate EQI at. Covers all the points
in the ranges.

``` r
newdata <- expand.grid(x.1 = seq(from=x_c_1_range[1], to=x_c_1_range[2], length.out = 100),
                       x.2 = seq(from=x_c_2_range[1], to=x_c_2_range[2], length.out = 100))
n_sample <- length(newdata[,1])
```

The next line checks which of the current design points exists in the
newdata. This is nessessary for tau_new function

``` r
des_rep <- design_repetitions(newdata, design_X)
```

The next line calculates the default tau_new if there were no
repetitions.

``` r
tau_new <- tau_new_func(MC_sample_size, noise_sd, n_sample)
```

Update the design locations that were repeated

``` r
if(sum(des_rep)!=0){
  tau_new[des_rep[,2],] <- cbind(tau1=sqrt(tau_eq_sqrd(noise.var$tau1[des_rep[,1]],noise.var$tau1[des_rep[,1]])),
                                 tau2=sqrt(tau_eq_sqrd(noise.var$tau2[des_rep[,1]],noise.var$tau2[des_rep[,1]])))
}
```

Add constraint info for objectives (currently set to no constraints).

``` r
ConstraintInfo <- NULL
# ConstraintInfo$ConstraintLimits<-matrix(c(2, 2),1,2)
# #Current observations to be compared against ConstraintLimits
# ConstraintInfo$y <- cbind(y1_new, y2_new)
```

Start the EQI loop

``` r
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

    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ## [1] 0.0003091727 0.0003524746 0.0003524746 0.0003524746 0.0003524746
    ## [6] 0.0003524746
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.02111622 3.505512 
    ##   - best initial criterion value(s) :  0.7520537 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=     -0.75205  |proj g|=       2.1981
    ## At iterate     1  f =       -1.376  |proj g|=       0.30054
    ## At iterate     2  f =      -1.4711  |proj g|=       0.39159
    ## At iterate     3  f =      -1.6061  |proj g|=        2.3362
    ## At iterate     4  f =       -1.657  |proj g|=        0.4253
    ## At iterate     5  f =       -1.668  |proj g|=       0.18614
    ## At iterate     6  f =      -1.6727  |proj g|=      0.060488
    ## At iterate     7  f =      -1.6728  |proj g|=      0.017588
    ## At iterate     8  f =      -1.6728  |proj g|=     0.0016245
    ## At iterate     9  f =      -1.6728  |proj g|=    0.00015981
    ## At iterate    10  f =      -1.6728  |proj g|=    7.0461e-06
    ## 
    ## iterations 10
    ## function evaluations 13
    ## segments explored during Cauchy searches 12
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 7.04607e-06
    ## final function value -1.67281
    ## 
    ## F = -1.67281
    ## final  value -1.672811 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ## [1] 0.0005990853 0.0005667531 0.0005667531 0.0005667531 0.0005667531
    ## [6] 0.0005667531
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01747439 2.772552 
    ##   - best initial criterion value(s) :  -0.2618488 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      0.26185  |proj g|=       1.1207
    ## At iterate     1  f =    -0.001129  |proj g|=        2.5014
    ## At iterate     2  f =     -0.19567  |proj g|=        2.6489
    ## At iterate     3  f =     -0.30516  |proj g|=        1.9313
    ## At iterate     4  f =     -0.49446  |proj g|=       0.70938
    ## At iterate     5  f =     -0.53665  |proj g|=       0.63629
    ## At iterate     6  f =     -0.61439  |proj g|=        1.8386
    ## At iterate     7  f =     -0.67066  |proj g|=        1.8507
    ## At iterate     8  f =     -0.70294  |proj g|=       0.53838
    ## At iterate     9  f =     -0.72179  |proj g|=       0.08775
    ## At iterate    10  f =     -0.72249  |proj g|=      0.069335
    ## At iterate    11  f =     -0.72251  |proj g|=      0.056702
    ## At iterate    12  f =     -0.72257  |proj g|=      0.014499
    ## At iterate    13  f =     -0.72257  |proj g|=     0.0012583
    ## At iterate    14  f =     -0.72257  |proj g|=    1.4456e-05
    ## At iterate    15  f =     -0.72257  |proj g|=    4.5012e-06
    ## 
    ## iterations 15
    ## function evaluations 20
    ## segments explored during Cauchy searches 17
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 0
    ## norm of the final projected gradient 4.50116e-06
    ## final function value -0.722575
    ## 
    ## F = -0.722575
    ## final  value -0.722575 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ## [1] 0.0003177960 0.0003091727 0.0003524746 0.0003524746 0.0003524746
    ## [6] 0.0003524746 0.0003524746
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01923072 3.253937 
    ##   - best initial criterion value(s) :  1.938315 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -1.9383  |proj g|=       1.8735
    ## At iterate     1  f =      -2.4532  |proj g|=        3.0841
    ## At iterate     2  f =      -2.5468  |proj g|=         3.069
    ## At iterate     3  f =       -2.714  |proj g|=        1.7267
    ## At iterate     4  f =      -2.7879  |proj g|=       0.60785
    ## At iterate     5  f =      -2.8697  |proj g|=       0.28165
    ## At iterate     6  f =      -2.9638  |proj g|=       0.38584
    ## At iterate     7  f =      -2.9959  |proj g|=      0.023758
    ## At iterate     8  f =       -2.996  |proj g|=      0.027312
    ## At iterate     9  f =       -2.996  |proj g|=     0.0058631
    ## At iterate    10  f =       -2.996  |proj g|=    0.00021275
    ## At iterate    11  f =       -2.996  |proj g|=    1.4421e-06
    ## At iterate    12  f =       -2.996  |proj g|=    1.0568e-09
    ## 
    ## iterations 12
    ## function evaluations 17
    ## segments explored during Cauchy searches 15
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 1.05685e-09
    ## final function value -2.99599
    ## 
    ## F = -2.99599
    ## final  value -2.995989 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ## [1] 0.0005974442 0.0005990853 0.0005667531 0.0005667531 0.0005667531
    ## [6] 0.0005667531 0.0005667531
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01761617 2.997113 
    ##   - best initial criterion value(s) :  2.22081 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -2.2208  |proj g|=       2.4431
    ## At iterate     1  f =      -2.6494  |proj g|=        2.1772
    ## At iterate     2  f =      -2.7441  |proj g|=        1.6144
    ## At iterate     3  f =      -2.8549  |proj g|=       0.24323
    ## At iterate     4  f =      -2.8675  |proj g|=       0.19962
    ## At iterate     5  f =      -2.8744  |proj g|=       0.29708
    ## At iterate     6  f =      -2.8951  |proj g|=       0.82572
    ## At iterate     7  f =      -2.9207  |proj g|=       0.77779
    ## At iterate     8  f =      -2.9266  |proj g|=       0.22718
    ## At iterate     9  f =      -2.9277  |proj g|=      0.030196
    ## At iterate    10  f =      -2.9277  |proj g|=     0.0011405
    ## At iterate    11  f =      -2.9277  |proj g|=    0.00014488
    ## At iterate    12  f =      -2.9277  |proj g|=    2.6153e-07
    ## 
    ## iterations 12
    ## function evaluations 15
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 2.61531e-07
    ## final function value -2.92769
    ## 
    ## F = -2.92769
    ## final  value -2.927687 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ## [1] 0.0003328102 0.0003177960 0.0003091727 0.0003524746 0.0003524746
    ## [6] 0.0003524746 0.0003524746 0.0003524746
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.02015232 3.436791 
    ##   - best initial criterion value(s) :  4.689599 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -4.6896  |proj g|=       2.6685
    ## At iterate     1  f =      -4.7395  |proj g|=        2.3467
    ## At iterate     2  f =      -5.0067  |proj g|=         1.604
    ## At iterate     3  f =      -5.2574  |proj g|=        1.6206
    ## At iterate     4  f =      -5.3481  |proj g|=        1.2232
    ## At iterate     5  f =      -5.5693  |proj g|=       0.30345
    ## At iterate     6  f =      -5.5713  |proj g|=      0.092234
    ## At iterate     7  f =      -5.5714  |proj g|=      0.046566
    ## At iterate     8  f =      -5.5715  |proj g|=     0.0012681
    ## At iterate     9  f =      -5.5715  |proj g|=    2.0122e-05
    ## At iterate    10  f =      -5.5715  |proj g|=    4.3263e-08
    ## 
    ## iterations 10
    ## function evaluations 16
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 4.32634e-08
    ## final function value -5.57145
    ## 
    ## F = -5.57145
    ## final  value -5.571454 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ## [1] 0.0006194281 0.0005974442 0.0005990853 0.0005667531 0.0005667531
    ## [6] 0.0005667531 0.0005667531 0.0005667531
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01777209 2.995284 
    ##   - best initial criterion value(s) :  3.827709 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -3.8277  |proj g|=       2.2916
    ## At iterate     1  f =      -4.2436  |proj g|=        1.0241
    ## At iterate     2  f =      -4.6136  |proj g|=       0.97255
    ## At iterate     3  f =      -4.7605  |proj g|=        2.8406
    ## At iterate     4  f =      -4.9376  |proj g|=        1.1674
    ## At iterate     5  f =      -5.1622  |proj g|=       0.21401
    ## At iterate     6  f =      -5.2513  |proj g|=       0.45092
    ## At iterate     7  f =      -5.2534  |proj g|=       0.19633
    ## At iterate     8  f =      -5.2538  |proj g|=       0.09275
    ## At iterate     9  f =      -5.2539  |proj g|=      0.017428
    ## At iterate    10  f =      -5.2539  |proj g|=    0.00060952
    ## At iterate    11  f =      -5.2539  |proj g|=    9.0622e-06
    ## 
    ## iterations 11
    ## function evaluations 14
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 9.06216e-06
    ## final function value -5.25394
    ## 
    ## F = -5.25394
    ## final  value -5.253940 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ## [1] 0.0003913708 0.0003328102 0.0003177960 0.0003091727 0.0003524746
    ## [6] 0.0003524746 0.0003524746 0.0003524746 0.0003524746
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01992834 3.502256 
    ##   - best initial criterion value(s) :  7.149582 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -7.1496  |proj g|=       2.2469
    ## At iterate     1  f =      -7.6383  |proj g|=       0.26207
    ## At iterate     2  f =      -7.6799  |proj g|=        0.1564
    ## At iterate     3  f =      -7.7347  |proj g|=        1.3669
    ## At iterate     4  f =      -7.7616  |proj g|=        0.5833
    ## At iterate     5  f =      -7.7684  |proj g|=       0.13042
    ## At iterate     6  f =      -7.7694  |proj g|=      0.013711
    ## At iterate     7  f =      -7.7694  |proj g|=     0.0011702
    ## At iterate     8  f =      -7.7694  |proj g|=    0.00019814
    ## At iterate     9  f =      -7.7694  |proj g|=    2.8261e-05
    ## 
    ## iterations 9
    ## function evaluations 12
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 2.82608e-05
    ## final function value -7.76944
    ## 
    ## F = -7.76944
    ## final  value -7.769442 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ## [1] 0.0005910930 0.0006194281 0.0005974442 0.0005990853 0.0005667531
    ## [6] 0.0005667531 0.0005667531 0.0005667531 0.0005667531
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01954891 3.370988 
    ##   - best initial criterion value(s) :  7.622131 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -7.6221  |proj g|=       1.0312
    ## At iterate     1  f =      -7.6404  |proj g|=       0.24686
    ## At iterate     2  f =       -7.648  |proj g|=       0.23832
    ## At iterate     3  f =      -7.6789  |proj g|=       0.54836
    ## At iterate     4  f =      -7.6979  |proj g|=       0.52351
    ## At iterate     5  f =      -7.7055  |proj g|=       0.20467
    ## At iterate     6  f =       -7.706  |proj g|=      0.025608
    ## At iterate     7  f =       -7.706  |proj g|=      0.001031
    ## At iterate     8  f =       -7.706  |proj g|=    3.4434e-05
    ## At iterate     9  f =       -7.706  |proj g|=    9.4516e-07
    ## 
    ## iterations 9
    ## function evaluations 11
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 9.45161e-07
    ## final function value -7.70603
    ## 
    ## F = -7.70603
    ## final  value -7.706035 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ##  [1] 0.0003353447 0.0003913708 0.0003328102 0.0003177960 0.0003091727
    ##  [6] 0.0003524746 0.0003524746 0.0003524746 0.0003524746 0.0003524746
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01772193 2.983496 
    ##   - best initial criterion value(s) :  10.03151 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -10.032  |proj g|=       2.8235
    ## At iterate     1  f =      -10.272  |proj g|=         1.498
    ## At iterate     2  f =      -10.307  |proj g|=        1.0687
    ## At iterate     3  f =       -10.35  |proj g|=       0.99318
    ## At iterate     4  f =      -10.381  |proj g|=        1.1842
    ## At iterate     5  f =      -10.402  |proj g|=       0.27292
    ## At iterate     6  f =      -10.404  |proj g|=       0.20944
    ## At iterate     7  f =      -10.409  |proj g|=      0.026435
    ## At iterate     8  f =      -10.409  |proj g|=     0.0045574
    ## At iterate     9  f =      -10.409  |proj g|=    0.00030976
    ## At iterate    10  f =      -10.409  |proj g|=    3.6723e-05
    ## 
    ## iterations 10
    ## function evaluations 13
    ## segments explored during Cauchy searches 12
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 3.67231e-05
    ## final function value -10.4086
    ## 
    ## F = -10.4086
    ## final  value -10.408629 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ##  [1] 0.0004985622 0.0005910930 0.0006194281 0.0005974442 0.0005990853
    ##  [6] 0.0005667531 0.0005667531 0.0005667531 0.0005667531 0.0005667531
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01918881 3.295652 
    ##   - best initial criterion value(s) :  8.239286 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -8.2393  |proj g|=       2.2034
    ## At iterate     1  f =       -8.754  |proj g|=        1.1372
    ## At iterate     2  f =       -9.292  |proj g|=        1.0807
    ## At iterate     3  f =      -9.5038  |proj g|=        2.4495
    ## At iterate     4  f =      -9.6199  |proj g|=        1.1359
    ## At iterate     5  f =      -9.7383  |proj g|=       0.97224
    ## At iterate     6  f =      -9.8988  |proj g|=       0.52391
    ## At iterate     7  f =      -9.9011  |proj g|=       0.25346
    ## At iterate     8  f =      -9.9018  |proj g|=      0.021251
    ## At iterate     9  f =      -9.9018  |proj g|=     0.0035576
    ## At iterate    10  f =      -9.9018  |proj g|=    0.00094245
    ## At iterate    11  f =      -9.9018  |proj g|=    4.3749e-05
    ## 
    ## iterations 11
    ## function evaluations 15
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 4.37485e-05
    ## final function value -9.90185
    ## 
    ## F = -9.90185
    ## final  value -9.901849 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ##  [1] 0.0003353447 0.0001866001 0.0003328102 0.0003177960 0.0003091727
    ##  [6] 0.0003524746 0.0003524746 0.0003524746 0.0003524746 0.0003524746
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01781098 3.000973 
    ##   - best initial criterion value(s) :  9.857797 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -9.8578  |proj g|=        2.868
    ## At iterate     1  f =      -10.131  |proj g|=        1.9964
    ## At iterate     2  f =      -10.217  |proj g|=        1.3668
    ## At iterate     3  f =      -10.316  |proj g|=         2.152
    ## At iterate     4  f =      -10.369  |proj g|=        1.2904
    ## At iterate     5  f =      -10.759  |proj g|=       0.29439
    ## At iterate     6  f =      -10.826  |proj g|=       0.45661
    ## At iterate     7  f =      -10.829  |proj g|=       0.10074
    ## At iterate     8  f =      -10.829  |proj g|=      0.070287
    ## At iterate     9  f =      -10.829  |proj g|=      0.042174
    ## At iterate    10  f =      -10.829  |proj g|=     6.935e-05
    ## At iterate    11  f =      -10.829  |proj g|=    1.1389e-05
    ## 
    ## iterations 11
    ## function evaluations 14
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 1.13888e-05
    ## final function value -10.8294
    ## 
    ## F = -10.8294
    ## final  value -10.829431 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ##  [1] 0.0004985622 0.0002848471 0.0006194281 0.0005974442 0.0005990853
    ##  [6] 0.0005667531 0.0005667531 0.0005667531 0.0005667531 0.0005667531
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01925487 3.309089 
    ##   - best initial criterion value(s) :  9.343609 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -9.3436  |proj g|=       3.1551
    ## At iterate     1  f =      -9.7098  |proj g|=        2.1669
    ## At iterate     2  f =      -9.7986  |proj g|=        1.6305
    ## At iterate     3  f =      -9.8516  |proj g|=        1.1205
    ## At iterate     4  f =      -9.8615  |proj g|=       0.34558
    ## At iterate     5  f =      -9.8687  |proj g|=       0.20072
    ## At iterate     6  f =      -9.8807  |proj g|=       0.10939
    ## At iterate     7  f =      -9.8808  |proj g|=      0.015908
    ## At iterate     8  f =      -9.8808  |proj g|=     0.0047354
    ## At iterate     9  f =      -9.8808  |proj g|=     0.0033139
    ## At iterate    10  f =      -9.8808  |proj g|=     0.0005072
    ## At iterate    11  f =      -9.8808  |proj g|=    4.0442e-05
    ## 
    ## iterations 11
    ## function evaluations 14
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 4.04418e-05
    ## final function value -9.88082
    ## 
    ## F = -9.88082
    ## final  value -9.880821 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ##  [1] 0.0003303053 0.0003353447 0.0001866001 0.0003328102 0.0003177960
    ##  [6] 0.0003091727 0.0003524746 0.0003524746 0.0003524746 0.0003524746
    ## [11] 0.0003524746
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01652765 2.774934 
    ##   - best initial criterion value(s) :  12.33268 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -12.333  |proj g|=      0.92373
    ## At iterate     1  f =      -12.583  |proj g|=        1.1762
    ## At iterate     2  f =      -12.813  |proj g|=        1.1827
    ## At iterate     3  f =      -13.577  |proj g|=       0.37041
    ## At iterate     4  f =      -13.592  |proj g|=       0.50236
    ## At iterate     5  f =      -13.598  |proj g|=       0.10832
    ## At iterate     6  f =      -13.598  |proj g|=      0.019917
    ## At iterate     7  f =      -13.598  |proj g|=     0.0053919
    ## At iterate     8  f =      -13.598  |proj g|=    0.00012904
    ## At iterate     9  f =      -13.598  |proj g|=    5.8076e-06
    ## 
    ## iterations 9
    ## function evaluations 12
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 5.80758e-06
    ## final function value -13.5982
    ## 
    ## F = -13.5982
    ## final  value -13.598210 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ##  [1] 0.0005659895 0.0004985622 0.0002848471 0.0006194281 0.0005974442
    ##  [6] 0.0005990853 0.0005667531 0.0005667531 0.0005667531 0.0005667531
    ## [11] 0.0005667531
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01907129 3.30781 
    ##   - best initial criterion value(s) :  11.86881 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -11.869  |proj g|=       3.0926
    ## At iterate     1  f =      -11.968  |proj g|=       0.25445
    ## At iterate     2  f =      -11.971  |proj g|=       0.16144
    ## At iterate     3  f =      -11.976  |proj g|=       0.27753
    ## At iterate     4  f =      -11.983  |proj g|=       0.52301
    ## At iterate     5  f =      -11.998  |proj g|=       0.70059
    ## At iterate     6  f =      -12.006  |proj g|=       0.30813
    ## At iterate     7  f =      -12.008  |proj g|=       0.10813
    ## At iterate     8  f =      -12.009  |proj g|=      0.046066
    ## At iterate     9  f =      -12.009  |proj g|=      0.008446
    ## At iterate    10  f =      -12.009  |proj g|=    0.00061439
    ## At iterate    11  f =      -12.009  |proj g|=    2.4758e-05
    ## 
    ## iterations 11
    ## function evaluations 14
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 2.47584e-05
    ## final function value -12.0086
    ## 
    ## F = -12.0086
    ## final  value -12.008551 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ##  [1] 0.0003540184 0.0003303053 0.0003353447 0.0001866001 0.0003328102
    ##  [6] 0.0003177960 0.0003091727 0.0003524746 0.0003524746 0.0003524746
    ## [11] 0.0003524746 0.0003524746
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01622371 2.800693 
    ##   - best initial criterion value(s) :  14.03079 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -14.031  |proj g|=       2.6507
    ## At iterate     1  f =      -14.306  |proj g|=        2.0235
    ## At iterate     2  f =      -14.845  |proj g|=        1.9481
    ## At iterate     3  f =      -15.814  |proj g|=       0.38336
    ## At iterate     4  f =      -16.033  |proj g|=        1.3121
    ## At iterate     5  f =      -16.185  |proj g|=        1.4981
    ## At iterate     6  f =      -16.254  |proj g|=       0.37695
    ## At iterate     7  f =      -16.259  |proj g|=       0.25495
    ## At iterate     8  f =      -16.266  |proj g|=      0.067577
    ## At iterate     9  f =      -16.266  |proj g|=     0.0045063
    ## At iterate    10  f =      -16.266  |proj g|=    8.5833e-05
    ## At iterate    11  f =      -16.266  |proj g|=    1.1002e-05
    ## 
    ## iterations 11
    ## function evaluations 15
    ## segments explored during Cauchy searches 16
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 1.10016e-05
    ## final function value -16.2664
    ## 
    ## F = -16.2664
    ## final  value -16.266417 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ##  [1] 0.0006134840 0.0005659895 0.0004985622 0.0002848471 0.0006194281
    ##  [6] 0.0005974442 0.0005990853 0.0005667531 0.0005667531 0.0005667531
    ## [11] 0.0005667531 0.0005667531
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01901343 3.360375 
    ##   - best initial criterion value(s) :  13.08921 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -13.089  |proj g|=       2.3138
    ## At iterate     1  f =      -14.123  |proj g|=       0.28716
    ## At iterate     2  f =       -14.17  |proj g|=       0.56604
    ## At iterate     3  f =      -14.185  |proj g|=       0.38602
    ## At iterate     4  f =       -14.21  |proj g|=       0.33486
    ## At iterate     5  f =      -14.229  |proj g|=       0.24771
    ## At iterate     6  f =       -14.23  |proj g|=      0.019731
    ## At iterate     7  f =       -14.23  |proj g|=      0.001028
    ## At iterate     8  f =       -14.23  |proj g|=    4.5178e-06
    ## 
    ## iterations 8
    ## function evaluations 11
    ## segments explored during Cauchy searches 10
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 4.51783e-06
    ## final function value -14.2298
    ## 
    ## F = -14.2298
    ## final  value -14.229813 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ##  [1] 0.0003109447 0.0003540184 0.0003303053 0.0003353447 0.0001866001
    ##  [6] 0.0003328102 0.0003177960 0.0003091727 0.0003524746 0.0003524746
    ## [11] 0.0003524746 0.0003524746 0.0003524746
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01651629 2.844259 
    ##   - best initial criterion value(s) :  17.91414 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -17.914  |proj g|=        2.163
    ## At iterate     1  f =      -18.321  |proj g|=       0.41215
    ## At iterate     2  f =      -18.475  |proj g|=       0.52502
    ## At iterate     3  f =      -18.754  |proj g|=        2.6601
    ## At iterate     4  f =      -18.832  |proj g|=       0.80839
    ## At iterate     5  f =      -18.855  |proj g|=       0.27308
    ## At iterate     6  f =      -18.859  |proj g|=      0.075696
    ## At iterate     7  f =      -18.859  |proj g|=      0.031202
    ## At iterate     8  f =      -18.859  |proj g|=    0.00099597
    ## At iterate     9  f =      -18.859  |proj g|=    3.8943e-05
    ## 
    ## iterations 9
    ## function evaluations 12
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 3.89432e-05
    ## final function value -18.8586
    ## 
    ## F = -18.8586
    ## final  value -18.858576 
    ## converged
    ## 
    ## optimisation start
    ## ------------------
    ## * estimation method   : MLE 
    ## * optimisation method : BFGS 
    ## * analytical gradient : used
    ## * trend model : ~1
    ## * covariance model : 
    ##   - type :  gauss 
    ##   - noise variances :
    ##  [1] 0.0006399401 0.0006134840 0.0005659895 0.0004985622 0.0002848471
    ##  [6] 0.0006194281 0.0005974442 0.0005990853 0.0005667531 0.0005667531
    ## [11] 0.0005667531 0.0005667531 0.0005667531
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01970439 3.505187 
    ##   - best initial criterion value(s) :  16.55798 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -16.558  |proj g|=       2.1561
    ## At iterate     1  f =      -16.678  |proj g|=        1.2226
    ## At iterate     2  f =      -16.709  |proj g|=       0.26596
    ## At iterate     3  f =      -16.715  |proj g|=      0.071478
    ## At iterate     4  f =      -16.715  |proj g|=      0.012379
    ## At iterate     5  f =      -16.715  |proj g|=     0.0025007
    ## At iterate     6  f =      -16.715  |proj g|=    0.00028581
    ## At iterate     7  f =      -16.715  |proj g|=    6.4704e-05
    ## 
    ## iterations 7
    ## function evaluations 10
    ## segments explored during Cauchy searches 10
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 6.47035e-05
    ## final function value -16.715
    ## 
    ## F = -16.715
    ## final  value -16.715042 
    ## converged

The following code reproduce the plots from the paper (!!!!!!!!!!!!!!!
ref !!!!!!!!!!!!!!!!!!) . Note that the legend is conditional on whether
there were repeated observations.

``` r
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

![](MO-E-EQI_files/figure-gfm/unnamed-chunk-25-1.png)<!-- --> New design
points plot.

``` r
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

![](MO-E-EQI_files/figure-gfm/unnamed-chunk-26-1.png)<!-- --> Simple
plot.

``` r
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

![](MO-E-EQI_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->
