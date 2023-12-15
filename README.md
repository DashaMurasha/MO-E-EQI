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
install.packages("devtools")
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

    ## Downloading GitHub repo DashaMurasha/MO-E-EQI@HEAD

    ## 
    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##      checking for file ‘/private/var/folders/rv/328k5chd2px4l66r9w03g4lm0000gp/T/RtmpZV8KmA/remotesa41660c75eb4/DashaMurasha-MO-E-EQI-c742450/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/rv/328k5chd2px4l66r9w03g4lm0000gp/T/RtmpZV8KmA/remotesa41660c75eb4/DashaMurasha-MO-E-EQI-c742450/DESCRIPTION’
    ##   ─  preparing ‘MOEEQI’:
    ##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
    ##      Warning: design_repetitions.Rd:59: unexpected END_OF_INPUT '
    ##    '
    ##      Warning: mean_obs.Rd:71: unexpected END_OF_INPUT '
    ##    '
    ##      Warning: pareto_front.Rd:84: unexpected END_OF_INPUT '
    ##    '
    ##      Warning: tau_eq_sqrd.Rd:63: unexpected END_OF_INPUT '
    ##    '
    ##   ─  checking for LF line-endings in source and make files and shell scripts
    ##   ─  checking for empty or unneeded directories
    ##      Omitted ‘LazyData’ from DESCRIPTION
    ##   ─  building ‘MOEEQI_0.1.0.tar.gz’
    ##      
    ## 

Note that the above instructions should only need running once in order
to install our package. After which we can just run:

``` r
library("MOEEQI")
```

    ## Loading required package: MASS

    ## Loading required package: DiceKriging

    ## Loading required package: prodlim

Next, we move on to the example accompanying the paper (!!!!!!!!!!!!!!!
ref !!!!!!!!!!!!!!!!!!). This example is available in the
‘test_function_EQI.R’ file.

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
    ## [1] 0.0003392235 0.0003392235 0.0003392235 0.0003392235 0.0003392235
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01775105 2.566296 
    ##   - best initial criterion value(s) :  -1.202881 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=       1.2029  |proj g|=       2.4857
    ## At iterate     1  f =      0.53036  |proj g|=       0.91993
    ## At iterate     2  f =      0.48378  |proj g|=       0.85468
    ## At iterate     3  f =      0.34604  |proj g|=         1.638
    ## At iterate     4  f =     -0.15933  |proj g|=       0.66323
    ## At iterate     5  f =     -0.17405  |proj g|=       0.03434
    ## At iterate     6  f =     -0.17411  |proj g|=      0.023027
    ## At iterate     7  f =     -0.17412  |proj g|=      0.018012
    ## At iterate     8  f =     -0.17412  |proj g|=    0.00012556
    ## At iterate     9  f =     -0.17412  |proj g|=    1.4813e-06
    ## 
    ## iterations 9
    ## function evaluations 12
    ## segments explored during Cauchy searches 10
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 1.48131e-06
    ## final function value -0.17412
    ## 
    ## F = -0.17412
    ## final  value -0.174120 
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
    ## [1] 0.0005671853 0.0005671853 0.0005671853 0.0005671853 0.0005671853
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01618822 2.580716 
    ##   - best initial criterion value(s) :  -1.098281 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=       1.0983  |proj g|=       1.1168
    ## At iterate     1  f =       1.0388  |proj g|=        1.0496
    ## At iterate     2  f =       1.0173  |proj g|=       0.54991
    ## At iterate     3  f =       1.0007  |proj g|=       0.41963
    ## At iterate     4  f =      0.95689  |proj g|=       0.39897
    ## At iterate     5  f =         0.91  |proj g|=       0.23299
    ## At iterate     6  f =      0.89033  |proj g|=        0.2328
    ## At iterate     7  f =       0.8865  |proj g|=       0.10415
    ## At iterate     8  f =      0.88527  |proj g|=      0.043169
    ## At iterate     9  f =      0.88507  |proj g|=      0.010983
    ## At iterate    10  f =      0.88507  |proj g|=     0.0012292
    ## At iterate    11  f =      0.88507  |proj g|=    6.3104e-05
    ## At iterate    12  f =      0.88507  |proj g|=    1.0137e-05
    ## 
    ## iterations 12
    ## function evaluations 16
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 0
    ## norm of the final projected gradient 1.01375e-05
    ## final function value 0.885065
    ## 
    ## F = 0.885065
    ## final  value 0.885065 
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
    ## [1] 0.0003469849 0.0003392235 0.0003392235 0.0003392235 0.0003392235
    ## [6] 0.0003392235
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.0148499 2.239366 
    ##   - best initial criterion value(s) :  0.9712027 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -0.9712  |proj g|=       2.0463
    ## At iterate     1  f =      -1.4872  |proj g|=       0.70905
    ## At iterate     2  f =      -1.8629  |proj g|=        1.0496
    ## At iterate     3  f =      -1.9792  |proj g|=        1.7797
    ## At iterate     4  f =      -2.1209  |proj g|=       0.64397
    ## At iterate     5  f =      -2.2152  |proj g|=        1.0647
    ## At iterate     6  f =      -2.2905  |proj g|=       0.49924
    ## At iterate     7  f =      -2.3118  |proj g|=      0.013526
    ## At iterate     8  f =      -2.3118  |proj g|=     0.0067672
    ## At iterate     9  f =      -2.3118  |proj g|=     0.0022992
    ## At iterate    10  f =      -2.3118  |proj g|=    7.6809e-06
    ## At iterate    11  f =      -2.3118  |proj g|=     6.547e-08
    ## 
    ## iterations 11
    ## function evaluations 16
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 6.54702e-08
    ## final function value -2.31179
    ## 
    ## F = -2.31179
    ## final  value -2.311790 
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
    ## [1] 0.0005610908 0.0005671853 0.0005671853 0.0005671853 0.0005671853
    ## [6] 0.0005671853
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01711059 2.610959 
    ##   - best initial criterion value(s) :  0.4984574 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=     -0.49846  |proj g|=      0.95795
    ## At iterate     1  f =     -0.77216  |proj g|=        2.4264
    ## At iterate     2  f =     -0.91538  |proj g|=       0.55739
    ## At iterate     3  f =     -0.96549  |proj g|=       0.37201
    ## At iterate     4  f =      -1.0967  |proj g|=       0.50253
    ## At iterate     5  f =      -1.1548  |proj g|=       0.36169
    ## At iterate     6  f =      -1.1576  |proj g|=       0.16488
    ## At iterate     7  f =      -1.1584  |proj g|=      0.016218
    ## At iterate     8  f =      -1.1585  |proj g|=      0.005251
    ## At iterate     9  f =      -1.1585  |proj g|=    0.00052973
    ## At iterate    10  f =      -1.1585  |proj g|=    3.0478e-05
    ## At iterate    11  f =      -1.1585  |proj g|=    2.5147e-06
    ## 
    ## iterations 11
    ## function evaluations 13
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 2.51475e-06
    ## final function value -1.15846
    ## 
    ## F = -1.15846
    ## final  value -1.158458 
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
    ## [1] 0.0003651903 0.0003469849 0.0003392235 0.0003392235 0.0003392235
    ## [6] 0.0003392235 0.0003392235
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01418728 2.20417 
    ##   - best initial criterion value(s) :  3.566497 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -3.5665  |proj g|=       2.0974
    ## At iterate     1  f =      -3.8088  |proj g|=        1.9136
    ## At iterate     2  f =       -3.892  |proj g|=        1.3335
    ## At iterate     3  f =       -4.075  |proj g|=        2.0363
    ## At iterate     4  f =      -4.1676  |proj g|=        2.0162
    ## At iterate     5  f =       -4.458  |proj g|=        1.9035
    ## At iterate     6  f =      -4.5151  |proj g|=       0.56121
    ## At iterate     7  f =      -4.5223  |proj g|=      0.034826
    ## At iterate     8  f =      -4.5224  |proj g|=      0.039177
    ## At iterate     9  f =      -4.5224  |proj g|=      0.013306
    ## At iterate    10  f =      -4.5224  |proj g|=    0.00013727
    ## At iterate    11  f =      -4.5224  |proj g|=    1.2598e-05
    ## 
    ## iterations 11
    ## function evaluations 16
    ## segments explored during Cauchy searches 14
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 1.25976e-05
    ## final function value -4.52239
    ## 
    ## F = -4.52239
    ## final  value -4.522394 
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
    ## [1] 0.0004803086 0.0005610908 0.0005671853 0.0005671853 0.0005671853
    ## [6] 0.0005671853 0.0005671853
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01426238 2.276984 
    ##   - best initial criterion value(s) :  3.212926 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -3.2129  |proj g|=       2.0122
    ## At iterate     1  f =      -3.3162  |proj g|=        1.1101
    ## At iterate     2  f =      -3.3389  |proj g|=       0.81663
    ## At iterate     3  f =      -3.3774  |proj g|=       0.54486
    ## At iterate     4  f =      -3.3915  |proj g|=       0.44729
    ## At iterate     5  f =      -3.4119  |proj g|=       0.14947
    ## At iterate     6  f =      -3.4124  |proj g|=      0.077798
    ## At iterate     7  f =      -3.4127  |proj g|=     0.0072894
    ## At iterate     8  f =      -3.4127  |proj g|=     0.0016862
    ## At iterate     9  f =      -3.4127  |proj g|=    0.00013469
    ## At iterate    10  f =      -3.4127  |proj g|=    9.3953e-06
    ## 
    ## iterations 10
    ## function evaluations 13
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 9.39526e-06
    ## final function value -3.41271
    ## 
    ## F = -3.41271
    ## final  value -3.412705 
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
    ## [1] 0.0003158729 0.0003651903 0.0003469849 0.0003392235 0.0003392235
    ## [6] 0.0003392235 0.0003392235 0.0003392235
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01564072 2.392768 
    ##   - best initial criterion value(s) :  5.904334 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -5.9043  |proj g|=       2.2672
    ## At iterate     1  f =      -6.0611  |proj g|=        1.9511
    ## At iterate     2  f =      -6.1457  |proj g|=        1.3419
    ## At iterate     3  f =      -6.2307  |proj g|=        2.2295
    ## At iterate     4  f =      -6.4845  |proj g|=       0.32164
    ## At iterate     5  f =      -6.8603  |proj g|=       0.22284
    ## At iterate     6  f =      -6.8609  |proj g|=       0.23208
    ## At iterate     7  f =      -6.8627  |proj g|=      0.063305
    ## At iterate     8  f =      -6.8628  |proj g|=      0.035459
    ## At iterate     9  f =      -6.8628  |proj g|=    0.00043302
    ## At iterate    10  f =      -6.8628  |proj g|=    5.4227e-06
    ## 
    ## iterations 10
    ## function evaluations 14
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 5.42266e-06
    ## final function value -6.86276
    ## 
    ## F = -6.86276
    ## final  value -6.862763 
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
    ## [1] 0.0005266163 0.0004803086 0.0005610908 0.0005671853 0.0005671853
    ## [6] 0.0005671853 0.0005671853 0.0005671853
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01503377 2.467766 
    ##   - best initial criterion value(s) :  4.480021 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=        -4.48  |proj g|=      0.86349
    ## At iterate     1  f =      -4.6607  |proj g|=        2.2384
    ## At iterate     2  f =      -4.8754  |proj g|=        2.2165
    ## At iterate     3  f =      -5.4525  |proj g|=        1.6334
    ## At iterate     4  f =      -5.4878  |proj g|=        0.7767
    ## At iterate     5  f =       -5.499  |proj g|=       0.28934
    ## At iterate     6  f =      -5.5011  |proj g|=      0.052875
    ## At iterate     7  f =      -5.5012  |proj g|=     0.0057432
    ## At iterate     8  f =      -5.5012  |proj g|=    0.00079671
    ## At iterate     9  f =      -5.5012  |proj g|=    5.9704e-05
    ## At iterate    10  f =      -5.5012  |proj g|=    2.4228e-06
    ## 
    ## iterations 10
    ## function evaluations 13
    ## segments explored during Cauchy searches 12
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 2.42281e-06
    ## final function value -5.5012
    ## 
    ## F = -5.5012
    ## final  value -5.501203 
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
    ## [1] 0.0003436513 0.0003158729 0.0003651903 0.0003469849 0.0003392235
    ## [6] 0.0003392235 0.0003392235 0.0003392235 0.0003392235
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01606935 2.614687 
    ##   - best initial criterion value(s) :  8.989272 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -8.9893  |proj g|=       0.7221
    ## At iterate     1  f =      -9.0649  |proj g|=       0.62978
    ## At iterate     2  f =      -9.1322  |proj g|=       0.69304
    ## At iterate     3  f =      -9.4292  |proj g|=        1.1339
    ## At iterate     4  f =      -9.4382  |proj g|=       0.12761
    ## At iterate     5  f =      -9.4385  |proj g|=      0.027726
    ## At iterate     6  f =      -9.4385  |proj g|=     0.0034871
    ## At iterate     7  f =      -9.4385  |proj g|=     0.0013334
    ## At iterate     8  f =      -9.4385  |proj g|=    9.4706e-06
    ## 
    ## iterations 8
    ## function evaluations 11
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 9.47063e-06
    ## final function value -9.43847
    ## 
    ## F = -9.43847
    ## final  value -9.438474 
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
    ## [1] 0.0006088488 0.0005266163 0.0004803086 0.0005610908 0.0005671853
    ## [6] 0.0005671853 0.0005671853 0.0005671853 0.0005671853
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.0149562 2.423476 
    ##   - best initial criterion value(s) :  6.271993 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=       -6.272  |proj g|=       2.2164
    ## At iterate     1  f =        -7.04  |proj g|=        1.8707
    ## At iterate     2  f =      -7.1812  |proj g|=         0.821
    ## At iterate     3  f =      -7.4425  |proj g|=       0.91106
    ## At iterate     4  f =      -7.7907  |proj g|=       0.44265
    ## At iterate     5  f =      -7.7989  |proj g|=      0.044879
    ## At iterate     6  f =       -7.799  |proj g|=       0.04054
    ## At iterate     7  f =      -7.7991  |proj g|=      0.024628
    ## At iterate     8  f =      -7.7991  |proj g|=    0.00088742
    ## At iterate     9  f =      -7.7991  |proj g|=     1.799e-05
    ## At iterate    10  f =      -7.7991  |proj g|=    1.0329e-06
    ## 
    ## iterations 10
    ## function evaluations 13
    ## segments explored during Cauchy searches 12
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 1.0329e-06
    ## final function value -7.79913
    ## 
    ## F = -7.79913
    ## final  value -7.799126 
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
    ##  [1] 0.0003344754 0.0003436513 0.0003158729 0.0003651903 0.0003469849
    ##  [6] 0.0003392235 0.0003392235 0.0003392235 0.0003392235 0.0003392235
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01637301 2.824797 
    ##   - best initial criterion value(s) :  11.0679 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -11.068  |proj g|=       1.5876
    ## At iterate     1  f =      -11.088  |proj g|=       0.17618
    ## At iterate     2  f =       -11.15  |proj g|=       0.12582
    ## At iterate     3  f =       -11.15  |proj g|=      0.026538
    ## At iterate     4  f =       -11.15  |proj g|=      0.012101
    ## At iterate     5  f =       -11.15  |proj g|=      0.012381
    ## At iterate     6  f =       -11.15  |proj g|=     0.0016214
    ## At iterate     7  f =       -11.15  |proj g|=    0.00014464
    ## 
    ## iterations 7
    ## function evaluations 11
    ## segments explored during Cauchy searches 8
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 0.000144637
    ## final function value -11.1502
    ## 
    ## F = -11.1502
    ## final  value -11.150214 
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
    ##  [1] 0.0006635865 0.0006088488 0.0005266163 0.0004803086 0.0005610908
    ##  [6] 0.0005671853 0.0005671853 0.0005671853 0.0005671853 0.0005671853
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01582483 2.58658 
    ##   - best initial criterion value(s) :  8.566319 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -8.5663  |proj g|=       2.4389
    ## At iterate     1  f =      -9.3568  |proj g|=        1.0441
    ## At iterate     2  f =      -9.6935  |proj g|=        0.8175
    ## At iterate     3  f =      -10.241  |proj g|=       0.67216
    ## At iterate     4  f =      -10.305  |proj g|=       0.45269
    ## At iterate     5  f =      -10.341  |proj g|=       0.33838
    ## At iterate     6  f =      -10.346  |proj g|=        0.2819
    ## At iterate     7  f =       -10.35  |proj g|=      0.060882
    ## At iterate     8  f =       -10.35  |proj g|=      0.016462
    ## At iterate     9  f =       -10.35  |proj g|=    0.00065819
    ## At iterate    10  f =       -10.35  |proj g|=    7.3178e-05
    ## 
    ## iterations 10
    ## function evaluations 13
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 7.31777e-05
    ## final function value -10.3497
    ## 
    ## F = -10.3497
    ## final  value -10.349671 
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
    ##  [1] 0.0003749209 0.0003344754 0.0003436513 0.0003158729 0.0003651903
    ##  [6] 0.0003469849 0.0003392235 0.0003392235 0.0003392235 0.0003392235
    ## [11] 0.0003392235
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01670866 2.805063 
    ##   - best initial criterion value(s) :  12.89376 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -12.894  |proj g|=      0.77862
    ## At iterate     1  f =      -12.942  |proj g|=       0.63152
    ## At iterate     2  f =      -13.094  |proj g|=        1.3284
    ## At iterate     3  f =      -13.214  |proj g|=         1.259
    ## At iterate     4  f =       -13.42  |proj g|=       0.21807
    ## At iterate     5  f =      -13.425  |proj g|=      0.094969
    ## At iterate     6  f =      -13.425  |proj g|=      0.015896
    ## At iterate     7  f =      -13.425  |proj g|=    0.00071227
    ## At iterate     8  f =      -13.425  |proj g|=    2.6452e-05
    ## 
    ## iterations 8
    ## function evaluations 10
    ## segments explored during Cauchy searches 9
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 2.64521e-05
    ## final function value -13.4251
    ## 
    ## F = -13.4251
    ## final  value -13.425085 
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
    ##  [1] 0.0005011100 0.0006635865 0.0006088488 0.0005266163 0.0004803086
    ##  [6] 0.0005610908 0.0005671853 0.0005671853 0.0005671853 0.0005671853
    ## [11] 0.0005671853
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01650937 2.777449 
    ##   - best initial criterion value(s) :  11.54925 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -11.549  |proj g|=       2.4802
    ## At iterate     1  f =      -12.264  |proj g|=        1.0554
    ## At iterate     2  f =      -12.455  |proj g|=       0.53784
    ## At iterate     3  f =      -12.707  |proj g|=       0.44618
    ## At iterate     4  f =      -12.846  |proj g|=        2.4589
    ## At iterate     5  f =      -12.983  |proj g|=       0.44612
    ## At iterate     6  f =      -13.006  |proj g|=       0.23586
    ## At iterate     7  f =      -13.011  |proj g|=       0.22048
    ## At iterate     8  f =      -13.014  |proj g|=      0.085832
    ## At iterate     9  f =      -13.015  |proj g|=      0.017194
    ## At iterate    10  f =      -13.015  |proj g|=     0.0042853
    ## At iterate    11  f =      -13.015  |proj g|=    0.00040626
    ## At iterate    12  f =      -13.015  |proj g|=     1.517e-06
    ## 
    ## iterations 12
    ## function evaluations 15
    ## segments explored during Cauchy searches 15
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 1.51699e-06
    ## final function value -13.0147
    ## 
    ## F = -13.0147
    ## final  value -13.014684 
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
    ##  [1] 0.0003311965 0.0003749209 0.0003344754 0.0003436513 0.0003158729
    ##  [6] 0.0003651903 0.0003469849 0.0003392235 0.0003392235 0.0003392235
    ## [11] 0.0003392235 0.0003392235
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01683845 2.94033 
    ##   - best initial criterion value(s) :  14.56848 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -14.568  |proj g|=       2.1051
    ## At iterate     1  f =      -15.008  |proj g|=       0.60617
    ## At iterate     2  f =      -15.218  |proj g|=        0.9523
    ## At iterate     3  f =      -15.661  |proj g|=         2.797
    ## At iterate     4  f =      -16.029  |proj g|=        1.8045
    ## At iterate     5  f =      -16.061  |proj g|=       0.87013
    ## At iterate     6  f =      -16.075  |proj g|=       0.32004
    ## At iterate     7  f =      -16.076  |proj g|=       0.15091
    ## At iterate     8  f =      -16.076  |proj g|=       0.05735
    ## At iterate     9  f =      -16.076  |proj g|=      0.007631
    ## At iterate    10  f =      -16.076  |proj g|=    3.3577e-05
    ## At iterate    11  f =      -16.076  |proj g|=    3.3454e-05
    ## 
    ## iterations 11
    ## function evaluations 14
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 3.34542e-05
    ## final function value -16.0761
    ## 
    ## F = -16.0761
    ## final  value -16.076108 
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
    ##  [1] 0.0005068585 0.0005011100 0.0006635865 0.0006088488 0.0005266163
    ##  [6] 0.0004803086 0.0005610908 0.0005671853 0.0005671853 0.0005671853
    ## [11] 0.0005671853 0.0005671853
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01769161 2.938481 
    ##   - best initial criterion value(s) :  14.33926 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -14.339  |proj g|=       2.6581
    ## At iterate     1  f =      -14.533  |proj g|=         1.666
    ## At iterate     2  f =      -14.709  |proj g|=       0.74114
    ## At iterate     3  f =      -14.763  |proj g|=        0.8925
    ## At iterate     4  f =      -14.782  |proj g|=       0.26788
    ## At iterate     5  f =      -14.784  |proj g|=        0.1078
    ## At iterate     6  f =      -14.785  |proj g|=       0.11896
    ## At iterate     7  f =      -14.785  |proj g|=      0.052008
    ## At iterate     8  f =      -14.785  |proj g|=    0.00044003
    ## At iterate     9  f =      -14.785  |proj g|=    1.2917e-05
    ## 
    ## iterations 9
    ## function evaluations 11
    ## segments explored during Cauchy searches 12
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 1.29169e-05
    ## final function value -14.7855
    ## 
    ## F = -14.7855
    ## final  value -14.785473 
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
    ##  [1] 0.0001547161 0.0003749209 0.0003344754 0.0003436513 0.0003158729
    ##  [6] 0.0003651903 0.0003469849 0.0003392235 0.0003392235 0.0003392235
    ## [11] 0.0003392235 0.0003392235
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.016901 2.952755 
    ##   - best initial criterion value(s) :  13.23158 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -13.232  |proj g|=       2.3931
    ## At iterate     1  f =      -13.908  |proj g|=        2.8103
    ## At iterate     2  f =       -14.19  |proj g|=        1.4192
    ## At iterate     3  f =      -14.768  |proj g|=        2.1177
    ## At iterate     4  f =      -15.725  |proj g|=       0.61822
    ## At iterate     5  f =      -15.745  |proj g|=       0.15007
    ## At iterate     6  f =      -15.747  |proj g|=       0.10511
    ## At iterate     7  f =      -15.747  |proj g|=      0.015848
    ## At iterate     8  f =      -15.747  |proj g|=    0.00066494
    ## At iterate     9  f =      -15.747  |proj g|=    7.2465e-06
    ## 
    ## iterations 9
    ## function evaluations 13
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 7.24652e-06
    ## final function value -15.7474
    ## 
    ## F = -15.7474
    ## final  value -15.747353 
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
    ##  [1] 0.0002824318 0.0005011100 0.0006635865 0.0006088488 0.0005266163
    ##  [6] 0.0004803086 0.0005610908 0.0005671853 0.0005671853 0.0005671853
    ## [11] 0.0005671853 0.0005671853
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01790448 2.97088 
    ##   - best initial criterion value(s) :  12.58215 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -12.582  |proj g|=       2.7031
    ## At iterate     1  f =       -13.41  |proj g|=        1.0326
    ## At iterate     2  f =      -14.009  |proj g|=       0.85273
    ## At iterate     3  f =      -14.832  |proj g|=       0.66646
    ## At iterate     4  f =      -15.064  |proj g|=        1.3996
    ## At iterate     5  f =      -15.282  |proj g|=       0.62242
    ## At iterate     6  f =      -15.405  |proj g|=       0.96307
    ## At iterate     7  f =      -15.546  |proj g|=       0.58379
    ## At iterate     8  f =      -15.583  |proj g|=       0.11342
    ## At iterate     9  f =      -15.585  |proj g|=      0.035065
    ## At iterate    10  f =      -15.585  |proj g|=       0.02012
    ## At iterate    11  f =      -15.585  |proj g|=    0.00047863
    ## At iterate    12  f =      -15.585  |proj g|=    4.8212e-06
    ## 
    ## iterations 12
    ## function evaluations 16
    ## segments explored during Cauchy searches 14
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 4.82119e-06
    ## final function value -15.585
    ## 
    ## F = -15.585
    ## final  value -15.584953 
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
    ##  [1] 0.0003428041 0.0001547161 0.0003749209 0.0003344754 0.0003436513
    ##  [6] 0.0003158729 0.0003651903 0.0003469849 0.0003392235 0.0003392235
    ## [11] 0.0003392235 0.0003392235 0.0003392235
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01738981 2.991609 
    ##   - best initial criterion value(s) :  16.90315 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -16.903  |proj g|=        2.406
    ## At iterate     1  f =      -18.463  |proj g|=       0.58954
    ## At iterate     2  f =      -18.538  |proj g|=       0.96635
    ## At iterate     3  f =      -18.592  |proj g|=        1.5454
    ## At iterate     4  f =      -18.603  |proj g|=       0.61169
    ## At iterate     5  f =      -18.605  |proj g|=      0.095302
    ## At iterate     6  f =      -18.606  |proj g|=      0.038965
    ## At iterate     7  f =      -18.606  |proj g|=      0.043338
    ## At iterate     8  f =      -18.606  |proj g|=      0.025708
    ## At iterate     9  f =      -18.606  |proj g|=     0.0037794
    ## At iterate    10  f =      -18.606  |proj g|=    0.00030972
    ## At iterate    11  f =      -18.606  |proj g|=     6.609e-06
    ## 
    ## iterations 11
    ## function evaluations 14
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 6.60905e-06
    ## final function value -18.6059
    ## 
    ## F = -18.6059
    ## final  value -18.605857 
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
    ##  [1] 0.0004349449 0.0002824318 0.0005011100 0.0006635865 0.0006088488
    ##  [6] 0.0005266163 0.0004803086 0.0005610908 0.0005671853 0.0005671853
    ## [11] 0.0005671853 0.0005671853 0.0005671853
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01816458 3.079909 
    ##   - best initial criterion value(s) :  17.28535 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -17.285  |proj g|=      0.45393
    ## At iterate     1  f =      -17.321  |proj g|=       0.34272
    ## At iterate     2  f =      -17.443  |proj g|=       0.96557
    ## At iterate     3  f =      -17.523  |proj g|=       0.17588
    ## At iterate     4  f =      -17.525  |proj g|=       0.13341
    ## At iterate     5  f =      -17.527  |proj g|=      0.017297
    ## At iterate     6  f =      -17.527  |proj g|=    0.00096585
    ## At iterate     7  f =      -17.527  |proj g|=    0.00014634
    ## 
    ## iterations 7
    ## function evaluations 11
    ## segments explored during Cauchy searches 9
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 0.000146338
    ## final function value -17.5271
    ## 
    ## F = -17.5271
    ## final  value -17.527070 
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

![](README_files/figure-gfm/unnamed-chunk-25-1.png)<!-- --> New design
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

![](README_files/figure-gfm/unnamed-chunk-26-1.png)<!-- --> Simple plot.

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

![](README_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->
