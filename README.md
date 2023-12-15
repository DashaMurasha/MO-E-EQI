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

    ## Skipping install of 'MOEEQI' from a github remote, the SHA1 (38912f5f) has not changed since last install.
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
    ## [1] 0.0003465255 0.0003465255 0.0003465255 0.0003465255 0.0003465255
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01778491 2.631419 
    ##   - best initial criterion value(s) :  0.1742082 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=     -0.17421  |proj g|=      0.65046
    ## At iterate     1  f =     -0.24888  |proj g|=       0.27877
    ## At iterate     2  f =     -0.29213  |proj g|=       0.34138
    ## At iterate     3  f =     -0.31354  |proj g|=       0.14961
    ## At iterate     4  f =     -0.31564  |proj g|=      0.099504
    ## At iterate     5  f =     -0.31578  |proj g|=      0.021781
    ## At iterate     6  f =     -0.31579  |proj g|=    0.00052609
    ## At iterate     7  f =     -0.31579  |proj g|=    4.0478e-05
    ## At iterate     8  f =     -0.31579  |proj g|=    3.4402e-06
    ## 
    ## iterations 8
    ## function evaluations 11
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 3.4402e-06
    ## final function value -0.315789
    ## 
    ## F = -0.315789
    ## final  value -0.315789 
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
    ## [1] 0.0006306149 0.0006306149 0.0006306149 0.0006306149 0.0006306149
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01757449 2.736529 
    ##   - best initial criterion value(s) :  -1.334428 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=       1.3344  |proj g|=       1.5772
    ## At iterate     1  f =      0.99062  |proj g|=        1.0013
    ## At iterate     2  f =      0.94046  |proj g|=       0.90683
    ## At iterate     3  f =      0.87429  |proj g|=        1.9684
    ## At iterate     4  f =      0.78332  |proj g|=        1.1875
    ## At iterate     5  f =       0.6944  |proj g|=        1.0999
    ## At iterate     6  f =      0.68852  |proj g|=       0.80636
    ## At iterate     7  f =      0.65306  |proj g|=       0.31758
    ## At iterate     8  f =      0.64889  |proj g|=       0.21871
    ## At iterate     9  f =      0.64796  |proj g|=      0.085896
    ## At iterate    10  f =       0.6475  |proj g|=      0.010361
    ## At iterate    11  f =      0.64749  |proj g|=    0.00060976
    ## At iterate    12  f =      0.64749  |proj g|=    4.6704e-05
    ## At iterate    13  f =      0.64749  |proj g|=    4.9223e-07
    ## 
    ## iterations 13
    ## function evaluations 16
    ## segments explored during Cauchy searches 14
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 4.92231e-07
    ## final function value 0.647488
    ## 
    ## F = 0.647488
    ## final  value 0.647488 
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
    ## [1] 0.0003466258 0.0003465255 0.0003465255 0.0003465255 0.0003465255
    ## [6] 0.0003465255
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01492895 2.290179 
    ##   - best initial criterion value(s) :  1.767922 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -1.7679  |proj g|=       2.1665
    ## At iterate     1  f =      -2.2036  |proj g|=        1.4007
    ## At iterate     2  f =       -2.285  |proj g|=       0.88562
    ## At iterate     3  f =      -2.3682  |proj g|=       0.23705
    ## At iterate     4  f =      -2.3939  |proj g|=       0.97373
    ## At iterate     5  f =      -2.4396  |proj g|=         1.769
    ## At iterate     6  f =      -2.5497  |proj g|=        1.4539
    ## At iterate     7  f =      -2.5653  |proj g|=       0.45802
    ## At iterate     8  f =      -2.5687  |proj g|=       0.08385
    ## At iterate     9  f =      -2.5691  |proj g|=      0.048199
    ## At iterate    10  f =      -2.5692  |proj g|=      0.017827
    ## At iterate    11  f =      -2.5692  |proj g|=    0.00085927
    ## At iterate    12  f =      -2.5692  |proj g|=    4.8409e-05
    ## At iterate    13  f =      -2.5692  |proj g|=    1.0946e-06
    ## 
    ## iterations 13
    ## function evaluations 17
    ## segments explored during Cauchy searches 15
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 1.09461e-06
    ## final function value -2.56918
    ## 
    ## F = -2.56918
    ## final  value -2.569185 
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
    ## [1] 0.0005914674 0.0006306149 0.0006306149 0.0006306149 0.0006306149
    ## [6] 0.0006306149
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01706911 2.546833 
    ##   - best initial criterion value(s) :  0.324634 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=     -0.32463  |proj g|=       2.2986
    ## At iterate     1  f =     -0.41237  |proj g|=       0.55003
    ## At iterate     2  f =     -0.41993  |proj g|=       0.37697
    ## At iterate     3  f =     -0.42725  |proj g|=       0.25366
    ## At iterate     4  f =     -0.43179  |proj g|=       0.26953
    ## At iterate     5  f =     -0.43506  |proj g|=      0.069095
    ## At iterate     6  f =     -0.43518  |proj g|=      0.010694
    ## At iterate     7  f =     -0.43519  |proj g|=    0.00062413
    ## At iterate     8  f =     -0.43519  |proj g|=    0.00062746
    ## At iterate     9  f =     -0.43519  |proj g|=    0.00091836
    ## At iterate    10  f =     -0.43519  |proj g|=     0.0017356
    ## At iterate    11  f =     -0.43519  |proj g|=     0.0019927
    ## At iterate    12  f =     -0.43519  |proj g|=     0.0012549
    ## At iterate    13  f =     -0.43519  |proj g|=    0.00021478
    ## At iterate    14  f =     -0.43519  |proj g|=    1.9895e-06
    ## 
    ## iterations 14
    ## function evaluations 17
    ## segments explored during Cauchy searches 15
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 0
    ## norm of the final projected gradient 1.98949e-06
    ## final function value -0.435187
    ## 
    ## F = -0.435187
    ## final  value -0.435187 
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
    ## [1] 0.0003171991 0.0003466258 0.0003465255 0.0003465255 0.0003465255
    ## [6] 0.0003465255 0.0003465255
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01414678 2.231315 
    ##   - best initial criterion value(s) :  3.827817 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -3.8278  |proj g|=       1.7318
    ## At iterate     1  f =      -3.9971  |proj g|=        2.0761
    ## At iterate     2  f =      -4.1477  |proj g|=        2.0439
    ## At iterate     3  f =      -4.2912  |proj g|=       0.26884
    ## At iterate     4  f =      -4.2968  |proj g|=       0.50737
    ## At iterate     5  f =      -4.2986  |proj g|=       0.30689
    ## At iterate     6  f =      -4.2997  |proj g|=      0.031768
    ## At iterate     7  f =      -4.2997  |proj g|=     0.0024857
    ## At iterate     8  f =      -4.2997  |proj g|=    2.1177e-05
    ## At iterate     9  f =      -4.2997  |proj g|=     5.754e-08
    ## 
    ## iterations 9
    ## function evaluations 12
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 5.75404e-08
    ## final function value -4.29973
    ## 
    ## F = -4.29973
    ## final  value -4.299729 
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
    ## [1] 0.0006506923 0.0005914674 0.0006306149 0.0006306149 0.0006306149
    ## [6] 0.0006306149 0.0006306149
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01436855 2.279124 
    ##   - best initial criterion value(s) :  1.552966 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=       -1.553  |proj g|=       2.4562
    ## At iterate     1  f =      -1.7392  |proj g|=        1.0263
    ## At iterate     2  f =      -1.7876  |proj g|=       0.71079
    ## At iterate     3  f =      -1.8057  |proj g|=       0.58124
    ## At iterate     4  f =       -1.829  |proj g|=       0.38903
    ## At iterate     5  f =      -1.8605  |proj g|=       0.14429
    ## At iterate     6  f =      -1.8663  |proj g|=        0.2301
    ## At iterate     7  f =      -1.8672  |proj g|=      0.058375
    ## At iterate     8  f =      -1.8672  |proj g|=     0.0042499
    ## At iterate     9  f =      -1.8672  |proj g|=     0.0031347
    ## At iterate    10  f =      -1.8672  |proj g|=    0.00012571
    ## At iterate    11  f =      -1.8672  |proj g|=    3.0192e-05
    ## 
    ## iterations 11
    ## function evaluations 14
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 0
    ## norm of the final projected gradient 3.01923e-05
    ## final function value -1.86718
    ## 
    ## F = -1.86718
    ## final  value -1.867177 
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
    ## [1] 0.0003549355 0.0003171991 0.0003466258 0.0003465255 0.0003465255
    ## [6] 0.0003465255 0.0003465255 0.0003465255
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01577166 2.446084 
    ##   - best initial criterion value(s) :  5.083944 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -5.0839  |proj g|=       2.5238
    ## At iterate     1  f =      -5.5511  |proj g|=        2.3491
    ## At iterate     2  f =      -6.0925  |proj g|=        2.2912
    ## At iterate     3  f =      -6.2638  |proj g|=        2.0857
    ## At iterate     4  f =      -6.5492  |proj g|=        2.0713
    ## At iterate     5  f =      -6.6179  |proj g|=        1.2999
    ## At iterate     6  f =      -6.6345  |proj g|=       0.25362
    ## At iterate     7  f =       -6.635  |proj g|=       0.12412
    ## At iterate     8  f =      -6.6354  |proj g|=      0.023846
    ## At iterate     9  f =      -6.6355  |proj g|=     0.0084713
    ## At iterate    10  f =      -6.6355  |proj g|=    0.00075284
    ## At iterate    11  f =      -6.6355  |proj g|=    3.1644e-05
    ## 
    ## iterations 11
    ## function evaluations 14
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 3.16443e-05
    ## final function value -6.63547
    ## 
    ## F = -6.63547
    ## final  value -6.635474 
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
    ## [1] 0.0006048420 0.0006506923 0.0005914674 0.0006306149 0.0006306149
    ## [6] 0.0006306149 0.0006306149 0.0006306149
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.0150663 2.46487 
    ##   - best initial criterion value(s) :  2.444147 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -2.4441  |proj g|=       2.3332
    ## At iterate     1  f =      -2.5014  |proj g|=         1.342
    ## At iterate     2  f =      -2.5128  |proj g|=        1.0319
    ## At iterate     3  f =      -2.5442  |proj g|=       0.44348
    ## At iterate     4  f =        -2.56  |proj g|=       0.41143
    ## At iterate     5  f =      -2.6071  |proj g|=        0.4081
    ## At iterate     6  f =      -2.6794  |proj g|=       0.56307
    ## At iterate     7  f =      -2.6817  |proj g|=      0.098852
    ## At iterate     8  f =      -2.6879  |proj g|=       0.70924
    ## At iterate     9  f =      -2.6905  |proj g|=       0.25966
    ## At iterate    10  f =      -2.6914  |proj g|=      0.048858
    ## At iterate    11  f =      -2.6914  |proj g|=     0.0020923
    ## At iterate    12  f =      -2.6914  |proj g|=    0.00014461
    ## At iterate    13  f =      -2.6914  |proj g|=    8.1044e-06
    ## 
    ## iterations 13
    ## function evaluations 21
    ## segments explored during Cauchy searches 14
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 8.10436e-06
    ## final function value -2.69143
    ## 
    ## F = -2.69143
    ## final  value -2.691429 
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
    ## [1] 0.0003747903 0.0003549355 0.0003171991 0.0003466258 0.0003465255
    ## [6] 0.0003465255 0.0003465255 0.0003465255 0.0003465255
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.0168879 2.773833 
    ##   - best initial criterion value(s) :  7.959529 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -7.9595  |proj g|=       2.4955
    ## At iterate     1  f =      -8.1333  |proj g|=       0.80239
    ## At iterate     2  f =      -8.7975  |proj g|=       0.34093
    ## At iterate     3  f =       -8.936  |proj g|=       0.55523
    ## At iterate     4  f =       -8.939  |proj g|=      0.093752
    ## At iterate     5  f =      -8.9391  |proj g|=      0.015207
    ## At iterate     6  f =      -8.9391  |proj g|=     0.0024865
    ## At iterate     7  f =      -8.9391  |proj g|=    0.00093865
    ## At iterate     8  f =      -8.9391  |proj g|=    0.00021093
    ## 
    ## iterations 8
    ## function evaluations 12
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 0.000210927
    ## final function value -8.93912
    ## 
    ## F = -8.93912
    ## final  value -8.939119 
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
    ## [1] 0.0006837032 0.0006048420 0.0006506923 0.0005914674 0.0006306149
    ## [6] 0.0006306149 0.0006306149 0.0006306149 0.0006306149
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01557586 2.523909 
    ##   - best initial criterion value(s) :  4.590515 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -4.5905  |proj g|=       2.2339
    ## At iterate     1  f =      -4.6717  |proj g|=       0.20456
    ## At iterate     2  f =      -4.6791  |proj g|=       0.10487
    ## At iterate     3  f =       -4.687  |proj g|=       0.11816
    ## At iterate     4  f =      -4.6872  |proj g|=       0.08476
    ## At iterate     5  f =      -4.6876  |proj g|=      0.010255
    ## At iterate     6  f =      -4.6876  |proj g|=     0.0041831
    ## At iterate     7  f =      -4.6876  |proj g|=    0.00030126
    ## At iterate     8  f =      -4.6876  |proj g|=    2.6117e-05
    ## 
    ## iterations 8
    ## function evaluations 11
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 2.6117e-05
    ## final function value -4.68762
    ## 
    ## F = -4.68762
    ## final  value -4.687618 
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
    ##  [1] 0.0003441376 0.0003747903 0.0003549355 0.0003171991 0.0003466258
    ##  [6] 0.0003465255 0.0003465255 0.0003465255 0.0003465255 0.0003465255
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01704668 2.955174 
    ##   - best initial criterion value(s) :  9.731737 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -9.7317  |proj g|=       2.0995
    ## At iterate     1  f =      -10.779  |proj g|=       0.42688
    ## At iterate     2  f =      -11.091  |proj g|=       0.85378
    ## At iterate     3  f =      -11.599  |proj g|=       0.86377
    ## At iterate     4  f =      -11.634  |proj g|=       0.10413
    ## At iterate     5  f =      -11.634  |proj g|=      0.044755
    ## At iterate     6  f =      -11.634  |proj g|=      0.023238
    ## At iterate     7  f =      -11.634  |proj g|=     0.0007049
    ## At iterate     8  f =      -11.634  |proj g|=    2.2245e-05
    ## 
    ## iterations 8
    ## function evaluations 11
    ## segments explored during Cauchy searches 10
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 2.22452e-05
    ## final function value -11.6342
    ## 
    ## F = -11.6342
    ## final  value -11.634161 
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
    ##  [1] 0.0005727601 0.0006837032 0.0006048420 0.0006506923 0.0005914674
    ##  [6] 0.0006306149 0.0006306149 0.0006306149 0.0006306149 0.0006306149
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01703123 2.782714 
    ##   - best initial criterion value(s) :  6.575822 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -6.5758  |proj g|=       2.1878
    ## At iterate     1  f =      -6.9297  |proj g|=        2.4018
    ## At iterate     2  f =      -7.0791  |proj g|=       0.75079
    ## At iterate     3  f =       -7.118  |proj g|=       0.32521
    ## At iterate     4  f =      -7.2267  |proj g|=       0.36389
    ## At iterate     5  f =      -7.2292  |proj g|=       0.23487
    ## At iterate     6  f =       -7.231  |proj g|=      0.059956
    ## At iterate     7  f =      -7.2314  |proj g|=      0.027118
    ## At iterate     8  f =      -7.2314  |proj g|=     0.0014011
    ## At iterate     9  f =      -7.2314  |proj g|=    1.6852e-05
    ## At iterate    10  f =      -7.2314  |proj g|=    8.4175e-08
    ## 
    ## iterations 10
    ## function evaluations 13
    ## segments explored during Cauchy searches 12
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 8.41746e-08
    ## final function value -7.23141
    ## 
    ## F = -7.23141
    ## final  value -7.231408 
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
    ##  [1] 0.0003068717 0.0003441376 0.0003747903 0.0003549355 0.0003171991
    ##  [6] 0.0003466258 0.0003465255 0.0003465255 0.0003465255 0.0003465255
    ## [11] 0.0003465255
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01746465 3.033991 
    ##   - best initial criterion value(s) :  12.29924 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -12.299  |proj g|=       2.7534
    ## At iterate     1  f =      -12.558  |proj g|=        1.1641
    ## At iterate     2  f =       -13.46  |proj g|=       0.49121
    ## At iterate     3  f =      -13.645  |proj g|=         2.854
    ## At iterate     4  f =        -13.7  |proj g|=        0.2917
    ## At iterate     5  f =      -13.706  |proj g|=        0.2164
    ## At iterate     6  f =      -13.706  |proj g|=      0.033105
    ## At iterate     7  f =      -13.706  |proj g|=      0.002545
    ## At iterate     8  f =      -13.706  |proj g|=     0.0014456
    ## At iterate     9  f =      -13.706  |proj g|=    4.2891e-05
    ## At iterate    10  f =      -13.706  |proj g|=    3.9495e-07
    ## 
    ## iterations 10
    ## function evaluations 14
    ## segments explored during Cauchy searches 13
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 3.94952e-07
    ## final function value -13.7064
    ## 
    ## F = -13.7064
    ## final  value -13.706402 
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
    ##  [1] 0.0006020274 0.0005727601 0.0006837032 0.0006048420 0.0006506923
    ##  [6] 0.0005914674 0.0006306149 0.0006306149 0.0006306149 0.0006306149
    ## [11] 0.0006306149
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01787253 2.914783 
    ##   - best initial criterion value(s) :  8.416896 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -8.4169  |proj g|=       2.2287
    ## At iterate     1  f =      -9.1598  |proj g|=        1.5822
    ## At iterate     2  f =      -9.2604  |proj g|=       0.41486
    ## At iterate     3  f =      -9.3217  |proj g|=       0.35139
    ## At iterate     4  f =      -9.4818  |proj g|=       0.27411
    ## At iterate     5  f =      -9.4846  |proj g|=        0.1869
    ## At iterate     6  f =      -9.4866  |proj g|=      0.039885
    ## At iterate     7  f =      -9.4867  |proj g|=     0.0095408
    ## At iterate     8  f =      -9.4867  |proj g|=    0.00098825
    ## At iterate     9  f =      -9.4867  |proj g|=    4.7878e-05
    ## At iterate    10  f =      -9.4867  |proj g|=    1.0162e-05
    ## 
    ## iterations 10
    ## function evaluations 13
    ## segments explored during Cauchy searches 12
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 1.01619e-05
    ## final function value -9.48667
    ## 
    ## F = -9.48667
    ## final  value -9.486669 
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
    ##  [1] 0.0003978713 0.0003068717 0.0003441376 0.0003747903 0.0003549355
    ##  [6] 0.0003171991 0.0003466258 0.0003465255 0.0003465255 0.0003465255
    ## [11] 0.0003465255 0.0003465255
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01840002 3.168734 
    ##   - best initial criterion value(s) :  14.39234 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -14.392  |proj g|=        2.925
    ## At iterate     1  f =      -14.551  |proj g|=       0.68376
    ## At iterate     2  f =      -15.616  |proj g|=       0.60068
    ## At iterate     3  f =      -16.253  |proj g|=        0.6221
    ## At iterate     4  f =      -16.268  |proj g|=       0.12894
    ## At iterate     5  f =      -16.269  |proj g|=       0.10083
    ## At iterate     6  f =       -16.27  |proj g|=     0.0073907
    ## At iterate     7  f =       -16.27  |proj g|=     0.0004362
    ## At iterate     8  f =       -16.27  |proj g|=    5.5326e-07
    ## 
    ## iterations 8
    ## function evaluations 12
    ## segments explored during Cauchy searches 10
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 5.53259e-07
    ## final function value -16.2701
    ## 
    ## F = -16.2701
    ## final  value -16.270108 
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
    ##  [1] 0.0004948190 0.0006020274 0.0005727601 0.0006837032 0.0006048420
    ##  [6] 0.0006506923 0.0005914674 0.0006306149 0.0006306149 0.0006306149
    ## [11] 0.0006306149 0.0006306149
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01814953 3.027089 
    ##   - best initial criterion value(s) :  12.06435 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -12.064  |proj g|=      0.29888
    ## At iterate     1  f =      -12.068  |proj g|=      0.050639
    ## At iterate     2  f =      -12.088  |proj g|=       0.12274
    ## At iterate     3  f =      -12.088  |proj g|=      0.024712
    ## At iterate     4  f =      -12.088  |proj g|=     0.0018823
    ## At iterate     5  f =      -12.088  |proj g|=    3.1579e-05
    ## At iterate     6  f =      -12.088  |proj g|=    1.7995e-06
    ## 
    ## iterations 6
    ## function evaluations 10
    ## segments explored during Cauchy searches 7
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 1.79949e-06
    ## final function value -12.0883
    ## 
    ## F = -12.0883
    ## final  value -12.088276 
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
    ##  [1] 0.0003260074 0.0003978713 0.0003068717 0.0003441376 0.0003747903
    ##  [6] 0.0003549355 0.0003171991 0.0003466258 0.0003465255 0.0003465255
    ## [11] 0.0003465255 0.0003465255 0.0003465255
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01871067 3.289407 
    ##   - best initial criterion value(s) :  18.38135 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -18.381  |proj g|=       1.6956
    ## At iterate     1  f =      -18.671  |proj g|=         1.061
    ## At iterate     2  f =      -18.743  |proj g|=       0.28049
    ## At iterate     3  f =      -18.793  |proj g|=       0.21727
    ## At iterate     4  f =      -18.794  |proj g|=      0.034373
    ## At iterate     5  f =      -18.794  |proj g|=      0.004999
    ## At iterate     6  f =      -18.794  |proj g|=    8.4265e-05
    ## At iterate     7  f =      -18.794  |proj g|=    5.2069e-06
    ## 
    ## iterations 7
    ## function evaluations 10
    ## segments explored during Cauchy searches 10
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 5.20693e-06
    ## final function value -18.7943
    ## 
    ## F = -18.7943
    ## final  value -18.794298 
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
    ##  [1] 0.0005329318 0.0004948190 0.0006020274 0.0005727601 0.0006837032
    ##  [6] 0.0006048420 0.0006506923 0.0005914674 0.0006306149 0.0006306149
    ## [11] 0.0006306149 0.0006306149 0.0006306149
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01904868 3.191344 
    ##   - best initial criterion value(s) :  13.76681 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -13.767  |proj g|=       1.4121
    ## At iterate     1  f =      -13.915  |proj g|=       0.32913
    ## At iterate     2  f =          -14  |proj g|=       0.39476
    ## At iterate     3  f =      -14.222  |proj g|=       0.95548
    ## At iterate     4  f =      -14.245  |proj g|=       0.08847
    ## At iterate     5  f =      -14.246  |proj g|=      0.041112
    ## At iterate     6  f =      -14.246  |proj g|=     0.0055195
    ## At iterate     7  f =      -14.246  |proj g|=    0.00014068
    ## At iterate     8  f =      -14.246  |proj g|=    5.3938e-06
    ## 
    ## iterations 8
    ## function evaluations 11
    ## segments explored during Cauchy searches 11
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 5.39382e-06
    ## final function value -14.2457
    ## 
    ## F = -14.2457
    ## final  value -14.245715 
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
    ##  [1] 0.0003260074 0.0003978713 0.0003068717 0.0001682426 0.0003747903
    ##  [6] 0.0003549355 0.0003171991 0.0003466258 0.0003465255 0.0003465255
    ## [11] 0.0003465255 0.0003465255 0.0003465255
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01862075 3.274728 
    ##   - best initial criterion value(s) :  16.91969 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=       -16.92  |proj g|=       2.9871
    ## At iterate     1  f =      -17.441  |proj g|=        1.2204
    ## At iterate     2  f =       -17.85  |proj g|=       0.58096
    ## At iterate     3  f =      -18.288  |proj g|=        1.1959
    ## At iterate     4  f =      -18.391  |proj g|=        1.5092
    ## At iterate     5  f =      -18.548  |proj g|=         1.436
    ## At iterate     6  f =       -18.66  |proj g|=       0.23145
    ## At iterate     7  f =      -18.662  |proj g|=       0.48897
    ## At iterate     8  f =      -18.664  |proj g|=      0.065752
    ## At iterate     9  f =      -18.664  |proj g|=     0.0064695
    ## At iterate    10  f =      -18.664  |proj g|=    0.00010384
    ## At iterate    11  f =      -18.664  |proj g|=    2.6418e-07
    ## 
    ## iterations 11
    ## function evaluations 15
    ## segments explored during Cauchy searches 16
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 2.64183e-07
    ## final function value -18.6644
    ## 
    ## F = -18.6644
    ## final  value -18.664394 
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
    ##  [1] 0.0005329318 0.0004948190 0.0006020274 0.0002987524 0.0006837032
    ##  [6] 0.0006048420 0.0006506923 0.0005914674 0.0006306149 0.0006306149
    ## [11] 0.0006306149 0.0006306149 0.0006306149
    ##   - parameters lower bounds :  1e-10 1e-10 
    ##   - parameters upper bounds :  3.141593 2 
    ##   - variance bounds :  0.01912528 3.205168 
    ##   - best initial criterion value(s) :  14.19667 
    ## 
    ## N = 3, M = 5 machine precision = 2.22045e-16
    ## At X0, 0 variables are exactly at the bounds
    ## At iterate     0  f=      -14.197  |proj g|=        2.962
    ## At iterate     1  f =      -14.388  |proj g|=       0.31041
    ## At iterate     2  f =      -14.454  |proj g|=       0.77112
    ## At iterate     3  f =      -14.518  |proj g|=       0.82504
    ## At iterate     4  f =      -14.529  |proj g|=       0.19707
    ## At iterate     5  f =       -14.53  |proj g|=      0.069533
    ## At iterate     6  f =       -14.53  |proj g|=    0.00063976
    ## At iterate     7  f =       -14.53  |proj g|=    2.2464e-05
    ## 
    ## iterations 7
    ## function evaluations 11
    ## segments explored during Cauchy searches 9
    ## BFGS updates skipped 0
    ## active bounds at final generalized Cauchy point 1
    ## norm of the final projected gradient 2.24642e-05
    ## final function value -14.5298
    ## 
    ## F = -14.5298
    ## final  value -14.529763 
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
