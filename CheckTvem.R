rm(list = ls());
library(mgcv);
library(refund);
library(boot);
source("C:\\Users\\JJD\\Documents\\tvem\\R\\simulate_tvem_example.R"); 
source("C:\\Users\\JJD\\Documents\\tvem\\R\\tvem.r"); 
source("C:\\Users\\JJD\\Documents\\tvem\\R\\print_tvem.r"); 
source("C:\\Users\\JJD\\Documents\\tvem\\R\\plot_tvem.r");   
set.seed(3172020);
# Set true parameters for data-generating model:
#     Sample size parameters:
nsim <- 100;   # number of simulations to run
nsub_values <- c(100,300);          # number of subjects
ntimes <- 60;   # number of potential times that could be observed
nboot <- 20;  # number of bootstrapped samples to generate per original sample
observe_rate <- .9;   # proportion of potential times on which there are actually observations;
start_time <- Sys.time();
# by design to reduce participant burden.
#     Timeline:
time_grid <- seq(-1,1,length=ntimes);  # vector of all possible times;
#     True coefficients:
gamma_int <- function(t) {return(1+t)}; # time-varying mean of outcome variable for the X=0 group
gamma_X <- function(t) {return(rep(1,length(t)))}; # time-varying effect of X on the outcome 
#gamma_X <- function(t) {return(t/2)}; # time-varying effect of X on the outcome 

sigma_y <- 1;
rho_y <- .5; 
answers <- NULL;

start_time <- Sys.time();
for (nsub in nsub_values) {
  for (this_sim in 1:nsim) {
    short_X <- rbinom(n=nsub, size=1, prob=.5);  # binary treatment assignment;
    # Simulate Y from X...
    autoreg_error <- matrix(0,nsub,ntimes);
    autoreg_error[,1] <- rnorm(n=nsub,mean=0,sd=sigma_y);
    for (j in 2:ntimes) {
      autoreg_error[,j] <- rho_y*autoreg_error[,j-1] +
        sqrt(1-rho_y^2)*rnorm(n=nsub,mean=0,sd=sigma_y);
    }
    all_Y <- matrix(0,nsub,ntimes);  # time-varying mediator; 
    for (i in 1:nsub) { 
      all_Y[i,] <- gamma_int(time_grid) + short_X[i]*gamma_X(time_grid);
    } 
    all_Y <- all_Y + autoreg_error;
    # Assemble random data into a long-form dataset:
    short_ID <- 1:nsub;
    long_data <- data.frame(
      ID = rep(short_ID, each=ntimes),
      time = rep(time_grid, times=nsub),
      X = rep(short_X, each=ntimes),
      Y = as.vector(t(all_Y)) );
    ans1 <- tvem(data=long_data,
                 formula=Y~X,
                 id=ID,
                 time=time,
                 family=gaussian(), 
                 num_knots=20);
    gamma_int_residuals <- ans1$grid_fitted_coefficients$`(Intercept)`[,"estimate"] - gamma_int(ans1$time_grid);
    gamma_int_se <- ans1$grid_fitted_coefficients$`(Intercept)`[,"standard_error"];
    gamma_X_residuals <- ans1$grid_fitted_coefficients$`X`[,"estimate"] - gamma_X(ans1$time_grid);
    gamma_X_se <- ans1$grid_fitted_coefficients$`X`[,"standard_error"];
    tvem_model_summary <- summary(ans1$back_end_model);
    current_time <- Sys.time();
    answers <- rbind(answers, 
                     # Create structures to hold results:  
                     # ... for the effect of X on M: 
                     c(nsub=nsub,
                       ntimes=ntimes, 
                       this_sim=this_sim,
                       gamma_int_estimate_bias =  mean(gamma_int_residuals), 
                       gamma_int_estimate_mse =  mean(gamma_int_residuals^2),
                       gamma_int_estimate_mean_std_err = mean(gamma_int_se),
                       gamma_int_estimate_coverage =  mean(abs(gamma_int_residuals)<1.96*gamma_int_se), 
                       gamma_int_estimate_family_coverage =  all(abs(gamma_int_residuals)<1.96*gamma_int_se), 
                       gamma_int_estimate_pvalue_time0 =  as.numeric(tvem_model_summary$p.pv["(Intercept)"]), 
                       gamma_int_estimate_pvalue_varying =  tvem_model_summary$s.table["s(time)","p-value"],
                       gamma_X_estimate_bias =  mean(gamma_X_residuals), 
                       gamma_X_estimate_mse =  mean(gamma_X_residuals^2),
                       gamma_X_estimate_mean_std_err = mean(gamma_X_se),
                       gamma_X_estimate_coverage =  mean(abs(gamma_X_residuals)<1.96*gamma_X_se), 
                       gamma_X_estimate_family_coverage =  all(abs(gamma_X_residuals)<1.96*gamma_X_se), 
                       gamma_X_estimate_pvalue_time0 =  as.numeric(tvem_model_summary$p.pv["X"]),   
                       gamma_X_estimate_pvalue_varying =  tvem_model_summary$s.table["s(time):X","p-value"]
                     ));
  }
}
print(apply(answers,2,mean));
time_required <- difftime(current_time,start_time);
print(time_required);
mean_answers_100 <- apply(answers[which(answers[,"nsub"]==100),],2,mean);
mean_answers_300 <- apply(answers[which(answers[,"nsub"]==300),],2,mean);
print(round(cbind(mean_answers_100,
                  mean_answers_300),3));

