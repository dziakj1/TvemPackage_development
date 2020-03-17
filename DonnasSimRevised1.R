rm(list = ls());
library(mgcv);
library(refund);
library(boot);
source("C:\\Users\\JJD\\Documents\\tvem\\R\\simulate_tvem_example.R"); 
source("C:\\Users\\JJD\\Documents\\tvem\\R\\tvem.r"); 
source("C:\\Users\\JJD\\Documents\\tvem\\R\\print_tvem.r"); 
source("C:\\Users\\JJD\\Documents\\tvem\\R\\plot_tvem.r");  
source("C:\\Users\\JJD\\Documents\\tvem\\R\\funreg_mediation.r");
source("C:\\Users\\JJD\\Documents\\tvem\\R\\print_funreg_mediation.r");
source("C:\\Users\\JJD\\Documents\\tvem\\R\\plot_funreg_mediation.r");
set.seed(32020);
# Set true parameters for data-generating model:
#     Sample size parameters:
nsim <- 20;   # number of simulations to run
nsub_values <- c(100,300);          # number of subjects
ntimes <- 60;   # number of potential times that could be observed
nboot <- 20;  # number of bootstrapped samples to generate per original sample
observe_rate <- .9;   # proportion of potential times on which there are actually observations;
start_time <- Sys.time();
# by design to reduce participant burden.
#     Timeline:
time_grid <- seq(0,1,length=ntimes);  # vector of all possible times;
#     True coefficients:
gamma_int <- function(t) {return(t^.5)}; # time-varying mean of mediator variable for the X=0 group
gamma_X <- function(t) {return(-(t/2)^.5)}; # time-varying effect of X on the mediator
alpha_M <- function(t) {exp(t)-1}; # functional (funreg) coefficient for cumulative effect of M on Y
alpha_int <- 0;  # mean of Y if the X is zero and M is the 0 function; 
alpha_X <- 0;  # direct effect of X on Y after adjusting for M;
beta_true <- mean(alpha_M(time_grid)*gamma_X(time_grid));
#     Error standard deviations:
sigma_m_intercept <- 1;
sigma_m_slope <- 1;
sigma_m_error <- 1;
rho_m_error <- .8; 
answers <- NULL;
start_time <- Sys.time();
for (nsub in nsub_values) {
for (this_sim in 1:nsim) {
  short_X <- rbinom(n=nsub, size=1, prob=.5);  # binary treatment assignment;
  # Simulate M from X...
  autoreg_error <- matrix(0,nsub,ntimes);
  autoreg_error[,1] <- rnorm(n=nsub,mean=0,sd=sigma_m_error);
  for (j in 2:ntimes) {
    autoreg_error[,j] <- rho_m_error*autoreg_error[,j-1] +
      sqrt(1-rho_m_error^2)*rnorm(n=nsub,mean=0,sd=sigma_m_error);
  }
  all_M <- matrix(0,nsub,ntimes);  # time-varying mediator; 
  for (i in 1:nsub) { 
    all_M[i,] <- gamma_int(time_grid) + short_X[i]*gamma_X(time_grid) +
      rnorm(n=1,mean=0,sd=sigma_m_intercept) + 
       rnorm(n=1,mean=0,sd=sigma_m_slope)*time_grid;
  } 
  all_M <- all_M + autoreg_error;
  # Simulate Y from M and X...
  eta <- rep(NA,nsub); # = E(Y|X,M);
  for (i in 1:nsub) {
    eta[i] <- alpha_int + mean(alpha_M(time_grid) * all_M[i,]) + alpha_X*short_X[i];
  }
  mu <- exp(eta)/(1+exp(eta));
  short_Y <- unlist(lapply(X=mu,FUN=rbinom,size=1,n=1));
  # Assemble random data into a long-form dataset:
  short_ID <- 1:nsub;
  long_data <- data.frame(
    ID = rep(short_ID, each=ntimes),
    time = rep(time_grid, times=nsub),
    X = rep(short_X, each=ntimes),
    M = as.vector(t(all_M)),
    Y = rep(short_Y, each=ntimes));
  fun_med_results <- funreg_mediation(data=long_data,
                                      treatment=X,
                                      mediator=M,
                                      outcome=Y,
                                      id=ID,
                                      time=time,
                                      logistic=TRUE,
                                      nboot=nboot); 
  tvem_model_summary <- summary(fun_med_results$original_results$tvem_XM_details$back_end_model);
  gamma_int_residuals <- fun_med_results$original_results$gamma_int_estimate - gamma_int(time_grid);
  gamma_int_se <- fun_med_results$original_results$gamma_int_se;
  gamma_X_residuals <- fun_med_results$original_results$gamma_X_estimate - gamma_X(time_grid);
  gamma_X_se <- fun_med_results$original_results$gamma_X_se;
  alpha_M_residuals <- fun_med_results$original_results$alpha_M_estimate - alpha_M(time_grid);
  alpha_M_se <- fun_med_results$original_results$alpha_M_se;
  norm_cover <- (fun_med_results$bootstrap_results$beta_boot_norm_lower < beta_true) & 
    (fun_med_results$bootstrap_results$beta_boot_norm_upper > beta_true);
  basic_cover <- (fun_med_results$bootstrap_results$beta_boot_basic_lower < beta_true) & 
    (fun_med_results$bootstrap_results$beta_boot_basic_upper > beta_true);
  perc_cover <- (fun_med_results$bootstrap_results$beta_boot_perc_lower < beta_true) & 
    (fun_med_results$bootstrap_results$beta_boot_perc_upper > beta_true);
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
                     gamma_int_estimate_pvalue_overall =  as.numeric(tvem_model_summary$p.coeff["(Intercept)"]), 
                     gamma_int_estimate_pvalue_varying =  tvem_model_summary$s.table["s(time)","p-value"],
                     gamma_X_estimate_bias =  mean(gamma_X_residuals), 
                     gamma_X_estimate_mse =  mean(gamma_X_residuals^2),
                     gamma_X_estimate_mean_std_err = mean(gamma_X_se),
                     gamma_X_estimate_coverage =  mean(abs(gamma_X_residuals)<1.96*gamma_X_se), 
                     gamma_X_estimate_family_coverage =  all(abs(gamma_X_residuals)<1.96*gamma_X_se), 
                     gamma_X_estimate_pvalue_overall =  as.numeric(tvem_model_summary$p.coeff["(Intercept)"]),   
                     gamma_X_estimate_pvalue_varying =  tvem_model_summary$s.table["s(time):treatment","p-value"],  
                     # ... for the joint effect of X and M on Y: 
                     alpha_int_bias =  fun_med_results$original_results$alpha_int_estimate - alpha_int,
                     alpha_int_mse = (fun_med_results$original_results$alpha_int_estimate - alpha_int)^2,
                     alpha_int_se =  fun_med_results$original_results$alpha_int_se,
                     alpha_int_coverage = abs(fun_med_results$original_results$alpha_int_estimate - alpha_int) < 1.96*fun_med_results$original_results$alpha_int_se,
                     alpha_X_bias =  fun_med_results$original_results$alpha_X_estimate - alpha_X,
                     alpha_X_mse = (fun_med_results$original_results$alpha_X_estimate - alpha_X)^2,
                     alpha_X_se =  fun_med_results$original_results$alpha_X_se,
                     alpha_X_coverage = abs(fun_med_results$original_results$alpha_X_estimate - alpha_X) < 1.96*fun_med_results$original_results$alpha_X_se,
                     alpha_M_estimate_bias =  mean(alpha_M_residuals), 
                     alpha_M_estimate_mse =  mean(alpha_M_residuals^2), 
                     alpha_M_estimate_coverage =  mean(abs(alpha_M_residuals)<1.96*alpha_M_se), 
                     alpha_M_estimate_family_coverage =  all(abs(alpha_M_residuals)<1.96*alpha_M_se), 
                     alpha_M_pvalue =   fun_med_results$original_results$alpha_M_pvalue,  
                     # ... for the direct effect of X on Y: 
                     delta_int_estimate =  fun_med_results$original_results$delta_int_estimate, 
                     delta_int_se =  fun_med_results$original_results$delta_int_se, 
                     delta_X_estimate =  fun_med_results$original_results$delta_X_estimate, 
                     delta_X_se =  fun_med_results$original_results$delta_X_se,  
                     # ... for the mediation effect of X on Y through M:  
                     beta_bias = fun_med_results$bootstrap_results$beta_boot_estimate - beta_true,
                     beta_mse = (fun_med_results$bootstrap_results$beta_boot_estimate - beta_true)^2,
                     beta_boot_se = fun_med_results$bootstrap_results$beta_boot_se, 
                     beta_boot_norm_lower = fun_med_results$bootstrap_results$beta_boot_norm_lower, 
                     beta_boot_norm_upper = fun_med_results$bootstrap_results$beta_boot_norm_upper, 
                     beta_boot_norm_coverage = norm_cover,
                     beta_boot_basic_lower = fun_med_results$bootstrap_results$beta_boot_basic_lower, 
                     beta_boot_basic_upper = fun_med_results$bootstrap_results$beta_boot_basic_upper, 
                     beta_boot_basic_coverage = basic_cover, 
                     beta_boot_perc_lower = fun_med_results$bootstrap_results$beta_boot_perc_lower,
                     beta_boot_perc_upper = fun_med_results$bootstrap_results$beta_boot_perc_upper,
                     beta_boot_perc_coverage = perc_cover
                   ));
  save.image("working.rdata");
}}
print(c(nsim, nboot));
time_required <- difftime(current_time,start_time);
print(time_required);
mean_answers_100 <- apply(answers[which(answers[,"nsub"]==100),],2,mean);
mean_answers_300 <- apply(answers[which(answers[,"nsub"]==300),],2,mean);
print(round(cbind(mean_answers_100,
            mean_answers_300),3));



