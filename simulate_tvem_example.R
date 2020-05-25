#' simulate_tvem_example:  Simulate a dataset for demonstrating the tvem function. 
#' 
#' By default, the data-generating model has a time-varying intercept,
#' and two time-varying covariates named x1 and x2. 
#' x1 has a time-varying effect and x2 has no effect.
#' 
#' @param n_subjects Number of subjects in dataset 
#' @param max_time The time point at the end of the simulated time interval
#' 
#' @export

simulate_tvem_example <- function(
  n_subjects = 300, 
  max_time = 7,
  simulate_binary = FALSE) {
  n_obs_possible <- 141;
  prop_obs_observed <- .3; 
  possible_observation_times <- seq(0,max_time,length=n_obs_possible);
  scaled_times <- possible_observation_times/max_time;
  n_obs_per_subject <- rbinom(n_subjects,n_obs_possible,prop_obs_observed);
    # Generate X1;
  mu_x1 <- 6 - .5*exp(scaled_times);
  x1_short_term_rho <- .7;
  x1_error_term_AR1 <- matrix(NA, n_subjects, n_obs_possible);
  x1_error_term_AR1[,1] <- rnorm(n_subjects);
  for (j in 2:n_obs_possible) {
    x1_error_term_AR1[,j] <- x1_short_term_rho*x1_error_term_AR1[,j-1] + 
      sqrt(1-x1_short_term_rho^2)*rnorm(n_subjects);
  }
  x1 <- t(apply(x1_error_term_AR1,1,"+",mu_x1));
  x1 <- round(x1,3);
  # Generate X2;
  mu_x2 <- 3 + scaled_times;
  x2_short_term_rho <- .7;
  x2_error_term_AR1 <- matrix(NA, n_subjects, n_obs_possible);
  x2_error_term_AR1[,1] <- rnorm(n_subjects);
  for (j in 2:n_obs_possible) {
    x2_error_term_AR1[,j] <- x2_short_term_rho*x2_error_term_AR1[,j-1] + 
      sqrt(1-x2_short_term_rho^2)*rnorm(n_subjects);
  }
  x2 <- t(apply(x2_error_term_AR1,1,"+",mu_x2));
  x2 <- round(x2,3);
  # Generate Y;
  sigma_y <- 2;
  #beta0_y <- 5 - exp(scaled_times);
  beta0_y <- -3+exp(scaled_times);
  beta1_y <- .1*exp(scaled_times);
  beta2_y <- rep(0,length(scaled_times));
  if (simulate_binary) {
    eta_y <- beta0_y + beta1_y*x1 + beta2_y*x2;
    mu_y <- plogis(eta_y);
    stopifnot((min(mu_y)<.75 | max(mu_y)>.25));
    y <- apply(mu_y,MARGIN=c(1,2),FUN=rbinom,n=1,size=1);
  } else {
    mu_y <- beta0_y + beta1_y*x1 + beta2_y*x2;
    y_short_term_rho <- .5;
    y_error_term_AR1 <- matrix(NA, n_subjects, n_obs_possible);
    y_error_term_AR1[,1] <- rnorm(n_subjects);
    for (j in 2:n_obs_possible) {
      y_error_term_AR1[,j] <- y_short_term_rho*y_error_term_AR1[,j-1] + 
        sqrt(1-y_short_term_rho^2)*rnorm(n_subjects);
    }
    y <- mu_y + sigma_y * y_error_term_AR1;
    y <- round(y,3);
  } 
  for (this_subject in 1:n_subjects) {
    unobserved_for_this_subject <- sample(1:ncol(y),size=n_obs_possible - n_obs_per_subject[this_subject]);
    y[this_subject,unobserved_for_this_subject] <- NA;
  }
  entire_long_dataset <- data.frame(subject_id=rep(1:n_subjects,each=n_obs_possible),
                                    time=rep(possible_observation_times,times=n_subjects),
                                    x1=as.vector(t(x1)),
                                    x2=as.vector(t(x2)),
                                    y=as.vector(t(y)));
  long_dataset <- entire_long_dataset[which(is.na(entire_long_dataset$y)==FALSE),]
  return(long_dataset);
  return(sim_data);
}