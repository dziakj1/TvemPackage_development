#' simulate_tvem_example:  Simulate a dataset for demonstrating the tvem function. 
#' 
#' By default, the data-generating model has a time-varying intercept,
#' and a single time-varying-effects covariate named tvx.  There are two
#' covariates named z1 and z2; the first has a time-invariant effect
#' and the second has no effect.
#' 
#' @param nsub Number of subjects in dataset
#' @param nobs Number of observations per subject
#' @param min_time The time point at the beginning of the simulated time interval
#' @param max_time The time point at the end of the simulated time interval
#' 
#' @export

simulate_tvem_example <- function(
  nsub = 500,
  nobs = 50,
  min_time = 0,
  max_time = 7) {
  subject_id <- rep(1:nsub,each=nobs);
  times <- as.vector(replicate(expr=sort(runif(nobs,min_time,max_time)),
                               n=nsub)); 
  true_b0 <- 5+ sin(2*pi*((times-min_time)/max_time));   
  true_b1 <-  3-cos(2*pi*((times-min_time)/max_time));
  true_b2 <- -8;
  tvx <- rnorm(nsub*nobs,mean=5); 
  t_by_tvx <- times*tvx;
  z1 <- rnorm(nsub*nobs,mean=2);
  z2 <- runif(nsub*nobs);
  mu <- true_b0 + true_b1*tvx + true_b2*z1;
  error_obs <- rep(rnorm(nsub),nobs);
  error_sub <- rnorm(nsub*nobs); 
  outcome <- round(mu + error_sub + error_obs, 2);
  sim_data <- data.frame(subject_id=subject_id,
                         times=times,
                         tvx=tvx,
                         z1=z1,
                         z2=z2,
                         y=outcome);
  return(sim_data)
}