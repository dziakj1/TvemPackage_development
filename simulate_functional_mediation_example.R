#' simulate_functional_mediation_example function
#' 
#' Simulates a dataset for demonstrating the functional_mediation function.
#' 
#' @param nsub Number of subjects
#' @param ntimes Number of potential times that could be observed on each subject
#' @param observe_rate Proportion of potential times on which there are actually 
#' observations. Not all times are observed; this assumed to be completely random and 
#' to be done by design to reduce participant burden.
#' @param gamma_int  Function representing the time-varying mean of mediator variable
#'  for the X=0 group
#' @param gamma_X Function representing the time-varying effect of X on the mediator
#' @param alpha_M Function representing the functional coefficient for cumulative 
#' (scalar-on-function) effect of M on Y adjusting for X
#' @param alpha_int  Mean of Y if the X is zero and M is the 0 function
#' @param alpha_X  Direct effect of X on Y after adjusting for M
#' @param sigma_Y Error standard deviation of the outcome Y 
#' @param sigma_M_error Error standard deviation of the mediator M
#' @param rho_M_error Autoregressive correlation coefficient of the mediator M from one observation to the next
#' @param simulate_binary_Y  Whether Y should be generated from a binary logistic (TRUE) or Gaussian (FALSE) model
#' 
#' @return A list with the following components:
#' \describe{
#' \item{time_grid}{The time grid for interpreting functional coefficients.}
#' \item{true_gamma_int}{True value of the time-varying gamma_int parameter, representing the time-specific mean of the mediator M when the treatment value X is 0.}
#' \item{true_gamma_X}{True value of the time-varying gamma_X parameter, representing the effect of X on M}
#' \item{true_alpha_int}{True value of the alpha_M parameter, representing the mean of the outcome Y when X=0 and M=0.}
#' \item{true_alpha_M}{True value of the alpha_M parameter, representing the functional effect of treatment on the outcome Y.}
#' \item{true_alpha_X}{True value of the alpha_X parameter, representing the effect of treatment on the outcome Y adjusting for the mediator.}
#' \item{true_beta}{True value of the beta parameter, representing the indirect (mediated) effect of treatment on the outcome Y.}
#' \item{dataset}{The simulated longitudinal dataset in long form.}
#' }   
#' 
#' @export

simulate_functional_mediation_example <- function(
  nsub = 250,
  ntimes = 100,
  observe_rate = .4, 
  gamma_int = function(t) {return(t^.5)}, # time-varying mean of mediator variable for the X=0 group;
  gamma_X = function(t) {return(-(t/2)^.5)}, # time-varying effect of X on the mediator;
  alpha_M = function(t) {(1/2)*(exp(t)-1)}, # functional (funreg) coefficient for cumulative effect of M on Y;
  alpha_int = 0,  # mean of Y if the X is zero and M is the 0 function; 
  alpha_X = .2,  # direct effect of X on Y after adjusting for M;
  sigma_Y = 1, 
  sigma_M_error = 2,
  rho_M_error = .8, 
  simulate_binary_Y=FALSE ) 
{
  time_grid <- (1:ntimes)/ntimes;  # vector of all possible times, scaled within 0 to 1;
  true_beta <- mean(alpha_M(time_grid)*gamma_X(time_grid));
  short_X <- rbinom(nsub,size=1,prob=.5); 
  # Simulate M from X...
  autoreg_error <- matrix(0,nsub,ntimes);
  autoreg_error[,1] <- rnorm(n=nsub,mean=0,sd=sigma_M_error);
  for (j in 2:ntimes) {
    autoreg_error[,j] <- rho_M_error*autoreg_error[,j-1] +
      sqrt(1-rho_M_error^2)*rnorm(n=nsub,mean=0,sd=sigma_M_error);
  }
  all_M <- matrix(0,nsub,ntimes);  # time-varying mediator; 
  for (i in 1:nsub) { 
    all_M[i,] <- gamma_int(time_grid) + short_X[i]*gamma_X(time_grid);
  } 
  all_M <- all_M + autoreg_error;
  if (simulate_binary_Y) {
    # Simulate Y from M and X...
    eta <- rep(NA,nsub); # = E(Y|X,M);
    for (i in 1:nsub) {
      eta[i] <- alpha_int + mean(alpha_M(time_grid) * all_M[i,]) + alpha_X*short_X[i];
    }
    mu <- exp(eta)/(1+exp(eta)); 
    short_Y <- unlist(lapply(X=mu,FUN=rbinom,size=1,n=1));
  } else {
    mu <- rep(NA,nsub); # = E(Y|X,M);
    for (i in 1:nsub) {
      mu[i] <- alpha_int + mean(alpha_M(time_grid) * all_M[i,]) + alpha_X*short_X[i];
    } 
    short_Y <- round(mu + rnorm(n=nsub,mean=0,sd=sigma_Y),5);
  }
  # Assemble simulated data into a long-form dataset:
  M <- all_M;
  for (i in 1:nsub) {
    which.missing.for.this.person <- which(rbinom(ntimes,1,1-observe_rate)==1);
    M[i,which.missing.for.this.person] <- NA;
  }
  temp <- data.frame(
    subject_id=rep(1:nsub,each=ntimes),
    t=rep(time_grid,times=nsub),
    X=rep(short_X,each=ntimes),
    M=as.vector(t(M)),
    Y=rep(short_Y,each=ntimes) 
  );
  long_simulated_data <- temp[which(!is.na(temp$M)),]; 
  return(list(time_grid=time_grid,
              true_gamma_int=gamma_int(time_grid),
              true_gamma_X=gamma_X(time_grid),
              true_alpha_int=alpha_int,
              true_alpha_M=alpha_M(time_grid),
              true_alpha_X=alpha_X,
              true_beta=true_beta,
              dataset=long_simulated_data));
}
