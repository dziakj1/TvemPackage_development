#' simulate_functional_mediation_example function
#' 
#' Simulates a dataset for demonstrating the functional_mediation function.
#' 
#' @param nsub Number of subjects
#' @param ntimes Number of potential times that could be observed on each subject
#' @param observe_rate proportion of potential times on which there are actually observations
#' @param gamma_int  Function representing the time-varying mean of mediator variable for the X=0 group
#' @param gamma_X Function representing the time-varying effect of X on the mediator
#' @param alpha_M Function representing the functional coefficient for cumulative (scalar-on-function) effect of M on Y adjusting for X
#' @param alpha_int  Mean of Y if the X is zero and M is the 0 function
#' @param alpha_X  Direct effect of X on Y after adjusting for M
#' @param sigma_Y Error standard deviation of the outcome Y 
#' 
#' @return A simulated longitudinal dataset in long form.
#' @export

simulate_functional_mediation_example <- function(
  nsub = 200,
  ntimes = 100,
  observe_rate = .2,
  gamma_int = function(t) {return(.5*cos(3.14159*(t/max(t))))},
  gamma_X = function(t) {return(sqrt(t/max(t)))},
  alpha_M = function(t) {3-exp(t/max(t))},
  alpha_int = 0, 
  alpha_X = .2,
  sigma_Y = 1,
  sigma_M_intercept = 2,
  sigma_M_slope = 2,
  sigma_M_error = 2,
  rho_M_error = .8 )
  {
  # by design to reduce participant burden.
  #     Timeline:
  time_grid <- 1:ntimes;  # vector of all possible times;
  beta.true <- mean(alpha_M(time_grid)*gamma_X(time_grid));
  #     Error standard deviations:
  # Simulate X...
  save.image("starting.rdata");
  X <- rbinom(n=nsub, size=1, prob=.5);  # binary treatment assignment;
  Z1 <- rpois(nsub,lambda=1);
  Z2 <- rpois(nsub,lambda=2);
  # Simulate M from X...
  AutoRegError <- matrix(0,nsub,ntimes);
  AutoRegError[,1] <- rnorm(n=nsub,mean=0,sd=sigma_M_error);
  for (j in 2:ntimes) {
    AutoRegError[,j] <- rho_M_error*AutoRegError[,j-1] +
      sqrt(1-rho_M_error^2)*rnorm(n=nsub,mean=0,sd=sigma_M_error);
  }
  all.M <- matrix(0,nsub,ntimes);  # time-varying mediator; 
  for (i in 1:nsub) { 
    all.M[i,] <- gamma_int(time_grid) + X[i]*gamma_X(time_grid) +
      rnorm(n=1,mean=0,sd=sigma_M_intercept) + 
      rnorm(n=1,mean=0,sd=sigma_M_slope/max(time_grid))*time_grid;
  } 
  all.M <- all.M + AutoRegError;
  # Simulate Y from M and X...
  mu <- rep(NA,nsub); # = E(Y|X,M);
  for (i in 1:nsub) {
    mu[i] <- alpha_int + mean(alpha_M(time_grid) * all.M[i,]) + alpha_X*X[i];
  } 
  Y <- mu + rnorm(n=nsub,mean=0,sd=sigma_Y);
  binY <- 1*(Y>mu);
  # Assemble random data into a long-form dataset:
  M <- all.M;
  for (i in 1:nsub) {
    which.missing.for.this.person <- which(rbinom(ntimes,1,1-observe_rate)==1);
    M[i,which.missing.for.this.person] <- NA;
  }
  temp <- data.frame(
    subject_id=rep(1:nsub,each=ntimes),
    time_value=rep(time_grid,times=nsub),
    X=rep(X,each=ntimes),
    M=as.vector(t(M)),
    Y=rep(Y,each=ntimes),
    Z1=rep(Z1,each=ntimes),
    Z2=rep(Z2,each=ntimes)
  );
  long_simulated_data <- temp[which(!is.na(temp$M)),]; 
  return(long_simulated_data);
}

