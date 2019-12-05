rm(list = ls());
library(refund);
library(boot); 
source("tvem.r");

set.seed(12022019);
nsub <- 300;          # number of subjects
ntimes <- 100;   # number of potential times that could be observed
nboot <- 199;  # number of bootstrapped samples to generate per original sample
observe.rate <- .2;   # proportion of potential times on which there are actually observations;
# by design to reduce participant burden.
#     Timeline:
time.grid <- 1:ntimes;  # vector of all possible times;
#     True coefficients:
gamma.int <- function(t) {return(.5*cos(3.14159*(t/max(t))))}; # time-varying mean of mediator variable for the X=0 group
gamma.X <- function(t) {return(sqrt(t/max(t)))}; # time-varying effect of X on the mediator
alpha.M <- function(t) {3-exp(t/max(t))}; # functional (funreg) coefficient for cumulative effect of M on Y
alpha.int <- 0;  # mean of Y if the X is zero and M is the 0 function; 
alpha.X <- .2;  # direct effect of X on Y after adjusting for M;
beta.true <- mean(alpha.M(time.grid)*gamma.X(time.grid));
#     Error standard deviations:
sigma.y <- 1; 
sigma.m.intercept <- 2;
sigma.m.slope <- 2;
sigma.m.error <- 2;
rho.m.error <- .8; 
# Simulate X...
save.image("starting.rdata");
X <- rbinom(n=nsub, size=1, prob=.5);  # binary treatment assignment;
# Simulate M from X...
AutoRegError <- matrix(0,nsub,ntimes);
AutoRegError[,1] <- rnorm(n=nsub,mean=0,sd=sigma.m.error);
for (j in 2:ntimes) {
  AutoRegError[,j] <- rho.m.error*AutoRegError[,j-1] +
    sqrt(1-rho.m.error^2)*rnorm(n=nsub,mean=0,sd=sigma.m.error);
}
all.M <- matrix(0,nsub,ntimes);  # time-varying mediator; 
for (i in 1:nsub) { 
  all.M[i,] <- gamma.int(time.grid) + X[i]*gamma.X(time.grid) +
    rnorm(n=1,mean=0,sd=sigma.m.intercept) + 
    rnorm(n=1,mean=0,sd=sigma.m.slope/max(time.grid))*time.grid;
} 
all.M <- all.M + AutoRegError;
# Simulate Y from M and X...
mu <- rep(NA,nsub); # = E(Y|X,M);
for (i in 1:nsub) {
  mu[i] <- alpha.int + mean(alpha.M(time.grid) * all.M[i,]) + alpha.X*X[i];
} 
Y <- mu + rnorm(n=nsub,mean=0,sd=sigma.y);
binY <- 1*(Y>mu);
# Assemble random data into a long-form dataset:
M <- all.M;
for (i in 1:nsub) {
  which.missing.for.this.person <- which(rbinom(ntimes,1,1-observe.rate)==1);
  M[i,which.missing.for.this.person] <- NA;
}
save.image("WorkingOnTest.rdata");
the.family <- "gaussian";
X[23] <- NA;
Y[32] <- NA;

windows();
source("MediationFunction.R");
#debug(MediationFunction);
answers1 <- MediationFunction(X=X,
                             M=M,
                             Y=Y,
                             binary=FALSE,
                             time=time.grid, 
                             nboot=99) ;
print(answers1$beta.estimate);
print(answers1$beta.boot.estimate);
print(answers1$beta.boot.se);
print(head(answers1$function.estimates));
windows();

answers2 <- MediationFunction(X=X,
                             M=M,
                             Y=binY,
                             binary=TRUE,
                             time=time.grid, 
                             nboot=99) 
print(answers2$beta.estimate);
print(answers2$beta.boot.estimate);
print(answers2$beta.boot.se);
print(head(answers2$function.estimates));

