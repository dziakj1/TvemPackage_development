rm(list = ls());
library(refund);
library(boot); 
source("tvem.r");
source("plot_tvem.r");
source("print_tvem.r");
set.seed(12022019);
nsub <- 200;          # number of subjects
ntimes <- 100;   # number of potential times that could be observed
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
Z1 <- rpois(nsub,lambda=1);
Z2 <- rpois(nsub,lambda=2);
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

temp <- data.frame(
  subject_id=rep(1:nsub,each=ntimes),
  time_value=rep(time.grid,times=nsub),
  X=rep(X,each=ntimes),
  M=as.vector(t(M)),
  Y=rep(Y,each=ntimes),
  Z1=rep(Z1,each=ntimes),
  Z2=rep(Z2,each=ntimes)
);
long_simulated_data <- temp[which(!is.na(temp$M)),];

save.image("WorkingOnTest.rdata");
source("NewFunregMediation.R");

source("TVEM.r");
the_tvem <- tvem(data=long_simulated_data,
                 formula=M~X+Z1,
                 invar_effects=~Z2,
                 id=subject_id,
                 time=time_value);
print_tvem(the_tvem); 
plot_tvem(the_tvem);
plot_tvem(the_tvem,diagnostics=TRUE);