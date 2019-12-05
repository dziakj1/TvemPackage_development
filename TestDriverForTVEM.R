 
rm(list = ls());
source("TVEM.r")
set.seed(123);
nsub <- 500;
nobs <- 50;
min.time <- 0;
max.time <- 7;
fit.df <- 20;
subject_id <- rep(1:nsub,each=nobs);
times <- as.vector(replicate(expr=sort(runif(nobs,min.time,max.time)),
                             n=nsub)); 
true.b0 <- 5+ sin(2*pi*((times-min.time)/max.time));   
true.b1 <-  3-cos(2*pi*((times-min.time)/max.time));
true.b2 <- -8;
tvx <- rnorm(nsub*nobs,mean=5);
# Note that if this distribution does not overlap with zero then
# the b0(t) function may have a very wide confidence interval.
# Maybe offer some kind of warning if min(x)>0 or max(x)<0.
# But perhaps it is okay for min(x) to be a little above zero,  
# as in Likert scale, but just not too far (like abs(mean(x)/sd(x))>6
# maybe? or min(x)-sd(x)>0 or max(x)+sd(x)<0?
t.by.tvx <- times*tvx;
z1 <- rnorm(nsub*nobs,mean=2);
z2 <- runif(nsub*nobs);
mu <- true.b0 + true.b1*tvx + true.b2*z1;
error.obs <- rep(rnorm(nsub),nobs);
error.sub <- rnorm(nsub*nobs); 
outcome <- round(mu + error.sub + error.obs, 2);
sim.data <- data.frame(subject_id=subject_id,
                       times=times,
                       tvx=tvx,
                       z1=z1,
                       z2=z2,
                       y=outcome);
ans <- tvem(data=sim.data,
            formula = y~tvx,
            id=subject_id,
            time=times,
            invar_effects = ~z1+z2);
b0 <- ans$fitted_coefficients[[1]];
b1 <- ans$fitted_coefficients[[2]];
grid_b0 <- ans$grid_fitted_coefficients[[1]];
grid_b1 <- ans$grid_fitted_coefficients[[2]];


##################################
# Draw plots and check accuracy versus data-generating
# model of the simulation;
##################################
par(mfrow=c(2,2));
plot(times,b0$estimate,col="red",ylab="est b0",ylim=c(min(true.b0)-1,max(true.b0)+1)); 
points(times,true.b0,col="blue",cex=.5);
points(times,b0$upper,col="green",cex=.5);
points(times,b0$lower,col="green",cex=.5);
plot(times,b1$estimate,col="red",ylab="est b1",ylim=c(min(true.b1)-1,max(true.b1)+1)); 
points(times,true.b1,col="blue",cex=.5);
points(times,b1$upper,col="green",cex=.5);
points(times,b1$lower,col="green",cex=.5);
plot(ans$time_grid,grid_b0$estimate,col="red",ylab="est b0",ylim=c(min(true.b0)-1,max(true.b0)+1)); 
points(times,true.b0,col="blue",cex=.5);
points(ans$time_grid,grid_b0$upper,col="green",cex=.5);
points(ans$time_grid,grid_b0$lower,col="green",cex=.5);
plot(ans$time_grid,grid_b1$estimate,col="red",ylab="est b1",ylim=c(min(true.b1)-1,max(true.b1)+1)); 
points(times,true.b1,col="blue",cex=.5);
points(ans$time_grid,grid_b1$upper,col="green",cex=.5);
points(ans$time_grid,grid_b1$lower,col="green",cex=.5);

