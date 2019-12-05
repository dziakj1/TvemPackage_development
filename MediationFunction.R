MediationFunction <- function(X,             # Vector of 0's and 1's giving treatment group membership for each subject
                              M,             # Matrix of observed variables of the mediator.  Values can be NA.  Rows correspond to subjects, and columns to observation times.
                              Y,             # Vector of numbers (for binary==FALSE) or 0's and 1's (for binary==TRUE) giving observed responses for each subject
                              time,     # Vector giving the value on the time variable for each column of M
                              binary,        # TRUE for a binary logistic model and FALSE for a normal linear model 
                              nboot,         # number of bootstrap samples to use for test;
                              plot=TRUE   # whether or not to draw plots;
) {
  time.grid <- time;
  do.plot <- plot;
  wide.data <- data.frame(Y=Y,X=X,M=M);
  # How are we handling time grid here?
  # TWO HELPER FUNCTIONS:
  # Function for analyzing data:
  analyze.data.outside.bootstrap <- function(the.wide.data) {
    local.Y <- the.wide.data[,"Y"];
    local.X <- the.wide.data[,"X"];
    local.M <- the.wide.data[,3:ncol(the.wide.data)];
    local.id <- row(local.M);
    long.data <- data.frame(M=as.vector(t(local.M)),
                            id=as.vector(t(local.id)),
                            time=rep(time.grid,times=nrow(local.M)),
                            X=rep(local.X,each=ncol(local.M)));
    for (i in 1:nsub) {
      long.data$X[which(long.data$id==i)] <- local.X[i];
    }
    long.usable.data <- long.data[which(!is.na(long.data$M)),]; 
    #--- EFFECT OF TREATMENT X ON MEDIATOR M ---; 
    tvem1 <- tvem(data=long.usable.data,
                  formula=M~X,
                  time=time,
                  id=id,
                  grid=time.grid);
    gamma.int.estimate <- tvem1$grid_fitted_coefficients[[1]]$estimate; 
    gamma.int.se <- tvem1$grid_fitted_coefficients[[1]]$standard_error;
    gamma.X.estimate <- tvem1$grid_fitted_coefficients[[2]]$estimate; 
    gamma.X.se <- tvem1$grid_fitted_coefficients[[2]]$standard_error;
    #--- EFFECT OF MEDIATOR M AND TREATMENT X ON OUTCOME Y ---;
    local.Y <- unlist(the.wide.data["Y"]);
    local.X <- unlist(the.wide.data["X"]);
    local.M <- as.matrix(the.wide.data[,3:ncol(the.wide.data)]);
    if (binary) {
      funreg.MY <- pfr(local.Y~lf(local.M,
                                  presmooth="interpolate")+local.X,
                       family="binomial");  # cannot use data= inside a function in pfr;
    } else {
      funreg.MY <- pfr(local.Y~lf(local.M,
                                  presmooth="interpolate")+local.X,
                       family="gaussian");  # cannot use data= inside a function in pfr;
    }
    alpha.int.estimate <- as.numeric(funreg.MY$coefficient["(Intercept)"]);
    alpha.int.se <- as.numeric(summary(funreg.MY)$se["(Intercept)"]);
    alpha.X.estimate <-  as.numeric(funreg.MY$coefficient["local.X"]);
    alpha.X.se <- as.numeric(summary(funreg.MY)$se["local.X"]);
    alpha.M.estimate <- as.numeric(coef(funreg.MY,n=ntimes)[,"value"]); 
    alpha.M.se <- coef(funreg.MY,n=ntimes)[,"se"]; 
    alpha.M.true <- alpha.M(time.grid);
    alpha.M.pvalue <- summary(funreg.MY)$s.table[1,"p-value"];
    #--- DIRECT EFFECT OF TREATMENT X ON OUTCOME Y ---;
    if (binary) {
      delta.effect.model <- glm(Y~X,data=the.wide.data,family=binomial);
    } else {
      delta.effect.model <- glm(Y~X,data=the.wide.data);
    }
    delta.int.estimate <- as.numeric(delta.effect.model$coefficients["(Intercept)"]);
    delta.int.se <- summary(delta.effect.model)$coefficients["(Intercept)","Std. Error"];
    delta.X.estimate <- as.numeric(delta.effect.model$coefficients["X"]);
    delta.X.se <- summary(delta.effect.model)$coefficients["X","Std. Error"];
    #--- MEDIATED EFFECT OF TREATMENT X THROUGH MEDIATOR M ON OUTCOME Y ---;
    beta.estimate <- mean(gamma.X.estimate*alpha.M.estimate);
    function.estimates <- cbind(time=time.grid,
                                gamma.int.est=round(gamma.int.estimate,4),
                                gamma.int.se=round(gamma.int.se,4),
                                gamma.X.est=round(gamma.X.estimate,4),
                                gamma.X.se=round(gamma.X.se,4),
                                alpha.M.est=round(alpha.M.estimate,4),
                                alpha.M.se=round(alpha.M.se,4));
    results <- list ( 
      alpha.X.estimate = alpha.X.estimate,
      alpha.X.se = alpha.X.se,
      delta.int.estimate = delta.int.estimate,
      delta.int.se = delta.int.se,
      delta.X.estimate = delta.X.estimate,
      delta.X.se = delta.X.se,
      beta.estimate = beta.estimate,
      function.estimates = function.estimates
    ); 
    return(results);
  }
  
  analyze.data.inside.bootstrap <- function(the.wide.data, indices) { 
    indices <- indices; 
    bootstrapped.data <- wide.data[indices,];
    local.Y <- bootstrapped.data[,"Y"];
    local.X <- bootstrapped.data[,"X"];
    local.M <- bootstrapped.data[,3:ncol(the.wide.data)];
    local.id <- row(local.M);
    long.data <- data.frame(M=as.vector(t(local.M)),
                            id=as.vector(t(local.id)),
                            time=rep(time.grid,times=nrow(local.M)),
                            X=rep(local.X,each=ncol(local.M)));
    for (i in 1:nsub) {
      long.data$X[which(long.data$id==i)] <- local.X[i];
    }
    long.usable.data <- long.data[which(!is.na(long.data$M)),]; 
    #--- EFFECT OF TREATMENT X ON MEDIATOR M ---; 
    tvem1 <- tvem(data=long.usable.data,
                  formula=M~X,
                  time=time,
                  id=id,
                  grid=time.grid);
    gamma.int.estimate <- tvem1$grid_fitted_coefficients[[1]]$estimate; 
    gamma.int.se <- tvem1$grid_fitted_coefficients[[1]]$standard_error;
    gamma.X.estimate <- tvem1$grid_fitted_coefficients[[2]]$estimate; 
    gamma.X.se <- tvem1$grid_fitted_coefficients[[2]]$standard_error;
    #--- EFFECT OF MEDIATOR M AND TREATMENT X ON OUTCOME Y ---;
    local.Y <- unlist(the.wide.data["Y"]);
    local.X <- unlist(the.wide.data["X"]);
    local.M <- as.matrix(the.wide.data[,3:ncol(the.wide.data)]);
    if (binary) {
      funreg.MY <- pfr(local.Y~lf(local.M,
                                  presmooth="interpolate")+local.X,
                       family="binomial");  # cannot use data= inside a function in pfr;
    } else {
      funreg.MY <- pfr(local.Y~lf(local.M,
                                  presmooth="interpolate")+local.X,
                       family="gaussian");  # cannot use data= inside a function in pfr; 
    }
    alpha.int.estimate <- as.numeric(funreg.MY$coefficient["(Intercept)"]);
    alpha.X.estimate <-  as.numeric(funreg.MY$coefficient["local.X"]);
    alpha.M.estimate <- coef(funreg.MY,n=ntimes)[,"value"];
    alpha.M.true <- alpha.M(time.grid);
    alpha.M.pvalue <- summary(funreg.MY)$s.table[1,"p-value"];
    #--- MEDIATED EFFECT OF TREATMENT X THROUGH MEDIATOR M ON OUTCOME Y ---;
    beta.estimate <- mean(gamma.X.estimate*alpha.M.estimate);
    return(beta.estimate);
  }
  # MAIN BODY OF FUNCTION:
  original.results <- analyze.data.outside.bootstrap(the.wide.data=wide.data); 
  # Bootstrap:
  before.boot <- Sys.time();
  boot1 <- boot(data=wide.data,statistic=analyze.data.inside.bootstrap,R=nboot);
  boot2 <- boot.ci(boot1,conf=.95,type="norm");
  boot3 <- boot.ci(boot1,conf=.95,type="basic");
  boot4 <- boot.ci(boot1,conf=.95,type="perc");
  after.boot <- Sys.time();
  boot.results <- list(beta.boot.estimate= norm.ci(boot1,conf=.001)[2],
                       beta.boot.se=sd(boot1$t),
                       beta.boot.norm.lower=boot2$normal[2],
                       beta.boot.norm.upper=boot2$normal[3],
                       beta.boot.basic.lower=boot3$basic[4],
                       beta.boot.basic.upper=boot3$basic[5],
                       beta.boot.perc.lower=boot4$percent[4],
                       beta.boot.perc.upper=boot4$percent[5],
                       boot1=boot1,
                       time.required=difftime(after.boot,before.boot));
  results <- c(original.results,
               boot.results);
  if (do.plot) {
    par(mfrow=c(2,2),mar=c(4,4,2,2),mgp=c(2,.5,0));
    hist(boot1$t,main="Histogram of bootstrap estimates",cex.lab=1.3,xlab=expression(beta[X]));
    abline(v=boot1$t0,lwd=2)
    #plot(time.grid,
    #     M[1,],
    #     ylim=c(min(M,na.rm=TRUE),max(M,na.rm=TRUE)),
    #     type="p",pch=16,
    #     main="Observed mediator distribution",
    #     xlab="t",ylab="M(t)");
    #for (i in 2:nrow(M)) {
    #  points(time.grid,M[i,],pch=16);
    #}
    plot(original.results$function.estimates[,"time"],
         original.results$function.estimates[,"gamma.int.est"],
         main=expression("Intercept gamma.int of mediator"),
         type="l",lwd=2,cex.lab=1.3,
         ylim=c(pmin(0,min(original.results$function.estimates[,"gamma.int.est"]-1.96*original.results$function.estimates[,"gamma.int.se"])),
                pmax(0,max(original.results$function.estimates[,"gamma.int.est"]+1.96*original.results$function.estimates[,"gamma.int.se"]))),
         xlab=expression(t),
         ylab=expression(alpha[int](t)));
    lines(original.results$function.estimates[,"time"],
          original.results$function.estimates[,"gamma.int.est"]-1.96*original.results$function.estimates[,"gamma.int.se"],lty="dashed");
    lines(original.results$function.estimates[,"time"],
          original.results$function.estimates[,"gamma.int.est"]+1.96*original.results$function.estimates[,"gamma.int.se"],lty="dashed");
    plot(original.results$function.estimates[,"time"],
         original.results$function.estimates[,"gamma.X.est"],
         ylim=c(pmin(0,min(original.results$function.estimates[,"gamma.X.est"]-1.96*original.results$function.estimates[,"gamma.X.se"])),
                pmax(0,max(original.results$function.estimates[,"gamma.X.est"]+1.96*original.results$function.estimates[,"gamma.X.se"]))),
         main=expression("Effect gamma.X of treatment on mediator"),
         type="l",lwd=2,cex.lab=1.3,cex.main=.9,
         xlab=expression(t),
         ylab=expression(gamma[X](t)));
    lines(original.results$function.estimates[,"time"],
          original.results$function.estimates[,"gamma.X.est"]-1.96*original.results$function.estimates[,"gamma.X.se"],lty="dashed");
    lines(original.results$function.estimates[,"time"],
          original.results$function.estimates[,"gamma.X.est"]+1.96*original.results$function.estimates[,"gamma.X.se"],lty="dashed");
    plot(original.results$function.estimates[,"time"],
         original.results$function.estimates[,"alpha.M.est"],
         main=expression("Effect alpha.M of mediator on response"),
         type="l",lwd=2,cex.lab=1.3,
         ylim=c(pmin(0,min(original.results$function.estimates[,"alpha.M.est"]-1.96*original.results$function.estimates[,"alpha.M.se"])),
                pmax(0,max(original.results$function.estimates[,"alpha.M.est"]+1.96*original.results$function.estimates[,"alpha.M.se"]))),
         xlab=expression(t),
         ylab=expression(alpha[M](t)));
    lines(original.results$function.estimates[,"time"],
          original.results$function.estimates[,"alpha.M.est"]-1.96*original.results$function.estimates[,"alpha.M.se"],lty="dashed");
    lines(original.results$function.estimates[,"time"],
          original.results$function.estimates[,"alpha.M.est"]+1.96*original.results$function.estimates[,"alpha.M.se"],lty="dashed");
  }
  return(results);
}
