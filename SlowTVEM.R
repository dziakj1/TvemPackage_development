library(splines); 

SlowTVEM <- function(dep,
                 id,
                 tcov,
                 time,
                 dist="normal", 
                 convergenceCriterion=1e-6, 
                 deg=3,
                 doPlot=TRUE,
                 gridSize=100,
                 maxIterations=5000, # Maximum number of iterations to attempt;
                 min.time=NA,
                 max.time=NA,
                 numInteriorKnots=10,
                 useAIC=FALSE,
                 xcov=NULL) {
  ##################################################################
  # SlowTVEM R macro Version 0.1
  # By John J. DZIAK and Runze LI
  # Based on MixTvem macro 1.1 by John J. Dziak, Xianming Tan and Runze Li
  # Fits a time-varying effects model.
  ##################################################################
  # Copyright:
  # (c) 2019 The Pennsylvania State University
  ##################################################################
  # License:
  # This program is free software; you can redistribute it and/or
  # modify it under the terms of the GNU General Public License as
  # published by the Free Software Foundation; either version 2 of
  # the License, or (at your option) any later version.
  # This program is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  # General Public License for more details.
  ##################################################################
  # Acknowledgments and references:
  # We fit a nonparametric varying-coefficient model, using
  # a penalized B-spline approach.  See
  #  Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing using
  #    B-splines and penalized likelihood. Statistical Science
  #    11(2): 89-121.
  #  Hastie, T., & Tibshirani, R. (1993). Varying-coefficient models.
  #    Journal of the Royal Statistical Society, Series B, 55, 757-796.
  #  Tan, X., Shiyko, M. P., Li, R., Li, Y., & Dierker, L. (2011, November 21).
  #    A time-varying effect model for intensive longitudinal data.
  #    Psychological Methods. Advance online publication.
  #    doi: 10.1037/a0025814.
  # The use of a sandwich formula to adjust for within-subject correlation
  # when calculating the standard errors is inspired by
  #  Liang, K.-Y., and Zeger, S. L. (1986). Longitudinal data analysis
  #    using generalized linear models. Biometrika, 73, 13-22.
  # The model fit criteria used are adapted versions of the standard AIC and BIC.
  # Akaike, H. (1973). Information theory and an extension of the maximum
  #     likelihood principle. In B. N. Petrov & F. Csaki (Eds.), Second
  #     international symposium on information theory (p. 267-281).
  #     Budapest, Hungary: Akademai Kiado.
  #   Schwarz, G. (1978). Estimating the dimension of a model. Annals of
  #     Statistics, 6, 461-464. 
  ##################################################################
  ## Define required helper functions:
  PenGEE <- function(   dist,
                        intId,
                        time,
                        X,
                        Y,
                        convergenceCriterion=1e-8,
                        maxIterations=1000,
                        penaltyMatrix=NULL,
                        roughnessPenalty=0
  ) {
    # Fit penalized working-independence GEE;
    # See Dziak & Li (2006) tech report.
    # The PGEE R library by Inan and Wang does something like this
    # but we are using the difference-based penalty of Eilers and Marx (1996)
    # instead of the SCAD penalty.  Also, for now we assume working independence;
    ## Prepare to begin loop;
    stopifnot(length(time)==length(intId));
    theta <- rep(0,ncol(X));
    thetaOld <- Inf;
    maxAbsDev <- Inf;
    numSubjects <- max(intId);
    numTotal <- length(intId);
    iteration <- 1;
    numObsBySub <- table(intId); # assumes that intId
    # consists of consecutive integers starting at 1;
    stopifnot(length(Y)==numTotal);
    stopifnot(length(unique(intId))==numSubjects);
    stopifnot(all.equal(unique(intId),1:numSubjects));
    stopifnot(as.integer(rownames(numObsBySub))==1:numSubjects);
    stopifnot(identical(intId,sort(intId)));
    stopifnot(length(numObsBySub)==numSubjects);
    stopifnot(length(intId)==numTotal);
    stopifnot(identical(intId,sort(intId)));
    dist <- tolower(dist);
    if (dist=="gaussian") {dist <- "normal";}
    if (dist=="binary" | dist=="binomial") {dist <- "logistic";}
    stopifnot((dist=="normal")|(dist=="logistic")|(dist=="poisson"));
    sigsq = 1;
    ## Begin loop;
    iteration <- 0;
    maxAbsDev <- Inf;
    thetaOld <- rep(0,ncol(X));
    while ((iteration<maxIterations)&(maxAbsDev>convergenceCriterion)) {
      ## Initial work;
      iteration <- iteration+1;
      hessian <- matrix(0,ncol(X),ncol(X));
      score <- matrix(0,ncol(X),1);
      for (i in 1:numSubjects) {
        these <- which(intId==i);
        Yi <- Y[these];
        Xi <- X[these,,drop=FALSE];
        etaOld <- as.vector(Xi %*% thetaOld);
        if (dist=="normal") {
          mui <- etaOld;
          derivs <- 1;
        }
        if (dist=="logistic") {
          mui <- exp(etaOld)/(1+exp(etaOld));
          derivs <- mui * (1-mui); 
        }
        if (dist=="poisson") {
          mui <- exp(etaOld);
          derivs <- mui;
        }
        hessian <- hessian +crossprod(Xi*derivs,Xi);
        score<- score + crossprod(Xi,Yi-mui);
      }
      theta <- thetaOld + solve(hessian + 
                                  roughnessPenalty*penaltyMatrix,  
                                score- roughnessPenalty*penaltyMatrix%*%thetaOld);
      etaHat <- X%*%theta;
      if (dist=="normal") {
        fittedY <- etaHat;
      }
      if (dist=="logistic") {
        etaHat[which(etaHat < -12)] <- -12;
        etaHat[which(etaHat > +12)] <- +12;
        fittedY <- exp(etaHat)/(1+exp(etaHat));
      }
      if (dist=="poisson") {
        fittedY <- exp(etaHat);
      }
      stopifnot(length(Y)==nrow(fittedY));
      maxAbsDev <- max(abs(theta-thetaOld)); 
      thetaOld <- as.vector(theta);
    }
    converged <- maxAbsDev<=convergenceCriterion; 
    ## Done with main loop
    # Calculate sigma;
    residsY <- Y-fittedY;
    if (dist=="normal") {
      ## Get new sigma estimates;
      f <- function(i) {
        return(sum(residsY[which(intId==i)]^2))
      };
      rssBySubject <- as.vector(sapply(1:numSubjects,f));
      sigsq <- sum(rssBySubject) / sum(numObsBySub);
    } else {
      sigsq <- NULL;
    }
    # Calculate effective number of parameters;
    enp <- 1 + sum(diag(solve(hessian+roughnessPenalty*penaltyMatrix, hessian)));
    np <- 1 + ncol(X) ; # the real total number of parameters; 
    # Calculate log-likelihood;
    if (dist=="normal") {
      logLikelihoodBySubject <- -(numObsBySub/2)*
        log(2*3.1415926535*(sigsq+1e-30)) -
        rssBySubject/(2*sigsq+1e-30);
      logLik <- sum(logLikelihoodBySubject);
    }
    if (dist=="logistic") {
      logLik <- sum(Y * log(fittedY)) + sum((1-Y) * log(1-fittedY));
    }
    if (dist=="poisson") {
      logLik <- sum(Y * log(fittedY) - fittedY);
      # Note:  this should also include a term +factorial(Y) but 
      # that is computationally awkward and is dependent only on 
      # the data, not on the model.	
    }
    return(list(aic=-2*logLik+2*enp,
                # A statistic similar to AIC (Akaike, 1973)
                # but using the effective (penalized) number of
                # parameters instead of a count of parameters
                # (see Eilers and Marx, 1996)
                bic=-2*logLik+log(numSubjects)*enp,
                # A statistic somewhat analogous to BIC
                # (Schwarz, 1978).  The effective (penalized)
                # number of parameters is used instead of a
                # count of parameters (see Eilers and Marx, 1996).
                converged=converged,
                # Indicates whether the EM algorithm converged
                dep=Y,
                # The observed response for each assessment,
                # copied here for convenience when doing
                # fit diagnostics, follow-up analyses, etc.
                enp=enp,
                # The effective number of regression
                # parameters (see Eilers and Marx, 1996)
                fittedY=fittedY,
                # The fitted values for the dependent
                # variable at each assessment.
                intId=intId,
                # The internally assigned subject identification
                # number for each assessment,
                # copied here for convenience.
                iteration=iteration,
                # The number of iterations run by the EM
                # algorithm.
                lambda=roughnessPenalty,
                # The current candidate value of the penalty
                # weight for the roughness penalty
                logLik=logLik,
                # The fitted log-likelihood
                np=np,
                # The counted number of parameters,
                # ignoring shrinkage imposed by the penalty
                residsY = residsY,
                # Residuals for each observation;
                sigsq=sigsq,
                # Estimated sigma squared (variance);
                theta=theta
                # Estimated regression coefficients for each
                # of the regression coefficients, including
                # spline basis coefficients 
    ));
  }
  PenGEECovMats <- function(dist,
                            intId,
                            mu,
                            numSubjects,
                            numTotal,
                            penalty,
                            sigsq=1,  # leave at default of 1 if not Gaussian;
                            theta,
                            X,
                            Y) {
    numthetas <- nrow(theta);
    stopifnot(numthetas==ncol(X));
    numParams <- numthetas;
    I1 <- matrix(0,numParams,numParams);
    I3 <- matrix(0,numParams,numParams);
    for (i in 1:numSubjects) {
      these <- which(intId==i);
      Xi <- as.matrix(X[which(intId==i),,drop=FALSE]);
      if (dist=="normal") {
        derivs <- 1/sigsq;
      }
      if (dist=="logistic") {
        etai <- as.vector(Xi %*%theta);
        mui <- exp(etai)/(1+exp(etai));
        derivs <- mui*(1-mui);
      }
      if (dist=="poisson") {
        etai <- as.vector(Xi %*%theta);
        mui <- exp(etai);
        derivs <- mui;
      } 
      ni <- sum(intId==i);
      contrib <- crossprod(Xi*derivs,Xi);
      I1 <- I1 + contrib;
    }
    #print("Check adding penalty"); 
    PenaltyRidge <- (1e-8)*mean(diag(penalty));
    diag(penalty) <- diag(penalty) + PenaltyRidge;   
    I1 <- I1 + penalty; 
    # Calculate matrix I3
    for (i in 1:numSubjects) {
      scorei <- matrix(0,numParams,1);
      ni <- sum(intId==i);
      xi <- X[which(intId==i),,drop=FALSE];
      yi <- Y[which(intId==i),drop=FALSE]; 
      mui <- mu[which(intId==i),drop=FALSE];
      hik <- (1/sigsq)*crossprod(xi,yi-mui);
      scorei <- hik - (1/numSubjects)*penalty%*%theta;
      I3 <- I3 + crossprod(t(scorei)); 
    }
    #print("I3");
    #print(str(scorei)); 
    #print(str(I3)); 
    # Calculate sandwichCovarianceTheta (in theta)
    if (kappa(I1)<Inf) {
      sandwichCovarianceTheta <-  solve(I1)%*%I3%*%solve(I1);  
      
    } else {
      sandwichCovarianceTheta <- matrix(NA,nrow(I1),ncol(I1));
    } 
    return(sandwichCovarianceTheta); 
  }
  #####################################################################
  ## Main body of TVEM function
  #print("Hello 1");
  ## Do listwise deletion; 
  if (!is.null(xcov)) {
    if (is.vector(xcov)) {xcov <- as.matrix(xcov);
    }
  }
  if (is.vector(tcov)) {tcov <- as.matrix(tcov);}
  NeedToDelete <- (apply(as.matrix(is.na(tcov)),1,sum)+
                     apply(as.matrix(is.na(xcov)),1,sum)+
                     is.na(dep)+is.na(time)+ 
                     is.na(id)) > 0;
  if (sum(NeedToDelete)>0) {
    if (deleteListwise==TRUE) {
      RowsNotToDelete <- (1:length(dep))[which(NeedToDelete==FALSE)];
      tcov <- tcov[RowsNotToDelete,];
      xcov <- xcov[RowsNotToDelete,];
      dep <- dep[RowsNotToDelete];
      time <- time[RowsNotToDelete];
      id <- id[RowsNotToDelete];
    } else {
      stop("Missing data detected");
    }
  }
  ## Process the subject ID's;
  if (!is.integer(id)) {
    temporary.id <- as.integer(as.factor(id));
  } else {
    temporary.id <- id;
  }
  if (is.na(min.time)) {min.time <- min(time);}
  if (is.na(max.time)) {max.time <- max(time);}
  numTotal <- length(id);
  intId <- rep(0,numTotal);
  stopifnot(is.integer(unique(temporary.id)));
  numSubjects <- length(unique(temporary.id));
  for (i in 1:numSubjects) {
    intId[which(temporary.id==unique(temporary.id)[i])] <- i;
  }
  #print("Hello 2");
  stopifnot(length(unique(intId))==numSubjects);
  stopifnot(all.equal(unique(intId),1:numSubjects));
  ## Process the time variable;
  stopifnot(!is.null(time));
  stopifnot((deg==1)|(deg==2)|(deg==3));
  time <- as.matrix(time);
  stopifnot(nrow(time)==numTotal);
  stopifnot(ncol(time)==1);
  stopifnot(is.numeric(time));
  if (is.null(colnames(time))) {
    colnames(time) <- "Time";
  }
  timeBasis <- list();
  designMatrix <- NULL;
  whichBeta <-  NULL;
  ## Construct the first part of the design matrix including ... ;
  ## ... the intercept column if there are no time-varying covariates;
  if (is.null(tcov)) {    interceptColumn <- matrix(1,numTotal,1);
  colnames(interceptColumn) <- "Intercept";
  designMatrix <- cbind(designMatrix,interceptColumn);
  whichBeta <- c(whichBeta,0);
  }
  ## ... and the covariates without time-varying effects;
  if (!is.null(xcov)) { 
    if (min(apply(xcov,2,var))<1e-10) {
      stop("Please do not include a constant column in xcov.");
    }
    #print(head(xcov));
    stopifnot(nrow(xcov)==numTotal);
    numXCov <- ncol(xcov);
    if (is.null(colnames(xcov))) {
      colnames(xcov) <- paste("X",1:ncol(xcov),sep="");
    };
    designMatrix <- cbind(designMatrix,xcov);
  } else {
    numXCov <- 0;
  }
  whichBeta <- c(whichBeta,rep(0,numXCov));
  ## Get the tcov matrix ready;
  if (is.null(tcov)) {
    numTCov <- 0;
    stop(paste("tcov is null.",
               "No time-varying betas or effects are in the model."));
  } else { 
    stopifnot(nrow(tcov)==numTotal);
    numTCov <- ncol(tcov);
    if (is.null(colnames(tcov))) {
      colnames(tcov) <- paste("TV",1:ncol(tcov),sep="");
    };
    ## Create the scale vector to be used with the penalty for the time-varying
    ## covariates;
    tcov.scale <- apply(tcov,2,sd);
    if (!identical(as.numeric(as.vector(tcov[,1])),as.numeric(rep(1,numTotal)))) {
      warning("tcov does not seem to contain an intercept (trajectory) column.");
    } else {
      tcov.scale[1] <- 1; # do not scale intercept column;
    }
    if (length(tcov.scale)>1) {
      if (min(tcov.scale[-1])<1e-10) {
        stop(paste("Please include the intercept column as the first",
                   "column in tcov, and do not include any other columns",
                   "which are constant across all subjects and assessments."));
      }
    }
    ## Create the basis functions;
    if (length(numInteriorKnots)!=1) {stop();}
    stopifnot(numInteriorKnots>0);
    num.intervals <- numInteriorKnots+1;
    dx <- (max.time-min.time)/num.intervals;
    all.knot.locations <- seq(min.time-(deg)*dx,
                              max.time+(deg)*dx,
                              by=dx);
    all.knot.locations[1] <- all.knot.locations[2]-1e-8;
    all.knot.locations[length(all.knot.locations)] <-
      all.knot.locations[length(all.knot.locations)-1]+1e-8;
    interior.knot.locations <- seq(min.time+dx,
                                   max.time-dx,
                                   by=dx);
    timeGrid <- seq(min.time,max.time,length=gridSize);
    timeBasis <- spline.des(all.knot.locations,
                            time,
                            deg+1,
                            rep(0,length(time)),
                            outer.ok=TRUE)$design;
    #print("Hello 3");
    # Adapted from Eilers and Marx, 1996;
    colnames(timeBasis) <- paste(colnames(time),".Spline.",
                                 1:ncol(timeBasis),
                                 sep="");
    timeBasisByGrid <- spline.des(all.knot.locations,
                                  timeGrid,
                                  deg+1,
                                  rep(0,length(timeGrid)),
                                  outer.ok=TRUE)$design;
    # Adapted from Eilers and Marx, 1996;
    colnames(timeBasisByGrid) <- paste(colnames(time),".Spline.",
                                       1:ncol(timeBasisByGrid),
                                       sep="");
    # Now generate regression matrix;
    for (j in 1:numTCov) {
      covariateTimesTimeBasis <- tcov[,j,drop=TRUE]*timeBasis;
      # Elementwise product;
      colnames(covariateTimesTimeBasis) <- paste(colnames(tcov)[j],
                                                 "times",
                                                 colnames(timeBasis),
                                                 sep=".");
      designMatrix <- cbind(designMatrix, covariateTimesTimeBasis);
      whichBeta <- c(whichBeta,rep(j,ncol(covariateTimesTimeBasis)));
    }
    stopifnot(sum(is.nan(designMatrix))==0);
    stopifnot(sum(is.null(designMatrix))==0);
    if(sum(is.na(designMatrix))>0) {
      stop("Missing data is not yet supported in this software.");
    }
  }
  ## Create the penalty weight matrix (see Eilers and Marx, 1996)
  diffMatrix <- crossprod(diff(diff(diag(rep(1,ncol(timeBasis))))));
  penaltyMatrix <- matrix(0,ncol(designMatrix),ncol(designMatrix));
  for (j in 1:numTCov) {
    indices <- numXCov+(j-1)*ncol(timeBasis)+(1:ncol(timeBasis));
    stopifnot(length(indices)==nrow(diffMatrix));
    stopifnot(length(indices)==ncol(diffMatrix));
    penaltyMatrix[indices,indices] <- diffMatrix/(tcov.scale[j]^2);
  }
  #######################################################################
  ## Find the best tuning parameter;
  f <- function(lambda) {
    #print(paste("Trying lambda=",round(lambda,5)));
    thisFit <- PenGEE(convergenceCriterion=convergenceCriterion,
                      dist=dist,
                      intId=intId,
                      maxIterations=maxIterations, 
                      penaltyMatrix=penaltyMatrix,
                      roughnessPenalty=lambda, 
                      time=time,
                      X=designMatrix,
                      Y=dep);
    if (useAIC==TRUE) { # Rather sloppy call to a global parameter;
      return(thisFit$aic);
    } else {
      return(thisFit$bic);
    }
  }
  lambda.search.results <- NULL;
  # Find the right order of magnitude;
  lambdasSearch1 <- 10^(seq(-4,ceiling(log10(numTotal))+1))*sd(dep);
  resultsSearch1 <- sapply(lambdasSearch1, f);
  bestLambda1 <- lambdasSearch1[which.min(resultsSearch1)];
  lambda.search.results <- rbind(lambda.search.results,
                             (cbind(lambdasSearch1,resultsSearch1)));
  # Find the right approximate value;
  lambdasSearch2 <- (10^seq(-1,1,by=.2))*bestLambda1;
  resultsSearch2 <- sapply(lambdasSearch2, f);
  bestLambda2 <- lambdasSearch2[which.min(resultsSearch2)];
  lambda.search.results <- rbind(lambda.search.results,
                             cbind(lambdasSearch2,resultsSearch2));
  # Zoom in closer;
  lambdasSearch3 <- (10^seq(-0.2,0.2,by=.02))*bestLambda2;
  resultsSearch3 <- sapply(lambdasSearch3, f);
  bestLambda3 <- lambdasSearch3[which.min(resultsSearch3)];
  lambda.search.results <- rbind(lambda.search.results,
                             cbind(lambdasSearch3,resultsSearch3));
  # Zoom in closer;
  lambdasSearch4 <- (10^seq(-0.02,0.02,by=.005))*bestLambda3;
  resultsSearch4 <- sapply(lambdasSearch4, f);
  bestLambda4 <- lambdasSearch4[which.min(resultsSearch4)];
  lambda.search.results <- rbind(lambda.search.results,
                             cbind(lambdasSearch4,resultsSearch4));
  # Good enough;
  lambda <- bestLambda4;
  ## Finally do the analysis;
  bestFit <- PenGEE(convergenceCriterion=convergenceCriterion,
                    dist=dist,
                    intId=intId,
                    maxIterations=maxIterations,
                    penaltyMatrix=penaltyMatrix,
                    roughnessPenalty=lambda,
                    time=time,
                    X=designMatrix,
                    Y=dep);
  fittedCoefficients <- list();
  for (j in 1:numTCov) {
    fittedCoefficients[[j]] <- timeBasis%*%
      bestFit$theta[which(whichBeta==j),];
  }
  estimate <- list();
  for (which.coef in 1:numTCov) {
    estimate[[which.coef]] <- fittedCoefficients[[which.coef]][order(time)];
  }
  ## Record the fitted values on the fit grid;
  fittedCoefficientsByGrid <- list();
  for (j in 1:numTCov) {
    fittedCoefficientsByGrid[[j]] <- timeBasisByGrid%*%
      bestFit$theta[which(whichBeta==j),];
  }
  estimate <- list();
  for (which.coef in 1:numTCov) {
    estimate[[which.coef]] <- list(); 
    estimate[[which.coef]] <-
      fittedCoefficientsByGrid[[which.coef]][order(timeGrid) ];
  }
  ## Get the standard errors;
  SandwichCovTheta <- PenGEECovMats(dist=dist,
                                    intId=intId,
                                    mu=bestFit$fittedY,
                                    numSubjects=numSubjects,
                                    numTotal=numTotal,
                                    penalty=lambda*penaltyMatrix,
                                    sigsq=ifelse(dist=="normal",bestFit$sigsq,1),
                                    theta=bestFit$theta,
                                    X=designMatrix,
                                    Y=dep);
  if (sum(is.na(SandwichCovTheta))>0) {
    warning("Could not compute sandwich standard errors.");
  }
  fittedCoefficientsStdErrSandwich <- list();
  fittedCoefficientsStdErrSandwichByGrid <- list(); 
  for (j in 1:numTCov) {
    thisCovmatSandwich <-
      SandwichCovTheta [which(whichBeta==j),which(whichBeta==j),drop=FALSE];
    fittedCoefficientsStdErrSandwich[[j]] <- sqrt(diag(
      as.matrix(timeBasis %*% thisCovmatSandwich %*% t(timeBasis))));
    fittedCoefficientsStdErrSandwichByGrid[[j]] <- sqrt(diag(
      as.matrix(timeBasisByGrid %*% thisCovmatSandwich %*% t(timeBasisByGrid))));
  }
  if (doPlot==TRUE) {
    stopifnot(numTCov<4);
    if (numTCov==2) {par(mfrow=c(1,2),
                         mar=c(4,4,4,2));}
    if (numTCov>2) {par(mfrow=c(2,2),
                        mar=c(4,4,4,2));}
    for (which.coef in 1:numTCov) {
      plot(x=timeGrid[order(timeGrid)],
           y=fittedCoefficientsByGrid[[which.coef]][order(timeGrid),],
           ylim=c(min(c(0,unlist(fittedCoefficientsByGrid[[which.coef]]))),
                  max(c(0,unlist(fittedCoefficientsByGrid[[which.coef]])))),
           type="l",
           xlab="Time",
           ylab=paste("Coefficient",which.coef));
      abline(h=0,lty="dashed");
      if (sum(is.na(fittedCoefficientsStdErrSandwichByGrid[[which.coef]]))==0) {
        lines(x=timeGrid[order(timeGrid)],
              y=fittedCoefficientsByGrid[[which.coef]][order(timeGrid)]+1.96*
                fittedCoefficientsStdErrSandwichByGrid[[which.coef]][order(timeGrid)],
              lty="dotted");
        lines(x=timeGrid[order(timeGrid)],
              y=fittedCoefficientsByGrid[[which.coef]][order(timeGrid)]-1.96*
                fittedCoefficientsStdErrSandwichByGrid[[which.coef]][order(timeGrid)],
              lty="dotted");
      }
    }
  }
  cat("TVEM R Function \n");
  cat(sprintf("Number of subjects: %29.0f \n",numSubjects));
  cat(sprintf("Total number of observations: %19.0f \n",numTotal));
  if (dist=="normal") {cat("Distribution type:  Normal")};
  if (dist=="logistic") {cat("Distribution type:  Binary, logistic link")};
  if (dist=="poisson") {cat("Distribution type:  Poisson, log link")};
  if(deg==1) {cat("Coefficient function change between knots treated as linear \n")};
  if(deg==2) {cat("Coefficient function change between knots treated as quadratic \n")};
  if(deg==3) {cat("Coefficient function change between knots treated as cubic \n")};
  if (!bestFit$converged) {
    warning("Penalized iteratively reweighted least squares algorithm did not converge");
  }
  cat(sprintf("Roughness penalty weight: %23.2f \n", lambda));
  cat(sprintf("Log pseudo-likelihood: %33.2f \n",bestFit$logLik));
  cat(sprintf("Pseudolikelihood AIC: %44.2f \n",bestFit$aic));
  cat(sprintf("Pseudolikelihood BIC: %44.2f \n",bestFit$bic));
  cat(sprintf("Count number of parameters: %21.2f \n", bestFit$np));
  cat(sprintf("Smoothed number of parameters: %18.2f \n", bestFit$enp));
  if(dist=="normal") { 
    cat(sprintf("RSS statistic: %25.2f \n", bestFit$RSS));cat("Standard deviation: \n");
    cat(round(sqrt(bestFit$sigsq),4)); cat("\n");
  }
  nontvem.reg.output <- NULL;
  if (!is.null(xcov)) {
    xcov.est <- bestFit$theta[which(whichBeta==0)];
    these.indices <- which(whichBeta==0);
    xcov.se.sandwich <- as.vector(sqrt(diag(SandwichCovTheta[these.indices, these.indices,drop=FALSE])));
    xcov.z.sandwich <- xcov.est/xcov.se.sandwich;
    xcov.p.sandwich <- 2*(1-pnorm(abs(xcov.z.sandwich)));
    nontvem.reg.output <- data.frame(Column=colnames(xcov),
                                     Estimate=xcov.est,
                                     SESandwich=xcov.se.sandwich,
                                     z.sandwich=round(xcov.z.sandwich,8),
                                     p.sandwich=round(xcov.p.sandwich,8));
    cat("Non-time-varying coefficients: \n");
    print(nontvem.reg.output);
  }
  # Return the answers;
  # Should we also return exp(beta) in addition to beta to get odds ratios in the binary case? 
  # Return the answers; 
  return(list( bestFit=bestFit,
               # Answer from the TVEMFitInner function,
               # as obtained for the best value of lambda.
               allBSplineKnots=all.knot.locations,
               # List of all knot locations for the spline
               beta=fittedCoefficients,
               # The fitted values for the coefficient functions
               # at each assessment time value
               betaByGrid=fittedCoefficientsByGrid,
               # The fitted values for the coefficient functions
               # at time value on a regular grid of GridSize points
               betaSESandwich=fittedCoefficientsStdErrSandwich,
               # The standard errors for the coefficient functions
               # at each assessment time value
               betaSESandwichByGrid=fittedCoefficientsStdErrSandwichByGrid,
               # The standard errors for the coefficient functions
               # at time value on a regular grid of GridSize points
               fittedValues=fitted,
               # The fitted values for the dependent variable
               # for each class, at each assessment time value.
               # Also included as a member of bestFit. 
               penalty=lambda*penaltyMatrix,
               timeBasis=timeBasis,
               # The basis matrix used to express the effect
               # of time.
               timeGrid=timeGrid,
               technical=list(LambdaSearchResults=lambda.search.results,
                              interiorKnots=interior.knot.locations)));
  # The time values used to arrange predictions;
  # on a regular set of points for graphing.    ;    
}