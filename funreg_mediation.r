library(mgcv);
funreg_mediation <- function(data,
                             treatment,
                             mediator,
                             outcome,
                             id,
                             time,
                             tve_covariates_on_mediator=NULL,
                             tie_covariates_on_mediator=NULL,
                             covariates_on_outcome=NULL,
                             logistic=FALSE, # FALSE for numerical outcome, TRUE for dichotomous 0/1;
                             nboot=49,
                             num_knots=20,  
                             # number of interior knots per function (not counting exterior knots);
                             spline_order=3,
                             # shape of function between knots, with default of 3 for cubic;
                             penalty_function_order=2,
                             smoothing_penalty_multiplier=1,
                             grid=100
) {    
  #-------------------------------------------;
  #--- PROCESSING OF INPUT -------------------;
  #-------------------------------------------;
  cat("Working")
  m <- match.call(expand.dots = FALSE);
  m$treatment <- NULL;
  m$mediator <- NULL;
  m$outcome <- NULL;
  m$tve_covariates_on_mediator <- NULL;
  m$tie_covariates_on_mediator <- NULL;
  m$covariates_on_outcome <- NULL;  
  m$logistic <- NULL;
  m$grid <- NULL;
  m$nboot <- NULL;
  if (is.matrix(eval.parent(m$data))) {
    m$data <- as.data.frame(data);
  }
  m[[1]] <- quote(stats::model.frame);
  m <- eval.parent(m);
  id_variable_name <- as.character(substitute(id));
  id_variable <- m[,id_variable_name];
  time_variable_name <- as.character(substitute(time));
  time_variable <- m[,time_variable_name]; 
  treatment_variable_name <- as.character(substitute(treatment));
  treatment_variable <- m[,treatment_variable_name]; 
  mediator_variable_name <- as.character(substitute(mediator));
  mediator_variable <- m[,mediator_variable_name]; 
  outcome_variable_name <- as.character(substitute(outcome));
  outcome_variable <- m[,outcome_variable_name]; 
  if (length(num_knots)>1) {
    stop("Please specify a single number for num_knots.")
  };   
  long_data_for_analysis <- as.data.frame(eval.parent(data));
  long_data_for_analysis$id_for_tvem_function <- id_variable;
  long_data_for_analysis$time_for_tvem_function <- time_variable;
  #-------------------------------------------;
  # ----- Process covariates -----------------;
  if (is.null (covariates_on_outcome)) {
    covariates_on_outcome_names <- NULL;
  } else {
    covariates_on_outcome_names <- attr(terms(covariates_on_outcome),"term.labels");
  }
  num_covariates_on_outcome <- length(covariates_on_outcome_names);
  covariates_on_outcome_data <- long_data_for_analysis[,covariates_on_outcome_names,drop=FALSE];
  if (is.null (tve_covariates_on_mediator)) {
    tve_covariates_on_mediator_names <- NULL;
  } else {
    tve_covariates_on_mediator_names <- attr(terms(tve_covariates_on_mediator),"term.labels");
  }
  num_tve_covariates_on_mediator <- length(tve_covariates_on_mediator_names); 
  tve_covariates_on_mediator_data <- long_data_for_analysis[,tve_covariates_on_mediator_names,drop=FALSE];
  if (is.null (tie_covariates_on_mediator)) {
    tie_covariates_on_mediator_names <- NULL;
  } else {
    tie_covariates_on_mediator_names <- attr(terms(tie_covariates_on_mediator),"term.labels");
  }
  num_tie_covariates_on_mediator <- length(tie_covariates_on_mediator_names);
  tie_covariates_on_mediator_data <- long_data_for_analysis[,tie_covariates_on_mediator_names,drop=FALSE];
  if (length(intersect(tve_covariates_on_mediator_names, tie_covariates_on_mediator_names))>0) {
    stop(paste("Please make sure that no covariate is specified as having",
               "both time-varying and time-invariant effects."));    
  }
  #-------------------------------------------;
  #--- CONVERT MEDIATOR M TO WIDE FORM -------;
  #-------------------------------------------;  
  if (max(unlist(lapply(split(treatment_variable,f=id_variable),var)))>1e-10) {
    stop("Please make sure that the subject-level treatment is constant within subject.")
  }
  if (max(unlist(lapply(split(outcome_variable,f=id_variable),var)))>1e-10) {
    stop("Please make sure that the subject-level outcome is constant within subject.")
  }
  observed_time_grid <- sort(unique(time_variable));
  wide_id <- sort(unique(id_variable));
  # Should the user provide the data as long or wide?
  nobs <- length(observed_time_grid);
  nsub <- length(wide_id);
  wide_mediator <- matrix(NA,nsub,nobs);
  wide_treatment <- rep(NA,nsub);
  wide_outcome <- rep(NA,nsub);
  if (num_covariates_on_outcome>0) {
    wide_covariates_on_outcome <- matrix(NA,nsub,num_covariates_on_outcome);
  }
  if (num_tve_covariates_on_mediator>0) {
    wide_tve_covariates_on_mediator <- matrix(NA,nsub,num_tve_covariates_on_mediator);
  }
  if (num_tie_covariates_on_mediator>0) {
    wide_tie_covariates_on_mediator <- matrix(NA,nsub,num_tie_covariates_on_mediator);
  }
  for (this_id in 1:length(wide_id)) {
    these_rows <- which(id_variable==wide_id[this_id]); 
    if (length(these_rows)>0) {
      if (min(treatment_variable[these_rows])!=
          max(treatment_variable[these_rows])) {
        stop("Please make sure the treatment is the same for each observation within subject.")
      }
      wide_treatment[this_id] <- treatment_variable[min(these_rows)];
      if (min(outcome_variable[these_rows])!=
          max(outcome_variable[these_rows])) {
        stop("Please make sure the outcome is the same for each observation within subject.")
      }
      wide_outcome[this_id] <- outcome_variable[min(these_rows)];
      if (num_covariates_on_outcome>0) {
        if (max(apply(covariates_on_outcome_data[these_rows,,drop=FALSE],2,var))>0) {
          stop("Please make sure that the covariates on the outcome do not vary within subject for this model.")
        }
      }
      if (num_tve_covariates_on_mediator>0) {
        if (max(apply(tve_covariates_on_mediator_data[these_rows,,drop=FALSE],2,var))>0) {
          stop("Please make sure that the covariates on the mediator do not vary within subject for this model.")
        }
      }
      if (num_tie_covariates_on_mediator>0) {
        if (max(apply(tie_covariates_on_mediator_data[these_rows,,drop=FALSE],2,var))>0) {
          stop("Please make sure that the covariates on the mediator do not vary within subject for this model.")
        }
      }
    }
    for (this_time in 1:length(observed_time_grid)) { 
      this_data <- which(id_variable==wide_id[this_id] &
                           time_variable==observed_time_grid[this_time]); 
      if (length(this_data)==1) {
        wide_mediator[this_id,this_time] <- long_simulated_data$M[this_data]; 
        if (num_covariates_on_outcome>0) {
          wide_covariates_on_outcome[this_id,] <- as.matrix(covariates_on_outcome_data[this_data,,drop=FALSE]);
        } 
        if (num_tve_covariates_on_mediator>0) {
          wide_tve_covariates_on_mediator[this_id,] <- as.matrix(tve_covariates_on_mediator_data[this_data,,drop=FALSE]);
        } 
        if (num_tie_covariates_on_mediator>0) {
          wide_tie_covariates_on_mediator[this_id,] <- as.matrix(tie_covariates_on_mediator_data[this_data,,drop=FALSE]);
        } 
      }
      if (length(this_data)==2) {
        stop("There seems to be more than one measurement on the same subject and time.")
      }
    }
  }  
  wide_data <- data.frame(cbind(wide_id=wide_id, 
                                wide_treatment=wide_treatment, 
                                wide_outcome=wide_outcome,
                                wide_mediator=wide_mediator));
  wide_id_column <- 1;
  wide_treatment_column <- 2;
  wide_outcome_column <- 3;
  wide_mediator_columns <- 4:ncol(wide_data); 
  if (num_covariates_on_outcome>0) {
    wide_data <- data.frame(cbind(wide_data,
                                  wide_covariates_on_outcome));
    wide_covariates_on_outcome_columns <- (ncol(wide_data)-num_covariates_on_outcome+1):ncol(wide_data);
  }
  if (num_tve_covariates_on_mediator>0) {
    wide_data <- cbind(wide_data,
                       wide_tve_covariates_on_mediator);
    wide_tve_covariates_on_mediator_columns <- (ncol(wide_data)-num_tve_covariates_on_mediator+1):ncol(wide_data);
  }
  if (num_tie_covariates_on_mediator>0) {
    wide_data <- data.frame(cbind(wide_data,
                                  wide_tie_covariates_on_mediator))
    wide_tie_covariates_on_mediator_columns <- (ncol(wide_data)-num_tie_covariates_on_mediator+1):ncol(wide_data);
  }
  #-------------------------------------------;
  #--- DEFINE MAIN ANALYSIS ------------------;
  #-------------------------------------------;  
  analyze_data_for_mediation <- function(wide_data, indices, get_details=FALSE) {
    local_wide_data <- wide_data[indices,];
    cat(".");
    #--- Take data frame apart into pieces to use with pfr function 
    wide_mediator <- local_wide_data[,wide_mediator_columns];
    wide_id <- local_wide_data[,wide_id_column];
    wide_treatment <- local_wide_data[,wide_treatment_column]; 
    wide_outcome <- local_wide_data[,wide_outcome_column];
    if (num_covariates_on_outcome>0) { 
      wide_covariates_on_outcome <- local_wide_data[,wide_covariates_on_outcome_columns,drop=FALSE];
     }  
    if (num_tve_covariates_on_mediator>0) { 
      wide_tve_covariates_on_mediator <- local_wide_data[,wide_tve_covariates_on_mediator_columns,drop=FALSE];
      }  
    if (num_tie_covariates_on_mediator>0) { 
      wide_tie_covariates_on_mediator <- local_wide_data[,wide_tie_covariates_on_mediator_columns,drop=FALSE];
      } 
    nobs <- length(wide_mediator_columns);
    nsub <- length(wide_id);
    #--- EFFECT OF MEDIATOR M AND TREATMENT X ON OUTCOME Y ---;
    pfr_formula <- wide_outcome~lf(wide_mediator,
                                   presmooth="interpolate")+wide_treatment;
    for (this_one in 1:num_covariates_on_outcome) {  
      assign(paste("wide_",covariates_on_outcome_names[this_one],sep=""),
             wide_covariates_on_outcome[,this_one]);
      new_one <- as.formula(paste("~.+wide_",covariates_on_outcome_names[this_one],sep=""));
      pfr_formula <- update(pfr_formula,new_one);
    }
    if (logistic) {
      funreg_MY <- pfr(pfr_formula,
                       family=binomial());
      # I wish I could let the user send the data and family in from outside the function,
      # but the pfr function does not allow this due to its unusual implementation
      # as a wrap-around for a hidden call to gam.
    } else {
      funreg_MY <- pfr(pfr_formula,
                       family=gaussian());
    } 
    alpha_int_estimate <- as.numeric(funreg_MY$coefficient["(Intercept)"]);
    alpha_int_se <- as.numeric(summary(funreg_MY)$se["(Intercept)"]);
    alpha_X_estimate <-  as.numeric(funreg_MY$coefficient["wide_treatment"]);
    alpha_X_se <- as.numeric(summary(funreg_MY)$se["wide_treatment"]);
    alpha_M_estimate <- as.numeric(coef(funreg_MY)[,"value"]); 
    alpha_M_se <- coef(funreg_MY)[,"se"]; 
    alpha_M_true <- alpha.M(time.grid);
    alpha_M_pvalue <- summary(funreg_MY)$s.table[1,"p-value"];
    #--- DIRECT EFFECT OF TREATMENT X ON OUTCOME Y ---;
    glm_formula <- wide_outcome~wide_treatment;
    for (this_one in 1:num_covariates_on_outcome) { 
      new_one <- as.formula(paste("~.+wide_",covariates_on_outcome_names[this_one],sep=""));
      glm_formula <- update(glm_formula,new_one);
    } 
    if (logistic) {
      model_for_direct_effect_XY <- glm(glm_formula,
                                        family=binomial);
    } else {
      model_for_direct_effect_XY <- glm(glm_formula);
    }
    delta_int_estimate <- as.numeric(model_for_direct_effect_XY$coefficients["(Intercept)"]);
    delta_int_se <- summary(model_for_direct_effect_XY)$coefficients["(Intercept)","Std. Error"];
    delta_X_estimate <- as.numeric(model_for_direct_effect_XY$coefficients["wide_treatment"]);
    delta_X_se <- summary(model_for_direct_effect_XY)$coefficients["wide_treatment","Std. Error"];
    #--- EFFECT OF TREATMENT X ON MEDIATOR M ---;   
    local_long_data <- data.frame(id=rep(wide_id, each=nobs),
                                  time=rep(observed_time_grid, times=nsub),
                                  outcome=rep(wide_outcome, each=nobs),
                                  treatment=rep(wide_treatment, each=nobs),
                                  mediator=as.vector(t(wide_mediator)));
    tvem_formula1 <- mediator~treatment;
    tvem_formula2 <- tie_covariates_on_mediator;
    for (this_one in 1:num_tve_covariates_on_mediator) {
      new_name <- tve_covariates_on_mediator_names[this_one];
      temp_data_frame <- data.frame(temp=rep(wide_tve_covariates_on_mediator[,this_one],each=nobs));
      colnames(temp_data_frame) <- new_name;
      local_long_data <- cbind(local_long_data, temp_data_frame); 
      new_term <- as.formula(paste("~.+",new_name,sep=""));
      tvem_formula1 <- update(tvem_formula1,new_term);
    }
    for (this_one in 1:num_tie_covariates_on_mediator) {
      new_name <- tie_covariates_on_mediator_names[this_one];
      temp_data_frame <- data.frame(temp=rep(wide_tie_covariates_on_mediator[,this_one],each=nobs));
      colnames(temp_data_frame) <- new_name;
      local_long_data <- cbind(local_long_data, temp_data_frame); 
    }
    local_long_data <- local_long_data[which(!is.na(local_long_data$mediator)),];
       # listwise deletion to remove empty observations; 
    tvem_XM <- tvem(data=local_long_data,
                    formula=tvem_formula1,
                    time=time,
                    id=id,
                    invar_effects=tvem_formula2,
                    grid=observed_time_grid); 
    gamma_int_estimate <- tvem_XM$grid_fitted_coefficients[[1]]$estimate; 
    gamma_int_se <- tvem_XM$grid_fitted_coefficients[[1]]$standard_error;
    gamma_X_estimate <- tvem_XM$grid_fitted_coefficients[[2]]$estimate; 
    gamma_X_se <- tvem_XM$grid_fitted_coefficients[[2]]$standard_error; 
    # #--- MEDIATED EFFECT OF TREATMENT X THROUGH MEDIATOR M ON OUTCOME Y ---;
    beta_estimate <- mean(gamma_X_estimate*alpha_M_estimate);
    if (get_details) {
      return(list(time_grid=observed_time_grid,
                  gamma_int_estimate=gamma_int_estimate,
                  gamma_int_se=gamma_int_se,
                  gamma_X_estimate=gamma_X_estimate,
                  gamma_X_se=gamma_X_se,
                  alpha_int_estimate=alpha_int_estimate,  
                  alpha_int_se=alpha_int_se, 
                  alpha_X_estimate=alpha_X_estimate, 
                  alpha_X_se=alpha_X_se,
                  alpha_M_estimate=alpha_M_estimate,   
                  alpha_M_se=alpha_M_se,   
                  alpha_M_pvalue=alpha_M_pvalue,
                  delta_int_estimate=delta_int_estimate,
                  delta_int_se=delta_int_se,
                  delta_X_estimate=delta_X_estimate,
                  delta_X_se=delta_X_se,
                  beta_estimate=beta_estimate,
                  tvem_XM_details=tvem_XM,
                  funreg_MY_details=funreg_MY,
                  direct_effect_details=model_for_direct_effect_XY));
    } else {
      return(beta_estimate);
    }
  } 
  #------------------------------------------;
  #---------------- Call function -----------;
  #------------------------------------------;
  original_results <- analyze_data_for_mediation(wide_data,
                                                 indices=1:nrow(wide_data),
                                                 get_details=TRUE); 
  before_boot <- Sys.time();  
  boot1 <- boot(data=wide_data,
                statistic=analyze_data_for_mediation, 
                R=nboot); 
  boot2 <- boot.ci(boot1,conf=.95,type="norm");
  boot3 <- boot.ci(boot1,conf=.95,type="basic");
  boot4 <- boot.ci(boot1,conf=.95,type="perc");
  after_boot <- Sys.time(); 
  print(paste("Bootstrapping completed in", 
              round(difftime(after_boot,before_boot,units="mins"),2),
              "minutes."));
  bootstrap_results <- list(beta.boot.estimate= norm.ci(boot1,conf=.001)[2],
                            beta.boot.se=sd(boot1$t),
                            beta.boot.norm.lower=boot2$normal[2],
                            beta.boot.norm.upper=boot2$normal[3],
                            beta.boot.basic.lower=boot3$basic[4],
                            beta.boot.basic.upper=boot3$basic[5],
                            beta.boot.perc.lower=boot4$percent[4],
                            beta.boot.perc.upper=boot4$percent[5],
                            boot1=boot1,
                            time.required=difftime(after_boot,before_boot));
  answer <- list(original_results=original_results,
                 bootstrap_results=bootstrap_results);
  class(answer) <- "funreg_mediation";
  return(answer);
}  