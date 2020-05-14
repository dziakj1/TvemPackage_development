#' funreg_mediation:  Fit funreg_mediation model.
#' 
#' Calculate indirect effect of a binary treatment on a scalar  
#' response as mediated by a longitudinal functional trajectory
#' (see Baron & Kenny, 1986; Lindquist, 2012; Coffman et al., 2019).
#' 
#' @note This function calls the tvem function in this package.
#' It also calls the pfr function in the refund package (see
#' Goldsmith et al., 2011) to perform penalized functional regression.
#' Some suggestions on interpreting the output from penalized functional
#' regression are given by Dziak et al. (2019). 
#' 
#' @references
#' Baron, R.M., & Kenny, D.A. (1986). The moderator-mediator variable
#' distinction in social psychological research: Conceptual, strategic, 
#' and statistical considerations. Journal of Personality & Social 
#' Psychology, 51: 1173-1182.
#' Coffman, D. L., Dziak, J. J., Li, R., & Piper, M. (2019). Functional 
#' regression mediation analysis with application to a smoking 
#' cessation intervention. Joint Statistical Meetings presentation, 
#' August 2019.
#' Dziak, J. J., Coffman, D. L., Reimherr, M., Petrovich, J., Li, R., 
#' Shiffman, S., & Shiyko, M. P. (2019). Scalar-on-function regression
#' for predicting distal outcomes from intensively gathered longitudinal
#' data: interpretability for applied scientists. Online publication, 
#' Statistics Surveys.
#' Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich, D.
#' (2011). Penalized functional regression. Journal of Computational 
#' and Graphical Statistics, 20(4), 830-851.
#' Lindquist, M. A. (2012). Functional Causal Mediation Analysis 
#' With an Application to Brain Connectivity. Journal of the American
#' Statistical Association, 107: 1297-1309.
#' 
#' @param data The dataset containing the data to be analyzed, in long format (one row per observation, multiple per individual).
#' @param treatment  The name of the variable containing the treatment assignment, assumed to be unidimensional (dichotomous or numerical; we recommend dichotomous with 0 for control and 1 for experimental).  The values of this variable should be the same for each row for a given subject.
#' @param mediator The name of the mediator variable. The values of this variable can (and should) vary within each subject.
#' @param outcome The name of the outcome variable. The values of this variable should be the same for each row for a given subject.
#' @param id The name of the variable identifying each subject.
#' @param time The name of the time variable. 
#' @param tve_covariates_on_mediator The time-varying-effects covariates, if any, to be included in the 
#' model predicting the mediator from the treatment.
#' @param tie_covariates_on_mediator The non-time-varying-effects covariates, if any, to be included in the 
#' model predicting the mediator from the treatment.
#' @param covariates_on_outcome The covariates, if any, to be included in the model predicting the outcome from the treatment.
#' They are assumed to be subject-level (non-time-varying both in value and in effect).
#' @param logistic Whether the outcome should be modeled as dichotomous with a logistic model (TRUE), 
#' or numerical with a normal model (FALSE). So far these are the only options we have implemented.
#' @param nboot Number of bootstrap samples for bootstrap significance test of the overall effect. This
#' test is done using the boot function from the boot package by Angelo Canty and Brian Ripley. 
#' It differs somewhat from the bootstrap approach used in a similar context by Lindquist (2012). 
#' @param boot_level One minus the nominal coverage to be attempted for the bootstrap confidence interval estimates.
#' @param tvem_spline_order Input to be passed on to the TVEM function
#' @param tvem_penalty_order Input to be passed on to the TVEM function
#' @param tvem_penalize Input to be passed on to the TVEM function
#' @param tvem_do_loop Whether to use a loop to select the number of knots with a pseudo-AIC or pseudo-BIC
#' @param tvem_num_knots If tvem_do_loop is FALSE, then tvem_num_knots is passed on to the tvem function as num_knots, 
#' the number of interior knots for the B-splines. If tvem_do_loop is TRUE then tvem_num_knots is reinterpreted as the highest number of interior knots to try. This can be either 
#' an integer or a vector if tvem_do_loop is FALSE, but must be an integer if tvem_do_loop is TRUE. 
#' @param tvem_use_bic This parameter only matters if tvem_do_loop is TRUE. If tvem_do_loop is TRUE
#' and tvem_use_bic is TRUE, then the information criterion used will be a pseudolikelihood version of BIC. 
#' If tvem_do_loop is TRUE and tvem_use_bic is FALSE, then the information criterion used will be
#' a pseudolikelihood version of AIC instead. If tvem_do_loop is FALSE then tvem_use_bic is ignored.
#' 
#' @return An object of type funreg_mediation. The components of an object of 
#' type funreg_mediation are as follows:
#' \describe{
#' \item{original_results}{The estimates from the fitted models for
#' predicting the mediator from the treatment, predicting the outcome from 
#' the mediator and treatment, and predicting the outcome from the treatment alone.}
#' \item{bootstrap_results}{The estimate and confidence interval of the indirect effect
#' using a bootstrap approach.} 
#' } 
#' 
#' The original_results component has these components within it:
#' \describe{
#' \item{time_grid}{Grid of time points on which the functional coefficients are estimated.}
#' \item{gamma_int_estimate}{Estimated intercept function (as a vector of estimates) from the TVEM regression of the mediator on the treatment}
#' \item{gamma_int_se}{Estimated pointwise standard errors associated with the above}
#' \item{gamma_X_estimate}{Estimated time-varying treatment effect from the TVEM regression of the mediator on the treatment}
#' \item{gamma_X_se}{Estimated pointwise standard errors associated with the above}
#' \item{alpha_int_estimate}{Estimated scalar intercept from the scalar-on-function regression of the outcome on the mediator and treatment}
#' \item{alpha_int_se}{Estimated standard error for the above}
#' \item{alpha_X_estimate}{Estimated scalar coefficient for the treatment, from the scalar-on-function regression of the outcome on the mediator and treatment}
#' \item{alpha_X_se}{Estimated standard error for the above}
#' \item{alpha_M_estimate}{Estimated functional coefficient for the mediator,  from the scalar-on-function regression of the outcome on the mediator and treatment}
#' \item{alpha_M_se}{Estimated pointwise standard errors associated with the above}
#' \item{alpha_M_pvalue}{The p-value for significance of mediator in predicting outcome, after adjusting for treatment}
#' \item{delta_int_estimate}{Intercept from simple model predicting outcome directly from treatment}
#' \item{delta_int_se}{Estimated standard error for the above}
#' \item{delta_X_estimate}{Coefficient for treatment in model predicting outcome directly from treatment}
#' \item{delta_X_se}{Estimated standard error for the above}
#' \item{beta_estimate}{Estimated indirect effect, calculated as the dot product of the effect of treatment on mediator and the treatment-adjusted effect of mediator on outcome.  It is a scalar, even though they are functions of time.}
#' \item{tvem_XM_details}{Detailed output from the tvem function for the time-varying-effect model predicting the mediator from the treatment}
#' \item{funreg_MY_details}{Detailed output from the refund::pfr function for the scalar-on-function functional regression predicting the outcome from the treatment and mediator}
#' \item{direct_effect_details}{Detailed output from the linear or generalized linear model predicting the outcome from the treatment alone, ignoring the mediator (direct effect)}
#' }
#' 
#' The bootstrap_results component has these components within it:
#' \describe{
#' \item{beta_boot_estimate}{Bootstrap point estimate of the indirect effect (average of bootstrap sample estimates).} 
#' \item{beta_boot_se}{Bootstrap standard error for the indirect effect (standard deviation of bootstrap sample estimates).}
#' \item{beta_boot_norm_lower}{Lower end of the bootstrap confidence interval using the normal method in boot.ci inthe boot package.}
#' \item{beta_boot_norm_upper}{Upper end of the bootstrap confidence interval using the normal method.}
#' \item{beta_boot_basic_lower}{Lower end of the bootstrap confidence interval using the basic method in boot.ci in the boot package.}
#' \item{beta_boot_basic_upper}{Upper end of the bootstrap confidence interval using the basic method.}
#' \item{beta_boot_perc_lower}{Lower end of the bootstrap confidence interval using the percentile method in boot.ci inthe boot package.}
#' \item{beta_boot_perc_upper}{Upper end of the bootstrap confidence interval using the percentile method.}
#' \item{boot_level}{The alpha level used for the bootstrap confidence interval.}
#' \item{boot1}{The output returned from the boot function.}
#' \item{time.required}{The amount of time spent doing the bootstrap test, including generating and analyzing all samples.}
#' } 
#' @export


funreg_mediation <- function(data,
                             treatment,
                             mediator,
                             outcome,
                             id,
                             time,
                             tve_covariates_on_mediator=NULL,
                             tie_covariates_on_mediator=NULL,
                             covariates_on_outcome=NULL,
                             tvem_penalize=TRUE,
                             tvem_penalty_order=1,
                             tvem_spline_order=3,
                             tvem_num_knots=3, 
                             tvem_do_loop=FALSE,
                             tvem_use_bic=FALSE,
                             logistic=FALSE, # FALSE for numerical outcome, TRUE for dichotomous 0/1;
                             nboot=199,
                             boot_level=.05) {    
  #-------------------------------------------;
  #--- PROCESSING OF INPUT -------------------;
  #-------------------------------------------; 
  m <- match.call(expand.dots = FALSE);
  m$treatment <- NULL;
  m$mediator <- NULL;
  m$outcome <- NULL;
  m$tve_covariates_on_mediator <- NULL;
  m$tie_covariates_on_mediator <- NULL;
  m$covariates_on_outcome <- NULL;  
  m$logistic <- NULL;
  m$tvem_num_knots <- NULL;
  m$tvem_penalty_order <- NULL;
  m$tvem_spline_order <- NULL;
  m$tvem_penalize <- NULL;
  m$tvem_do_loop <- NULL;
  m$grid <- NULL;
  m$nboot <- NULL;
  if (is.matrix(eval.parent(m$data))) {
    m$data <- as.data.frame(data);
  }
  m[[1]] <- quote(stats::model.frame);
  m <- eval.parent(m);
  id_variable_name <- as.character(substitute(id));
  time_variable_name <- as.character(substitute(time));
  treatment_variable_name <- as.character(substitute(treatment));
  mediator_variable_name <- as.character(substitute(mediator));
  outcome_variable_name <- as.character(substitute(outcome));
  long_data_for_analysis <- as.data.frame(eval.parent(data)); 
  id_variable <- long_data_for_analysis[,id_variable_name];
  time_variable <- long_data_for_analysis[,time_variable_name]; 
  treatment_variable <- long_data_for_analysis[,treatment_variable_name]; 
  mediator_variable <- long_data_for_analysis[,mediator_variable_name]; 
  outcome_variable <- long_data_for_analysis[,outcome_variable_name]; 
  if (sum(is.na(id_variable>0))) {
    stop("Please remove data rows which have missing data on the id variable.");
  }
  #-------------------------------------------;
  # ----- Process covariates -----------------;
  if (is.null(covariates_on_outcome)) {
    covariates_on_outcome_names <- NULL;
  } else {
    covariates_on_outcome_names <- attr(terms(covariates_on_outcome),"term.labels");
  }
  num_covariates_on_outcome <- length(covariates_on_outcome_names);
  covariates_on_outcome_data <- long_data_for_analysis[,covariates_on_outcome_names,drop=FALSE];
  if (is.null(tve_covariates_on_mediator)) {
    tve_covariates_on_mediator_names <- NULL;
  } else {
    tve_covariates_on_mediator_names <- attr(terms(tve_covariates_on_mediator),"term.labels");
  }
  num_tve_covariates_on_mediator <- length(tve_covariates_on_mediator_names); 
  tve_covariates_on_mediator_data <- long_data_for_analysis[,tve_covariates_on_mediator_names,drop=FALSE];
  if (is.null(tie_covariates_on_mediator)) {
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
  if (max(unlist(lapply(split(treatment_variable,f=id_variable),var, na.rm=TRUE)))>1e-10) {
    stop("Please make sure that the subject-level treatment is constant within subject.")
  }
  if (max(unlist(lapply(split(outcome_variable,f=id_variable),var, na.rm=TRUE)))>1e-10) {
    stop("Please make sure that the subject-level outcome is constant within subject.")
  }
  observed_time_grid <- sort(unique(time_variable));
  wide_id <- sort(unique(id_variable[which(!is.na(id_variable))]));
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
  stopifnot(length(wide_id)>0);
  if (!is.null(covariates_on_outcome_data)) {
    stopifnot(length(id_variable)==nrow(covariates_on_outcome_data));
  }
  if (!is.null(tie_covariates_on_mediator_data)) {
    stopifnot(length(id_variable)==nrow(tie_covariates_on_mediator_data));
  }
  if (!is.null(tve_covariates_on_mediator_data)) {
    stopifnot(length(id_variable)==nrow(tve_covariates_on_mediator_data));
  }
  for (this_id in 1:length(wide_id)) {
    these_rows <- which(id_variable==wide_id[this_id]); 
    if (length(these_rows)>0) {
      if (min(treatment_variable[these_rows], na.rm=TRUE)!=
          max(treatment_variable[these_rows], na.rm=TRUE)) {
        stop("Please make sure the treatment is the same for each observation within subject.")
      }
      wide_treatment[this_id] <- treatment_variable[min(these_rows)];
      if (min(outcome_variable[these_rows], na.rm=TRUE)!=
          max(outcome_variable[these_rows], na.rm=TRUE)) {
        stop("Please make sure the outcome is the same for each observation within subject.")
      }
      wide_outcome[this_id] <- outcome_variable[min(these_rows)];
      if (num_covariates_on_outcome>0) {
        if (max(apply(covariates_on_outcome_data[these_rows,,drop=FALSE],2,var, na.rm = TRUE), na.rm = TRUE)>0) {
          print(paste("Possible problem in covariates_on_outcome_data for participant",wide_id[this_id]));
          print(covariates_on_outcome_data[these_rows,,drop=FALSE]);
          stop("Please make sure that the covariates on the outcome do not vary within subject for this model.")
        }
      }
      if (num_tve_covariates_on_mediator>0) {
        if (max(apply(tve_covariates_on_mediator_data[these_rows,,drop=FALSE],2,var, na.rm = TRUE), na.rm = TRUE)>0) {
          print(paste("Possible problem in tve_covariates_on_mediator_data for participant",wide_id[this_id]));
          print(tve_covariates_on_mediator_data[these_rows,,drop=FALSE]);
          stop("Please make sure that the covariates on the mediator do not vary within subject for this model.")
        }
      }
      if (num_tie_covariates_on_mediator>0) {
        if (max(apply(tie_covariates_on_mediator_data[these_rows,,drop=FALSE],2,var, na.rm = TRUE), na.rm = TRUE)>0) {
          print(paste("Possible problem in tie_covariates_on_mediator_data for participant",wide_id[this_id]));
          print(tie_covariates_on_mediator_data[these_rows,,drop=FALSE]);
          stop("Please make sure that the covariates on the mediator do not vary within subject for this model.")
        }
      }
    }
    stopifnot(length(observed_time_grid)>0);
    for (this_time in 1:length(observed_time_grid)) { 
      this_data <- which(id_variable==wide_id[this_id] &
                           time_variable==observed_time_grid[this_time]); 
      if (length(this_data)==1) {
        wide_mediator[this_id,this_time] <- data[this_data,mediator_variable_name]; 
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
  analyze_data_for_mediation <- function(wide_data, 
                                         indices, 
                                         get_details=FALSE) {
    local_wide_data <- wide_data[indices,];
    #--- Take data frame apart into pieces to use with pfr function 
    wide_mediator <- as.matrix(local_wide_data[,wide_mediator_columns]);
    wide_id <- unlist(local_wide_data[,wide_id_column]);
    wide_treatment <- unlist(local_wide_data[,wide_treatment_column]); 
    wide_outcome <- unlist(local_wide_data[,wide_outcome_column]);
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
    if (num_covariates_on_outcome>0) {
      for (this_one in 1:num_covariates_on_outcome) {  
        assign(paste("wide_",covariates_on_outcome_names[this_one],sep=""),
               wide_covariates_on_outcome[,this_one]);
        new_one <- as.formula(paste("~.+wide_",covariates_on_outcome_names[this_one],sep=""));
        pfr_formula <- update(pfr_formula,new_one);
      }
    }
    ###save.image("working.rdata");
    if (logistic) {
      funreg_MY <- try(pfr(pfr_formula,
                           scale=1,
                           family=binomial()));
      # I wish I could let the user send the data and family in from outside the function,
      # but the pfr function does not allow this due to its unusual implementation
      # as a wrap-around for a hidden call to gam.;
    } else {
      funreg_MY <- try(pfr(pfr_formula,
                           family=gaussian()));
    } 
    if (any(class(funreg_MY)=="try-error")) {
      print(pfr_formula);
      print(funreg_MY);
      stop("Error in running pfr.")
    } 
	alpha_int_estimate <- as.numeric(funreg_MY$coefficient["(Intercept)"]);
    alpha_int_se <- as.numeric(summary(funreg_MY)$se["(Intercept)"]);
    alpha_X_estimate <-  as.numeric(funreg_MY$coefficient["wide_treatment"]);
    alpha_X_se <- as.numeric(summary(funreg_MY)$se["wide_treatment"]);
    temp_coefs <- coef(funreg_MY, coords=list(observed_time_grid));
    alpha_M_estimate <- as.numeric(temp_coefs[,"value"]); 
    alpha_M_se <- as.numeric(temp_coefs[,"se"]);
    time_grid_for_fitting <- observed_time_grid;
    alpha_M_pvalue <- summary(funreg_MY)$s.table[1,"p-value"];
    #--- DIRECT EFFECT OF TREATMENT X ON OUTCOME Y ---;
    glm_formula <- wide_outcome~wide_treatment;
    if (num_covariates_on_outcome>0) {
      for (this_one in 1:num_covariates_on_outcome) { 
        new_one <- as.formula(paste("~.+wide_",covariates_on_outcome_names[this_one],sep=""));
        glm_formula <- update(glm_formula,new_one);
      } 
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
    if (num_tve_covariates_on_mediator>0) {
      for (this_one in 1:num_tve_covariates_on_mediator) {
        new_name <- tve_covariates_on_mediator_names[this_one];
        temp_data_frame <- data.frame(temp=rep(wide_tve_covariates_on_mediator[,this_one],each=nobs));
        colnames(temp_data_frame) <- new_name;
        local_long_data <- cbind(local_long_data, temp_data_frame); 
        new_term <- as.formula(paste("~.+",new_name,sep=""));
        tvem_formula1 <- update(tvem_formula1,new_term);
      }
    }
    if (num_tie_covariates_on_mediator>0) {
      for (this_one in 1:num_tie_covariates_on_mediator) {
        new_name <- tie_covariates_on_mediator_names[this_one];
        temp_data_frame <- data.frame(temp=rep(wide_tie_covariates_on_mediator[,this_one],each=nobs));
        colnames(temp_data_frame) <- new_name;
        local_long_data <- cbind(local_long_data, temp_data_frame); 
      }
    }
    local_long_data <- local_long_data[which(!is.na(local_long_data$mediator)),];
    # listwise deletion to remove empty observations; 
    if (tvem_do_loop) {
      tvem_results_list <- list();
      max_knots <- tvem_num_knots;
      IC_values <- rep(Inf,max_knots+1);
      for (this_num_knots in 0:max_knots) {
        tvem_XM <- suppressWarnings(tvem(data=local_long_data,
                                         formula=tvem_formula1,
                                         time=time,
                                         id=id,
                                         invar_effects=tvem_formula2,
                                         spline_order=tvem_spline_order,
                                         penalty_function_order=tvem_penalty_order,
                                         penalize=tvem_penalize,
                                         num_knots=this_num_knots,
                                         grid=time_grid_for_fitting));
        IC_values[1+this_num_knots] <- ifelse(tvem_use_bic,
                                              tvem_XM$model_information$pseudo_bic,
                                              tvem_XM$model_information$pseudo_aic);
        tvem_results_list[[1+this_num_knots]] <- tvem_XM;
      }
      IC_table <- data.frame(0:max_knots, IC_values);
      colnames(IC_table) <- c("Number_Of_Interior_Knots",ifelse(tvem_use_bic,
                                                                "Pseudo_BIC",
                                                                "Pseudo_AIC"));
      tvem_XM <- tvem_results_list[[which.min(IC_values)]];   
      
    } else {
      if (get_details) {
        tvem_XM <- tvem(data=local_long_data,
                        formula=tvem_formula1,
                        time=time,
                        id=id,
                        invar_effects=tvem_formula2,
                        spline_order=tvem_spline_order,
                        penalty_function_order=tvem_penalty_order,
                        penalize=tvem_penalize,
                        num_knots=tvem_num_knots,
                        grid=time_grid_for_fitting);
      } else {
        tvem_XM <- suppressWarnings(tvem(data=local_long_data,
                                         formula=tvem_formula1,
                                         time=time,
                                         id=id,
                                         invar_effects=tvem_formula2,
                                         spline_order=tvem_spline_order,
                                         penalty_function_order=tvem_penalty_order,
                                         penalize=tvem_penalize,
                                         num_knots=tvem_num_knots,
                                         grid=time_grid_for_fitting));
      }
    }
    gamma_int_estimate <- tvem_XM$grid_fitted_coefficients[[1]]$estimate; 
    gamma_int_se <- tvem_XM$grid_fitted_coefficients[[1]]$standard_error;
    gamma_X_estimate <- tvem_XM$grid_fitted_coefficients[[2]]$estimate; 
    gamma_X_se <- tvem_XM$grid_fitted_coefficients[[2]]$standard_error; 
    # #--- MEDIATED EFFECT OF TREATMENT X THROUGH MEDIATOR M ON OUTCOME Y ---;
    if (length(gamma_X_estimate)!=length(alpha_M_estimate)) {
      stop("Dimension error in functional mediation function;")
    }
    beta_estimate <- mean(gamma_X_estimate*alpha_M_estimate);
    if (get_details) {
      answer_list <- list(time_grid=time_grid_for_fitting,
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
                          direct_effect_details=model_for_direct_effect_XY); 
      if (tvem_do_loop) {
        answer_list$tvem_IC_table <- IC_table; 
      }
      return(answer_list);
    } else { 
	  this_bootstrap <<- this_bootstrap + 1;
	  cat(this_bootstrap);cat(".")
      if (tvem_do_loop) {
        ICs_table_from_bootstraps <<- rbind(ICs_table_from_bootstraps,IC_table[,2]);
      }
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
  ICs_table_from_bootstraps <- NULL;
  cat("Ran original results.\n");
  cat("Working on bootstrap results:\n");
  this_bootstrap <- 0;
  boot1 <- boot(data=wide_data,
                statistic=analyze_data_for_mediation, 
                R=nboot);   
  boot2 <- boot.ci(boot1,conf=1-boot_level,type="norm");
  boot3 <- boot.ci(boot1,conf=1-boot_level,type="basic");
  boot4 <- boot.ci(boot1,conf=1-boot_level,type="perc");
  after_boot <- Sys.time(); 
  bootstrap_results <- list(beta_boot_estimate= norm.ci(boot1,conf=.001)[2],
                            beta_boot_se=sd(boot1$t),
                            beta_boot_norm_lower=boot2$normal[2],
                            beta_boot_norm_upper=boot2$normal[3],
                            beta_boot_basic_lower=boot3$basic[4],
                            beta_boot_basic_upper=boot3$basic[5],
                            beta_boot_perc_lower=boot4$percent[4],
                            beta_boot_perc_upper=boot4$percent[5],
                            boot_level=boot_level,
                            boot1=boot1,
                            time_required=difftime(after_boot,before_boot));
  if (tvem_do_loop) {
    colnames(ICs_table_from_bootstraps) <- paste(0:(ncol(ICs_table_from_bootstraps)-1),"InteriorKnots",sep="");
    bootstrap_results$ICs_table_from_bootstraps <- ICs_table_from_bootstraps;
    bootstrap_results$num_knots_from_bootstraps <- table(paste(apply(ICs_table_from_bootstraps,1,which.min)+1,"knots"));
  }
  answer <- list(observed_time_grid_for_debug=observed_time_grid,
                 wide_data_for_debug=wide_data,
                 original_results=original_results,
                 bootstrap_results=bootstrap_results);
  class(answer) <- "funreg_mediation";
  return(answer);
}  