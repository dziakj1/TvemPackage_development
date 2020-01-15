library(mgcv);
tvem <- function(data,
                 formula,
                 id,
                 time,
                 invar_effects=NULL,
                 family=gaussian(), # use binomial() for binary;
                 num_knots=20,  
                 # number of interior knots per function (not counting exterior knots);
                 spline_order=3,
                 # shape of function between knots, with default of 3 for cubic;
                 penalty_function_order=2,
                 smoothing_penalty_multiplier=1,
                 grid=100,
                 alpha=.05
) {   
  ##################################
  # Process the input;
  ################################## 
  m <- match.call(expand.dots = FALSE);
  m$formula <- NULL;
  m$invar_effects <- NULL;
  m$family <- NULL;
  m$grid <- NULL;
  if (is.matrix(eval.parent(m$data))) {
    m$data <- as.data.frame(data);
  }
  m[[1]] <- quote(stats::model.frame);
  m <- eval.parent(m); 
  id_variable_name <- as.character(substitute(id));
  id_variable <- m[,id_variable_name];
  time_variable_name <- as.character(substitute(time));
  time_variable <- m[,time_variable_name]; 
  if (is.character(family)) {
    # Handle the possibility that the user specified family
    # as a string instead of an object of class family.
    family <- tolower(family);
    if (family=="normal" | 
        family=="gaussian" | 
        family=="linear" | 
        family=="numerical") {
      family <- gaussian();
    }
    if (family=="binary" | 
        family=="bernoulli" | 
        family=="logistic" | 
        family=="binomial") {
      family <- binomial();
    }
    if (family=="count" | 
        family=="poisson"|
        family=="loglinear") {
      family <- poisson();
    }
  }
  the_terms <- terms(formula);
  whether_intercept <- attr(the_terms,"intercept");
  if (whether_intercept != 1) {
    stop("An intercept function is required in the current version of this function")
  }
  formula_variable_names <- attr(attr(the_terms,"factors"),"dimnames")[[1]];
  response_name <- formula_variable_names[[1]];
  num_varying_effects <- length(formula_variable_names)-1; # not including intercept;
  if (num_varying_effects>0) {
    varying_effects_names <- c(formula_variable_names[2:(1+num_varying_effects)]);
  }
  if (length(num_knots)>1) {
    stop("Please specify a single number for num_knots.")
  };
  crit_value <- qnorm(1-alpha/2);
  ##################################################
  # Construct regular grid for plotting fitted coefficient functions;
  ##################################################
  if (is.null(grid)) {
    grid <- 100;
  }
  if (length(grid)==1) {
    grid_size <- as.integer(grid);
    min_time <- min(time_variable, na.rm = TRUE);
    max_time <- max(time_variable, na.rm = TRUE);
    time_grid <- seq(min_time, max_time, length=grid_size);
  }
  if (length(grid)>1) {
    time_grid <- grid;
  }
  ##################################
  # Construct formula to send in to the back-end computation function.
  # This function is bam (big generalized additive models) in the  
  # mgcv (Mixed GAM Computation Vehicle) package by Simon Wood of R Project.
  ##################################
  if (is.null (invar_effects)) {
    invar_effects_names <- NULL;
  } else {
    invar_effects_names <- attr(terms(invar_effects),"term.labels");
  }
  num_invar_effects <- length(invar_effects_names);
  bam_formula <- as.formula(paste(response_name," ~ ", "1"));
  
  if (length(intersect(invar_effects_names,
                       varying_effects_names))>0) {
    stop(paste("Please do not specify the same variable",
               "as having time-varying and time-invariant effects."));
  }
  if (num_varying_effects>0) {
    for (i in 1:num_varying_effects) {
      new_text <- paste("~ . +",varying_effects_names[i]); 
      bam_formula <- update(bam_formula,as.formula(new_text));
      # add linear effect of each time-varying-effect to the formula
    }
  }
  if (num_invar_effects>0) {
    for (i in 1:num_invar_effects) {
      new_text <- paste("~ . +",invar_effects_names[i]);
      bam_formula <- update(bam_formula,as.formula(new_text));
      # add effect of each non-time-varying-effect to the formula
    }
  }
  new_text <- paste("~ . + s(",time_variable_name,",bs='ps',by=NA,pc=0,",
                    "m=c(",
                    spline_order-1,
                    ",",
                    penalty_function_order,
                    "),",
                    "k=",
                    num_knots+spline_order+1,
                    ")",sep=""); 
  # for time-varying intercept; 
  bam_formula <- update(bam_formula,as.formula(new_text)); 
  # for time-varying intercept;
  if (num_varying_effects>0) {
    for (i in 1:num_varying_effects) {
      this_covariate_name <- varying_effects_names[i];
      new_text <- paste("~ . + s(",time_variable_name,",bs='ps',by=",
                        this_covariate_name,
                        ", pc=0,",
                        "m=c(",
                        spline_order-1,
                        ",",
                        penalty_function_order,"),",
                        "k=",
                        num_knots+spline_order+1,")",
                        sep=""); 
      bam_formula <- update(bam_formula,as.formula(new_text));
    }
  } 
  data_for_analysis <- as.data.frame(eval.parent(data)); 
  model1 <- bam(bam_formula,
                data=data_for_analysis,
                gamma=smoothing_penalty_multiplier);
  ##################################
  # Extract coefficient estimates;
  ##################################
  grand_intercept <- model1$coefficients["(Intercept)"];
  bam_coefs <- predict.bam(model1,type="terms"); 
  estimated_b0 <- grand_intercept + 
    bam_coefs[,paste("s(",time_variable_name,")",sep="")];
  estimated_b <- NULL;
  if (num_varying_effects>0) {
    for (i in 1:num_varying_effects) {
      this_covariate_name <- varying_effects_names[i];
      this_covariate_values <- model1$model[,this_covariate_name];
      this_spline_term <- paste("s(",time_variable_name,"):",
                                this_covariate_name,sep=""); 
      this_covariate_function_values <-  (bam_coefs[,this_covariate_name]+
                                            bam_coefs[,this_spline_term])/this_covariate_values;
      estimated_b <- cbind(estimated_b,
                           this_covariate_function_values);
    }
    colnames(estimated_b) <- varying_effects_names;
  }
  ##################################
  # Calculate confidence intervals using sandwich formula
  # to better take into account the existence of within-subject
  # correlation;
  ##################################
  # Consult design matrix of spline bases on grid of observed times:
  design <- predict.bam(model1,type="lpmatrix");
  # Construct design matrix of spline bases on regular grid of times:
  temp_data <- data.frame(as.matrix(0*model1$model[1,])%x%rep(1,length(time_grid)));
  colnames(temp_data) <- colnames(as.data.frame(model1$model));
  temp_data[,time_variable_name] <- time_grid;
  grid_design <- predict.bam(model1,type="lpmatrix",newdata = temp_data); 
  # Working independence variance estimates (will make sandwich):
  bread <- model1$Vc;
  npar <- length(model1$coefficients);
  subject_id <- id_variable; 
  meat <- matrix(0,npar,npar);  # apologies to vegetarians -- 
  # peanut butter is also fine!
  for (i in unique(subject_id)) {
    these <- which(subject_id==i);
    if (length(these)>0) {
      meat <- meat + tcrossprod(crossprod(
        design[these,],model1$residuals[these]));
    }
  } 
  sandwich <- bread %*% meat %*% bread;  
  general_term_name <- sub(" .*","", 
                           gsub(x=names(model1$coefficients),
                                pattern="[.]",
                                replacement=" "));
  # removes the .1, .2, .3, etc., from after the term name in the names
  # of coefficients.
  indices_for_b0 <- which((general_term_name=="(Intercept)") | 
                            (general_term_name == paste("s(",
                                                        time_variable_name,
                                                        ")",
                                                        sep="")));
  basis_for_b0 <-  design[,indices_for_b0];
  estimated_b0_hard_way <-  basis_for_b0%*%model1$coefficients[indices_for_b0];  
  stopifnot(max(abs(estimated_b0_hard_way-estimated_b0), na.rm =TRUE)<1e-10); #just for double checking;
  grid_basis_for_b0 <- grid_design[,indices_for_b0]; 
  grid_estimated_b0 <- grid_basis_for_b0 %*% model1$coefficients[indices_for_b0];
  cov_mat_for_b0 <- sandwich[indices_for_b0,indices_for_b0];
  # Fitted b0 on observed times:
  temp_function <- function(i){return(tcrossprod(basis_for_b0[i,] %*%
                                                   cov_mat_for_b0, 
                                                 basis_for_b0[i,]))};
  standard_error_b0 <- drop(sqrt(sapply(X=1:nrow(basis_for_b0),
                                        FUN=temp_function)));
  upper_b0 <- estimated_b0 + crit_value*standard_error_b0;
  lower_b0 <- estimated_b0 - crit_value*standard_error_b0;
  fitted_coefficients <- list("(Intercept)"=data.frame(estimate = estimated_b0,
                                                       standard_error = standard_error_b0,
                                                       upper = upper_b0,
                                                       lower = lower_b0));
  # Fitted b0 on regular grid:
  temp_function <- function(i){return(tcrossprod(grid_basis_for_b0[i,] %*%
                                                   cov_mat_for_b0, 
                                                 grid_basis_for_b0[i,]))};
  grid_standard_error_b0 <- drop(sqrt(sapply(X=1:nrow(grid_basis_for_b0),
                                             FUN=temp_function)));
  grid_upper_b0 <- grid_estimated_b0 + crit_value*grid_standard_error_b0;
  grid_lower_b0 <- grid_estimated_b0 - crit_value*grid_standard_error_b0;
  grid_fitted_coefficients <- list("(Intercept)"=data.frame(estimate = grid_estimated_b0,
                                                            standard_error = grid_standard_error_b0,
                                                            upper = grid_upper_b0,
                                                            lower = grid_lower_b0));
  if (num_varying_effects>0) {
    for (i in 1:num_varying_effects) {
      this_covariate_name <- varying_effects_names[i];
      this_covariate_values <- model1$model[,this_covariate_name];
      this_spline_term <- paste("s(",time_variable_name,"):",
                                this_covariate_name,sep="");
      
      indices_for_this_b <- which((general_term_name == this_covariate_name) | 
                                    (general_term_name == this_spline_term)); 
      estimated_this_b <- (basis_for_b0%*%model1$coefficients[indices_for_this_b]);
      grid_estimated_this_b <- (grid_basis_for_b0%*%model1$coefficients[indices_for_this_b]);
      # Using the b0 basis for b1 makes code simpler but assumes that number of knots, 
      # spline order, and coefficient order are the same between 
      # intercept and substantive effect; perhaps remove this assumption
      # in some later version;   
      estimated_this_b <- estimated_b[,i]; 
      estimated_this_b_hard_way <- as.vector(basis_for_b0%*%model1$coefficients[indices_for_this_b]);
      stopifnot(max(abs(estimated_this_b_hard_way-estimated_b[,i]), na.rm=TRUE)<1e-10); #just for double checking;
      cov_mat_for_this_b <- sandwich[indices_for_this_b,indices_for_this_b];
      temp_function <- function(i){return(tcrossprod(basis_for_b0[i,] %*%
                                                       cov_mat_for_this_b, 
                                                     basis_for_b0[i,]))};
      standard_error_this_b <- drop(sqrt(sapply(X=1:nrow(basis_for_b0),
                                                FUN=temp_function)));
      upper_this_b <- estimated_this_b + crit_value*standard_error_this_b;
      lower_this_b <- estimated_this_b - crit_value*standard_error_this_b;
      temp_list <- list(data.frame(estimate = estimated_this_b,
                                   standard_error = standard_error_this_b,
                                   upper = upper_this_b,
                                   lower = lower_this_b));
      names(temp_list) <- this_covariate_name;
      fitted_coefficients <- c(fitted_coefficients,
                               temp_list);
      temp_function <- function(i){return(tcrossprod(grid_basis_for_b0[i,] %*%
                                                       cov_mat_for_this_b, 
                                                     grid_basis_for_b0[i,]))};
      grid_standard_error_this_b <- drop(sqrt(sapply(X=1:nrow(grid_basis_for_b0),
                                                     FUN=temp_function)));
      grid_upper_this_b <- grid_estimated_this_b + crit_value*grid_standard_error_this_b;
      grid_lower_this_b <- grid_estimated_this_b - crit_value*grid_standard_error_this_b;
      temp_list_grid <- list(data.frame(estimate = grid_estimated_this_b,
                                        standard_error = grid_standard_error_this_b,
                                        upper = grid_upper_this_b,
                                        lower = grid_lower_this_b));
      names(temp_list_grid) <- this_covariate_name;
      grid_fitted_coefficients <- c(grid_fitted_coefficients,
                                    temp_list_grid);
    }
    names(fitted_coefficients[2:length(fitted_coefficients)]) <- varying_effects_names;
    names(grid_fitted_coefficients[2:length(grid_fitted_coefficients)]) <- varying_effects_names;
  }
  invar_effects_estimates <- NULL;
  if (num_invar_effects>0) {
    invar_effects_indices <- which(general_term_name %in% invar_effects_names);
    invar_effects_est <- coef(model1)[invar_effects_names];
    invar_effects_se <- sqrt(diag(sandwich[invar_effects_indices,
                                           invar_effects_indices,
                                           drop=FALSE]));
    invar_effects_estimates <- data.frame(estimate=invar_effects_est,
                                          standard_error=invar_effects_se);
  }
  ##################################
  # Return answers;
  ##################################
  model_information <- list(whether_intercept=whether_intercept,
                            response_name=response_name,
                            num_varying_effects=num_varying_effects,
                            varying_effects_names=varying_effects_names,
                            num_invar_effects=num_invar_effects,
                            invar_effects_names=invar_effects_names);
  answer <- list(time_grid=time_grid,
                 grid_fitted_coefficients=grid_fitted_coefficients,
                 invar_effects_estimates=invar_effects_estimates,
                 model_information=model_information,
                 back_end_model=model1);
  class(answer) <- "tvem";
  return(answer);
}  