#' print.funreg_mediation:  Print output from a model that was fit by the funreg_mediation function.
#' 
#' @param fitted_funmed_model The funreg_mediation object (output of the funreg_mediation function)
#' 
#' @export
#' @method print funreg_mediation

print.funreg_mediation <- function(fitted_funmed_model) {
  cat("======================================================= \n");
  cat("Functional Regression Mediation Function Output \n");
  cat("======================================================= \n");
  cat("Indirect effect bootstrap estimate:\n");
  cat(fitted_funmed_model$bootstrap_results$beta_boot_estimate);
  cat("\nIndirect effect bootstrap confidence interval:");
  cat("\n... by normal method:\n");
  cat(c(round(fitted_funmed_model$bootstrap_results$beta_boot_norm_lower,4), ", ",
        round(fitted_funmed_model$bootstrap_results$beta_boot_norm_upper,4)));
  cat("\n... by percentile method:\n");
  cat(c(round(fitted_funmed_model$bootstrap_results$beta_boot_perc_lower,4), ", ",
          round(fitted_funmed_model$bootstrap_results$beta_boot_perc_upper,4))); 
  cat("\nComputation time:\n");
  print(fitted_funmed_model$bootstrap_results$time_required);
  cat("======================================================= \n");
  cat("TVEM Model for Predicting Mediator from Treatment:\n");
  print.tvem(fitted_funmed_model$original_results$tvem_XM_details, 
             ornate=FALSE);
  cat("======================================================= \n");
  cat("Functional Mediation Model for Predicting Mediator from Treatment: \n");
  print(fitted_funmed_model$original_results$funreg_MY_details);
  cat("Scalar terms:\n")
  temp <- fitted_funmed_model$original_results$funreg_MY_details$coefficients;
  print(temp[-grep(x=names(temp),pattern=":")]);
  cat("======================================================= \n");
  cat("Parametric model for Predicting Outcome from Treatment and Mediator: \n");
  print(summary(fitted_funmed_model$original_results$direct_effect_details));
  cat("======================================================= \n");
  if(!is.null(fitted_funmed_model$original_results$tvem_IC_table)) {
    cat("ICs table for selecting number of interior knots in TVEM:\n");
    print(fitted_funmed_model$original_results$tvem_IC_table);
    cat("======================================================= \n");
  }
  if(!is.null(fitted_funmed_model$bootstrap_results$num_knots_from_bootstraps)) {
    cat("Numbers of interior knots selected in bootstrap samples:\n");
    print(fitted_funmed_model$bootstrap_results$num_knots_from_bootstraps);
    cat("======================================================= \n");
  }  
}