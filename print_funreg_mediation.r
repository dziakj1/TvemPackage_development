print_funreg_mediation <- function(fitted_funmed_model) {
  cat("======================================================= \n");
  cat("Functional Regression Mediation Function Output \n");
  cat("======================================================= \n");
  cat("Indirect effect bootstrap estimate:\n");
  cat(fitted_funmed_model$bootstrap_results$beta.boot.estimate);
  cat("\nIndirect effect bootstrap confidence interval:");
  cat("\n... by normal method:\n")
  cat(c(round(answers1$bootstrap_results$beta.boot.norm.lower,4), ", ",
        round(answers1$bootstrap_results$beta.boot.norm.upper,4)));
  cat("\n... by percentile method:\n")
  cat(c(round(answers1$bootstrap_results$beta.boot.perc.lower,4), ", ",
          round(answers1$bootstrap_results$beta.boot.perc.upper,4))); 
  cat("\nComputation time:\n")
  print(answers1$bootstrap_results$time.required);
  cat("======================================================= \n");
  cat("TVEM Model for Predicting Mediator from Treatment:\n");
  print_tvem(fitted_funmed_model$original_results$tvem_XM_details, 
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
}