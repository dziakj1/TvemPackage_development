#' print.tvem:  Print output from a model that was fit by the tvem function.
#' 
#' @param the_tvem The tvem object (output of the tvem function)
#' @param ornate Whether to print lines between different sections of the output for easier reading.
#' 
#' @export
#' @method print tvem




print.tvem <- function(the_tvem, ornate=TRUE) {
  if (ornate) {
    divider <- "======================================================= \n";
  } else {
    divider <- "\n";
  }
  if (ornate) {
    cat(divider);
    cat("Time-Varying Effects Modeling (TVEM) Function Output \n");
    cat(divider);
  }
  cat(paste("Response variable:  ",
            the_tvem$model_information$response_name,
            "\n"));
  cat(paste("Time interval:  ",
            round(min(the_tvem$time_grid),4),
            "to",
            round(max(the_tvem$time_grid),4),"\n"));
  cat("Number of subjects:  ");
  cat(paste(the_tvem$model_information$n_subjects));
  cat("\nEffects specified as time-varying:  ");
  cat(paste(names(the_tvem$grid_fitted_coefficients),sep=" ",collapse=", ")); 
  cat("\nYou can use the plot_tvem function to view their plots.\n");
  if (!is.null(the_tvem$invar_effects_estimates)) {
    cat(divider)
    cat("Effects specified as non-time-varying: \n");
    print(the_tvem$invar_effects_estimates);
  }
  cat(divider);
  cat("Back-end model fitted in mgcv::bam function: \n")
  cat(paste("Method",the_tvem$back_end_model$method))
  cat("\nFormula:\n");
  print(the_tvem$back_end_model$formula); 
  cat(paste("Pseudolikelihood AIC:",round(the_tvem$model_information$pseudo_aic,2)));
  cat(paste("\nPseudolikelihood BIC:",round(the_tvem$model_information$pseudo_bic,2),"\n"));
  if (the_tvem$model_information$used_listwise_deletion) {
    cat("Note: Used listwise deletion for missing data.\n");
  }
  cat(divider);
}
