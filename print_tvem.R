#' print.tvem:  Print output from a model that was fit by the tvem function.
#' 
#' @param the_tvem The tvem object (output of the tvem function)
#' @param ornate Whether to print lines between different sections of the output for easier reading.
#' 
#' @export
#' @S3method print tvem




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
            min(the_tvem$time_grid),
            "to",
            max(the_tvem$time_grid)),"\n");
  cat("Effects specified as time-varying: \n       ");
  cat(paste(names(the_tvem$grid_fitted_coefficients),sep=" ",collapse=", ")); 
  cat("\n (use plot_tvem function to view their plots)");
  if (!is.null(the_tvem$invar_effects_estimates)) {
    cat("\nEffects specified as non-time-varying: \n");
    print(the_tvem$invar_effects_estimates);
  }
  cat(divider);
  cat("Back-end model fitted in mgcv::bam function: \n")
  print(summary(the_tvem$back_end_model));
  cat("Note:  
      Tests of 'Parametric coefficients' for TVEM coefficient functions here
      can be interpreted as testing whether the relationship is significant at time zero.  
      Tests of 'Approximate significance of smooth terms' refer to whether it differs 
      at other times relative to time zero. \n"); 
  cat(divider);
}
