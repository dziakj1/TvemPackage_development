select_tvem <- function(maxKnots=5,
                        keepGoingIfTooFew=TRUE,
                        use_bic=FALSE,
                        penalize=FALSE,
                        ...) { 
  #' select_tvem:  Select number of knots for an unpenalized TVEM.
  #' @param maxKnots  The maximum number of knots to try (0 through maxKnots)
  #' @param keepGoingIfTooFew  Whether to continue in a stepwise fashion if the maxKnots 
  #' does not seem to be high enough
  #' @param use_bic Whether to use BIC (FALSE) instead of AIC (TRUE) when selecting the best 
  #' model.  Note that both of these IC's are calculated from the working-independence 
  #' pseudolikelihood rather than the unknown true likelihood.  However, for BIC, the
  #' sample size is taken to be the number of subjects, not the number of observations.
  #' @param penalize Whether to also include a penalty function in estimation
  #' @param  ... Other inputs to be sent along to each call to the tvem function. 
  #' 
  #' @return A list with two components:
  #' \describe{
  #'   \item{ICsTable}{The different numbers of knots which were tried.}
  #'   \item{bestFit}{The tvem object containing the results of the fitted model.}
  #'   }
  args1 <- match.call();  
  num_knots_values <- 0:maxKnots;  # will try at least each of these values for num_knots;
  IC_values <- rep(NA,length(num_knots_values)); 
  done_looking <- FALSE;
  num_knots_values_index <- 0;
  while (done_looking==FALSE) {
    num_knots_values_index <- num_knots_values_index + 1;
    this_num_knots <- num_knots_values[num_knots_values_index];
    print(paste("Trying with",this_num_knots,"knots."));
    more_args <- as.list(args1)[-1];
    more_args$num_knots <- this_num_knots;
    more_args$use_naive_se <- TRUE;
    more_args$print_gam_formula <- FALSE; 
    ans1 <- suppressWarnings(do.call(tvem, more_args));  
    if (use_bic==FALSE) {
      # use AIC;
      IC_values[num_knots_values_index] <- ans1$back_end_model$deviance + 
        2*sum(ans1$back_end_model$edf);
    }
    if (use_bic==TRUE) {
      # use BIC;
      IC_values[num_knots_values_index] <- ans1$back_end_model$deviance + 
        nsub*sum(ans1$back_end_model$edf);
    }
    if (num_knots_values_index==length(num_knots_values)) {
      # If you have come to end of list of values to try for number of knots
      if ((which.min(IC_values)==num_knots_values_index) &
          (keepGoingIfTooFew==TRUE)) {
        # If there still might be more knots needed then try more
        num_knots_values <- c(num_knots_values, max(num_knots_values)+1);
        IC_values <- c(IC_values, NA);
      }  else {
        done_looking <- TRUE;
      }
    }
  } 
  best_num_knots <- num_knots_values[which.min(IC_values)];  
  more_args <- as.list(args1)[-1];
  more_args$num_knots <- best_num_knots;
  more_args$use_naive_se <- FALSE; 
  ans1 <- do.call(tvem, more_args);
  ICsTable <- cbind(knots=num_knots_values, ic=IC_values);
  print(ICsTable);
  return(list(ICsTable=ICsTable,
         bestFit=ans1));
}