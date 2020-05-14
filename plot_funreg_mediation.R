#' plot.funreg_mediation:  Draw plots for a funreg_mediation model.
#' 
#' Produces plots from a funreg_mediation object produced by 
#' the funreg_mediation function.  These plots will be shown on the default
#' output device (likely the screen);  they can of course be 
#' written to a file instead, by preceding the call to plot.tvem 
#' with a call to png(), pdf(), or other R graphic file output functions.
#' 
#' @param fitted_funmed_model The funreg_mediation object to be plotted.
#' @param use_panes Whether to plot multiple coefficient functions in a single image.
#' @param what_plot One of "pfr","coefs", or "tvem."  These options are as follows:
#' \describe{
#' \item{pfr}{For a "pfr" plot, the functional coefficient for predicting the outcome from the mediator (conditional on treatment), is shown, by calling the plot method for the penalized functional regression results.  See the documentation for plot.gam() in the mgcv function for more information, because the pfr() function uses gam() as a back end.}
#' \item{coefs}{For a "coefs" plot, the three important functional coefficients in the model (intercept for predicting mediator, treatment effect on mediator, and mediator effect on outcome adjusting for treatment) are plotted one after another.  That is, the plots are shown for gamma_int_estimate, gamma_X_estimate, and alpha_M_estimate, each as a function of time_grid.}
#' \item{tvem}{For a "tvem" plot, the functional coefficients in the TVEM model predicting mediator from treatment are displayed.}
#' }
#' 
#' @export
#' @method plot funreg_mediation
#' 
plot.funreg_mediation <- function(fitted_funmed_model,
                                  use_panes=TRUE, 
                                  what_plot=c("pfr",
                                              "coefs",
                                              "tvem")) {
  if (use_panes==FALSE) {par(mfrow=c(1,1));}
  what <- match.arg(what_plot)
  if (what=="pfr") {
    if (use_panes) {par(mfrow=c(1,1));}
    plot(fitted_funmed_model$original_results$funreg_MY_details,
         xlab="Time",
         ylab="Coefficient function",
         main="Functional effect of mediator on outcome"); # plot.gam function;
  } 
  if (what=="coefs") {
    if (use_panes) {par(mfrow=c(2,2));}
    plot(fitted_funmed_model$original_results$time_grid,
         fitted_funmed_model$original_results$gamma_int_estimate,
         xlab=expression(t),
         ylab=expression(gamma[int](t)),
         main="Intercept for Mediator");
    plot(fitted_funmed_model$original_results$time_grid,
         fitted_funmed_model$original_results$gamma_X_estimate,
         xlab=expression(t),
         ylab=expression(gamma[X](t)),
         main="Treatment on Mediator");
    plot(fitted_funmed_model$original_results$time_grid,
         fitted_funmed_model$original_results$alpha_M_estimate,
         xlab=expression(t),
         ylab=expression(alpha[M](t)),
         main="Mediator on Outcome");  
    plot(x=0:5,y=0:5,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="");
    text(x=1.1,y=4,pos=4,labels=paste("Mediator intercept:"));
    est1 <- round(fitted_funmed_model$original_results$alpha_int_estimate,2); 
    text(x=1.2,y=3,pos=4,est1);
    text(x=1.1,y=2,pos=4,labels=paste("Indirect effect:"));
    est2 <- round(fitted_funmed_model$bootstrap_results$beta_boot_estimate,2); 
    text(x=1.2,y=1,pos=4,labels=est2);
  }
  if (what=="tvem") {
    plot.tvem(fitted_funmed_model$original_results$tvem_XM_details);
  }
}