plot_funreg_mediation <- function(fitted_funmed_model,
                                  use_panes=TRUE, 
                                  what_plot=c("pfr",
                                              "coefs",
                                              "tvem")) {
  if (use_panes==FALSE) {par(mfrow=c(1,1));}
  what <- match.arg(what_plot)
  if (what=="pfr") {
    par(mfrow=c(1,1));
    plot(fitted_funmed_model$original_results$funreg_MY_details,
         xlab="Time",
         ylab="Coefficient function",
         main="Functional effect of mediator on outcome");
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
    est2 <- round(fitted_funmed_model$bootstrap_results$beta.boot.estimate,2); 
    text(x=1.2,y=1,pos=4,labels=est2);
  }
  if (what=="tvem") {
    plot_tvem(fitted_funmed_model$original_results$tvem_XM_details);
  }
}