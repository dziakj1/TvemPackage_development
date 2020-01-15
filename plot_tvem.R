plot_tvem <- function(the_tvem,
                      use_panes=TRUE,
                      which_plot=NULL,
                      diagnostics=FALSE) {
  if (diagnostics) {
    if (use_panes) {par(mfrow=c(1,2))};
    hist(the_tvem$back_end_model$fitted,
         main="Residuals",
         xlab="Residual");
    plot(the_tvem$back_end_model$fitted,
         the_tvem$back_end_model$residuals,
         main="Fitted versus residuals",
         xlab="Fitted",
         ylab="Residuals");
  } else {
    num_tv_coefs <- length(the_tvem$grid_fitted_coefficients);
    if (use_panes) {
      if (num_tv_coefs==1) {par(mfrow=c(1,1));}
      if (num_tv_coefs==2) {par(mfrow=c(1,2));}
      if (num_tv_coefs==3 | num_tv_coefs==4) {par(mfrow=c(2,2));}
      if (num_tv_coefs==5 | num_tv_coefs==6) {par(mfrow=c(3,2));}
      if (num_tv_coefs==7 | num_tv_coefs==8) {par(mfrow=c(4,2));}
      if (num_tv_coefs==9) {par(mfrow=c(3,3));}
      if(num_tv_coefs>10) {stop("Too many functions to plot in panes.");}  
    }
    the_grid <- the_tvem$time_grid;
    temp_plot_function <- function(the_var_name,
                                   the_grid,
                                   the_coef) {
      plot(x=the_grid,
           y=the_coef$estimate,
           main=paste("TVEM coefficient:\n",the_var_name),
           xlab=expression(t),
           ylab=expression(beta(t)),
           ylim=c(ymin,ymax),
           type="l",
           lty="solid",
           lwd=2, 
           mgp=c(2,1,0));
      abline(h=0,lty="dotted");
      lines(x=the_grid,
            y=the_coef$lower);
      lines(x=the_grid,
            y=the_coef$upper);
      return(0);
    }
    if (is.null(which_plot)) {
      for (which_plot in 1:num_tv_coefs) {
        the_coef <- the_tvem$grid_fitted_coefficients[[which_plot]];
        the_var_name <- names(the_tvem$grid_fitted_coefficients)[which_plot];
        ymin <- min(0,min(the_coef$lower));
        ymax <- max(0,max(the_coef$upper));
        temp_plot_function(the_var_name,
                           the_grid,
                           the_coef);
      }
    } else {
      the_coef <- the_tvem$grid_fitted_coefficients[[which_plot]];
      the_var_name <- names(the_tvem$grid_fitted_coefficients)[which_plot];
      ymin <- min(0,min(the_coef$lower));
      ymax <- max(0,max(the_coef$upper));
      temp_plot_function(the_var_name,
                         the_grid,
                         the_coef);
    }
  }
}