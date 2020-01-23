library(mgcv);
library(refund);
library(boot);
 
source("simulate_tvem_example.R");
long_simulated_data_1 <- simulate_tvem_example();

source("tvem.r");
ans1 <- tvem(data=long_simulated_data_1,
            formula = y~tvx,
            id=subject_id,
            time=times,
            invar_effects = ~z1+z2);

source("print_tvem.r");
print.tvem(ans1);

source("plot_tvem.r");
plot.tvem(ans1);	

source("simulate_functional_mediation_example.R"); 
long_simulated_data_2 <- simulate_functional_mediation_example(); 

source("funreg_mediation.r");
ans2 <- funreg_mediation(data=long_simulated_data_2,
                             treatment=X,
                             mediator=M,
                             outcome=Y,
                             id=subject_id,
                             time=time_value,
                             tve_covariates_on_mediator=~Z1,
                             tie_covariates_on_mediator=~Z2,
                             covariates_on_outcome=~Z1+Z2,
                             nboot=4); 

source("print_funreg_mediation.r");
print.funreg_mediation(ans2);	

source("plot_funreg_mediation.r");
plot.funreg_mediation(ans2,what_plot="pfr");
plot.funreg_mediation(ans2,what_plot="coefs");
plot.funreg_mediation(ans2,what_plot="tvem");
