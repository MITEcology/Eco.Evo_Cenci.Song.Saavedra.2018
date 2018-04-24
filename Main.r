rm(list = ls())
suppressMessages(require(Rcpp))
suppressMessages(require(RcppArmadillo))
source('Function.r')
sourceCpp('Cfun.cpp')
require(parallel)
##########################################################################################
################################################ Main file 
#### Choose the number of realizations
ensamble = 1:1000
#### Choose the perturbation
choice = c('Matrix_random_perturbation', 'Growth_random_perturbation')
title_list = c('Random perturbation of the interactions', 'Random perturbation of the growth rates')
### The methods are:
### ODE Solve L-V explicitely (slow)
### Inverse solve the stationary state
### NOTA IMPORTANTE if you use method 'ODE' than you need to adjust the initial conditions of 'Integra' in the main function
### Use 'inverse' for clarity
method = c('ODE', 'inverse')
#####################################
#####################################
max_step = 1000; 
#### make it a multiple of 3 just to make a complete modular network
network_size = 21
#### Choose the number of cores to use for parallel computing
Lavoratori = detectCores() - 1
cl <- makeCluster(Lavoratori, type = "FORK")
cat('###### ', network_size, 'species with', length(ensamble), 'realizations\n')
ForFigure = list()
for(chc in 1:length(choice)){
  number_of_not_survival = unlist(parLapply(cl, ensamble, main_function, choice[chc], max_step, method[2], network_size))
  number_of_not_survival = number_of_not_survival[number_of_not_survival!= 3]
  ForFigure[[chc]] = 1 - number_of_not_survival/length(ensamble)
}
stopCluster(cl)
##########################################################################################################
### Plot the result. This is the histogram of the number fo the realization without an extinction
plt_ = TRUE
save_figure = TRUE
if(plt_ == TRUE){
  PltCombined(ForFigure, title_list, save_figure)
}

