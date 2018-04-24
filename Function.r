library(deSolve)

#### Solve LV
#### The parameter exponent allow to perturb the form of the model itself
#### The two methods are:  
#### 1) 'ODE': integrate the dynamical system
#### 2) 'inverse: inverte the interaction matrix at stationarity
Integra <- function(alpha, r, method, N0 = NULL, dt = 0.1, MaxTime = 1000){
  if(method == 'ODE'){
    parms <- list(r = r, alpha = alpha)
    time_step <- seq(0, MaxTime, by = dt)
    model <- function(t,N,parms){ dN <- N * (parms$r - parms$alpha %*% N); list(dN)}
    sol <- ode(N0,time_step,model,parms)
  }
    else{
      sol = FixedPoint(alpha, r)
    }
    return(sol)
}

Check_ForCoexistence <- function(X, method){
  ### Return 0 if one species get extinct
  ### Return 1 if not
  if(method == 'ODE'){
    X = X[nrow(X),]
  }else{
    X = X
  }
    X = X[X < .Machine$double.eps]
  if(length(X) > 0){
    ### Extinction
    return(0)
  } else{
    ### Coexistence
    return(1)
  }
}

######### Make Random Stable Network
Make_X <- function(rown = 6, coln = 9, ct = 0.255, mu = 0.1){
  ### the number of column is going to be the linear dimension of your final matrix
  ### Make the binary matrix
  esci = 0
  while(esci == 0){
    X <- matrix(1, rown, coln)
    C = sum(X)/length(X)
    while(C > ct){
      X[sample(rown, 1), sample(coln, 1)] = 0
      C = sum(X)/length(X)
    }
    if (prod(rowSums(X))*prod(colSums(X)) == 0){
    }else{esci = 1}}
  X = t(X) %*% X
  X = lapply(1:nrow(X), function(i, X) return(X[,i]/sum(X[,i])), X)
  X <- matrix(unlist(X), ncol = length(X), byrow = TRUE)
  ### Set the diagonal
  X = X - diag(diag(X)) + diag(nrow(X))
  ### Normalize
  D <- diag(1/sqrt(diag(t(X) %*% X)))
  X <- X %*% D
  return(X)
}
################## Perturbation in the interactions
#### Perturbations
GrowthRandomPerturbation <- function(x, y){
  return(abs(x + y))
}
GrowthDirectedPerturbation <- function(x, i, rd){
  x[i] = x[i] + rd
  return(abs(x))
}

main_function <- function(i, ch, mx_stp, metodo, size){
  RN = Make_X(coln = size)
  controller = mu = 0
  while(controller == 0){ mu = mu + 0.05;  MN = ModularNetwork(mu, size); if(mean(MN) > mean(RN)){ controller = 1} else if(mu > 4){ RN = Make_X(coln = size); mu = 0}}
  controller = mu = 0
  while(controller == 0){ mu = mu + 0.05;  NN = NestedNetwork(mu, size); if(mean(NN) > mean(RN)){ controller = 1} else if(mu > 4){ RN = Make_X(coln = size); mu = 0}}
  r_MN = Centroid(MN)
  r_NN = Centroid(NN)
  step = 0
  ck = 1
  #### In case everyone dies or two dies and one doesn't or none of them die then ignore it by returning three
  ToReturn = 3
  if(ch == 'Matrix_random_perturbation'){
    ### Take the binary structure of the interaction matrix
    ### So that you only perturb the existing interactions
    Aux_MN = ifelse(MN > 0, 1, 0)
    Aux_NN = ifelse(NN > 0, 1, 0)
    while(ck == 1 & step < mx_stp){
      pert = replicate(nrow(MN), rnorm(nrow(RN), 0, 0.1*(sd(MN) + sd(NN))/2))
      MN = abs(MN + Aux_MN * pert)
      NN = abs(NN + Aux_NN * pert)
      ck_MN = Check_ForCoexistence(Integra(MN, r_MN, metodo), metodo)
      ck_NN = Check_ForCoexistence(Integra(NN, r_NN, metodo), metodo)
      if(ck_MN == 1 & ck_NN == 1){ 
        ### continue or return 3
        ck = 1; step = step + 1
      } else if(ck_MN == 0 & ck_NN == 1){ 
        ### return -1, i.e. failure in the modular network
        ck = 0; ToReturn = -1
      } else if (ck_MN == 1 & ck_NN == 0){
        ### return 0, i.e. failure in Nested Network
        ck = 0; ToReturn = 1
      }  else {
        ck = 0; 
      } 
    }
  }else if(ch == 'Growth_random_perturbation'){
    while(ck == 1 & step < mx_stp){
      pert = rnorm(length(r_MN), 0, 0.05*(mean(r_MN) + mean(r_NN))/2)
      r_MN = abs(GrowthRandomPerturbation(r_MN, pert))
      r_NN = abs(GrowthRandomPerturbation(r_NN, pert))
      ck_MN = Check_ForCoexistence(Integra(MN, r_MN, metodo), metodo)
      ck_NN = Check_ForCoexistence(Integra(NN, r_NN, metodo), metodo)
      step = step + 1
      if(ck_MN == 1 & ck_NN == 1){ 
        ### continue or return 3
        ck = 1; 
      } else if(ck_MN == 0 & ck_NN == 1){ 
        ### return -1, i.e. failure in the modular network
        ck = 0; ToReturn = -1
      } else if (ck_MN == 1 & ck_NN == 0){
        ck = 0; ToReturn = 1
      } else{
        ck = 0;
      }
    }
  }
  return(ToReturn)
}
#############################################################################################
#############################################################################################
#############################################################################################
GiveModularStructure <- function(size = 9, number_of_modules = 3){
  number_of_modules = size/3
  X = matrix(0, size,size)
  nblck = i = 1
  nmod = number_of_modules
  for(j in 1:ncol(X)){
    if(j < (nmod*nblck + 1)){
      for(k in 0:(nmod - 1)){
        X[i + k,j] = 1
      }
    }else{     
      i = i + nmod;
      for(k in 0:(nmod - 1)){
        X[i + k,j] = 1
      }
      nblck = nblck + 1; }
  }
  for(nb in 2:nblck){
    X[nb*nmod - 1, (nb-1)*nmod - 1] = 1
    X[(nb-1)*nmod - 1, nb*nmod - 1] = 1
  }
  return(X)
}
GiveNestedStructure <- function(size = 9){
  M = matrix(0, size, size)
  for(j in 1:ncol(M)){
    for(i in 1:(nrow(M) - j + 1)){
      M[i,j] = 1
    }
  }
  return(M)
}
############## Create the modular and nested networks
ModularNetwork <- function(mu, size){
  X = 0.1*GiveModularStructure(size = size)
  DStability = 0
  while(DStability == 0){
    ## more randomly interaction strengths
    Aux = ifelse(X > 0, 1, 0)
    X = Aux * replicate(nrow(X), rnorm(nrow(X), 0.4, 0.3))
    X = mu*X
    X = X - diag(diag(X)) + diag(nrow(X))
    D <- diag(1/sqrt(diag(t(X) %*% X)))
    X <- X %*% D
    if(check_if_stable(0.5*(X + t(X))) == 1){ DStability = 1}
  }
  return(X)
}
NestedNetwork <- function(mu, size){
  X = 0.1*GiveNestedStructure(size = size)
  DStability = 0
  while(DStability == 0){
    ### more randomly interaction strengths
    Aux = ifelse(X > 0, 1, 0)
    X = Aux * replicate(nrow(X), rnorm(nrow(X), 0.5, 0.3))
    X = mu*X
    X = X - diag(diag(X)) + diag(nrow(X))
    D <- diag(1/sqrt(diag(t(X) %*% X)))
    X <- X %*% D
    if(check_if_stable(0.5*(X + t(X))) == 1){ DStability = 1}
  }
  return(X)
}
#############################################################################################
#############################################################################################
#############################################################################################
##### Now print, plot and do all the other stuff that is not part of the simulation
PrintHistOnFile <- function(X, nome){
  sink(nome)
  x = c(-1, 1)
  y = X[[2]]; y = y[y!=0]
  for(i in 1:2){
    cat(x[i], y[i], '\n')
  }
  sink()
}

#### To plot the Results ona a pdf file
PltCombined <- function(X, ...){
  if(save_figure == TRUE){
    nome_plot = paste(network_size, '_species_.pdf', sep = '')
    pdf(nome_plot, width = 10, height = 12, family = "Helvetica")
  }
  par(mfrow=c(1,2))
  for(i in 1:length(X)){
    if(i == 1){
      hist(X[[i]], col = 'red',  xaxt='n', main = NULL, 
           ylab = 'number of realizations without extinctions', xlab = '')
    }else{
      hist(X[[i]], col = 'red', xaxt='n', main = NULL, xlab = '', ylab = '')
    }
    title(title_list[i])
    ticks = c(0.999,1.0010)
    axis(side = 1, at = ticks, labels = c('Modular Network', 'Nested Network'))
  }
  if(save_figure == TRUE){
    dev.off()
    comando = paste('open', nome_plot, sep = ' ')
    system(comando)
  }
}
