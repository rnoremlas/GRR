rm(list=ls())
source("functions.txt")
source("GRR2.txt")
source("GRR3.txt")

#####################################################################
## Data 

names(longley)
n = nrow(longley) 

attach(longley) 
   Y = (Employed - mean(Employed))/sqrt((n-1)*var(Employed))
   X = cbind(GNP.deflator, GNP, Population) # con esta opción gana k_l
   X = cbind(GNP.deflator, GNP, Unemployed, Armed.Forces, Population)
   X = standardize(X)
detach(longley)

p = ncol(X)
   
#####################################################################
# OLS estimation
  
  reg = lm(Y~X+0)
  summary(reg) # results of table 7
  
  beta = as.double(reg$coef)
  sigma2 = as.double(summary(reg)[[6]]^2)
  GoF_OLS = as.double(summary(reg)[[8]])

  # decomposition
  
  XX = crossprod(X)
  descomposition = eigen(XX)
  Landa = diag(descomposition[[1]])
  Gamma = descomposition[[2]]
  alfa = t(X)%*%Y
  delta = t(Gamma)%*%alfa
  xi = t(Gamma)%*%beta

#####################################################################  
    
  ite = 10000 # bootstrap
  
  # estimates of parameter k: regular case
  
  k_OLS = 0 # OLS
  k_HKB = as.double((p*sigma2)/crossprod(beta))  # Hoerl, Kennard and Baldwing's k value
  k_HK = sigma2/max(xi^2) # Hoerl and Kennard's k value
  k_minimum = kmin_calculation(leap=0.00001,p,Landa,K,sigma2,xi) # k that truly minimize MSE
  
  # estimates of parameter k: generalized case
  
  # first part of GRR
  
    kls = kl_calculation(sigma2,xi,p,Gamma,Landa,K_mins,delta,Y) # function in "GRR3.txt"
      kls # results of Table 8
      pos_kl = which.min(kls[,2]) # k_l for MSE min
      k_l = kls[pos_kl,1]
  
  # second part of GRR (actual paper)
    
    landas = diag(Landa) # k that mitigate multicollinearity
    landas
    NC = sqrt(max(landas)/min(landas)) 
    NC
  
    landas_ord = sort(landas)
    NC_min = sqrt(max(landas_ord)/landas_ord[2]) 
    NC_min # minimum CN
    
    knc_min_a = landas_ord[2]-min(landas)
    knc_min_b = max(landas)-min(landas)
    c(knc_min_a, knc_min_b) # interval for k_5 in which CN is equal to its minimum value (3.041291)
  
  #
  
  var = which.min(landas) # variable to be modified (coincide with pos_kl?)
  
    # condition number for regular and generalized ridge regression (function in 'GRR2.txt') for differents value of k

    CN(landas, var, k_OLS) # for k = 0
    CN(landas, var,  knc_min_a) # for k = knc_min_a (the generalised case must coincide with NC_min)
    CN(landas, var,  knc_min_b) # for k = knc_min_b (the generalised case must coincide with NC_min)
  
  # model estimation
  
  kas = c(k_OLS, k_HKB, k_HK, k_minimum, k_l, knc_min_a, knc_min_b) 
  GRR_kas = c(0,0,0,0,1,2,2) # value 0 for regular ridge regression, 1/2 for particular case of generalized ridge regression
  
  table0 = "k values"
  table1 = "beta 1"
  for (i in 2:p) table1 = c(table1, paste("beta",i))
  table2 = "MSE"
  table3 = "GoF"
  table4 = "CN"
  table_boots = c()
  j = 0
  for (k in kas){
    j = j + 1
    if (GRR_kas[j] == 0) K = diag(k,p)
    if (GRR_kas[j] == 1) {
      ks = array(0, p)
      ks[pos_kl] = k
      K = diag(ks) 
    }
    if (GRR_kas[j] == 2) {
      ks = array(0, p)
      ks[var] = k
      K = diag(ks) 
    }
    print(K)
    GRR_est = GRR(Gamma,Landa,K,delta,p,sigma2,xi,Y)
    table0 = cbind(table0, k)
    table1 = cbind(table1, round(GRR_est[[1]], digits=4))
    table2 = cbind(table2, round(GRR_est[[2]], digits=4))   
    table3 = cbind(table3, round(GRR_est[[3]], digits=4))
    cn_index = 2
    if (GRR_kas[j] == 0) cn_index = 1
    table4 = cbind(table4, round(CN(landas, var, k)[[cn_index]], digits=4)) 
    #
    boots = GRR_bootstrap(ite,n,p,beta,K,seed=0,estimation=1,option=1) # function in "GRR3.txt"
    table_boots = rbind(table_boots, "", round(boots, digits=4))
  }  
  table = rbind(table0, table1, table2, table3, table4)
  colnames(table) = c("", "k_OLS", "k_HKB", "k_HK", "k_min", "k_l", "k_NCa", "k_NCb")
  table # results of Table 9
  table_boots # results of Table 9

# results of Table 9 for LaTeX
  
  for (i in 1:nrow(table)){
    tabla = c(table[i,1])
    for (j in 2:ncol(table)) {
      tabla = c(tabla, "&", table[i,j])
    }
    cat(tabla, "\\\\", "\n", sep=" ", file="_resultados_8_2.txt", append=TRUE)  
    if ((i > 1)&(i<=(p+1))) {
      tabla_boots = c()
      for (j in seq(i,nrow(table_boots),p+2)) {
        tabla_boots = c(tabla_boots, "&", table_boots[j,])
      }
      cat(tabla_boots, "\\\\", "\n", sep=" ", file="_resultados_8_2.txt", append=TRUE)  
    }
    if (i == p+3) {
      tabla_boots = c()
      for (j in seq(i-1,nrow(table_boots),p+2)) {
        tabla_boots = c(tabla_boots, "&", table_boots[j,])
      }
      cat(tabla_boots, "\\\\", "\n", sep=" ", file="_resultados_8_2.txt", append=TRUE)  
    }
  }
  
#####################################################################

## comparison with lmridge and ridge packages
  
  kas_RR = c(k_OLS, k_HKB, k_HK)
  kas_RR

  names(longley)
  attach(longley)   
    data_frame = data.frame(Employed, GNP.deflator, GNP, Unemployed, Armed.Forces, Population)
  detach(longley)
  
  # lmridge
  
    library("lmridge") # install.packages("lmridge")
    # The scaling option "sc" scales the predictors to correlation form, such that the correlation matrix has unit diagonal elements

      reg_lmridge = lmridge(Employed~GNP.deflator+GNP+Unemployed+Armed.Forces+Population+0, data = data_frame, K = kas_RR, scaling="sc")   
      summary(reg_lmridge) # results of Table 10 to 12
        kest(reg_lmridge) 
        rstats1(reg_lmridge) 

  
  # ridge
  
    library("ridge") # install.packages("ridge")
    # The scaling option "corrform" (the default) scales the predictors to correlation form, such that the correlation matrix has unit diagonal. 
  
      reg_linearRidge = linearRidge(Employed~GNP.deflator+GNP+Unemployed+Armed.Forces+Population+0, data = data_frame, lambda = k_OLS, scaling="corrForm") 
        summary(reg_linearRidge)# results of Table 10
        
      reg_linearRidge = linearRidge(Employed~GNP.deflator+GNP+Unemployed+Armed.Forces+Population+0, data = data_frame, lambda = k_HKB, scaling="corrForm")
        summary(reg_linearRidge) # results of Table 11
        
      reg_linearRidge = linearRidge(Employed~GNP.deflator+GNP+Unemployed+Armed.Forces+Population+0, data = data_frame, lambda = k_HK, scaling="corrForm")
        summary(reg_linearRidge) # results of Table 12
  