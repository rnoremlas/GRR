rm(list=ls())

#####################################################################
## Gorman and Toman data

# data

  library(ridge)
  data(Gorman)
  head(Gorman)

  y = Gorman[,1]
  Y = y - mean(y)
  n = length(y)
  x = Gorman[,-1]
  X = cbind(rep(1, n), x)
  p = ncol(X)
  
# OLS estimation
  
  reg = lm(Y~x)
  summary(reg) 
  
  beta = as.double(reg$coef)
  sigma2 = as.double(summary(reg)[[6]]^2)
  GoF_OLS = as.double(summary(reg)[[8]])
  
#####################################################################
## functions for the calculation of 
##      beta(K) from expressions (3) and (5), 
##      MSE from expressions (10)-(11) and (14)
##      GoF from expressions (16) and (21)
  
  source("8.1_functions.txt")
  
## k_l calculation from k_l = \sigma^2/xi^2 (kl_calculation function)
## bootstrap inference (GRR_bootstrap function for ridge and bootstrap.lasso function for lasso)
  
  source("8.2_functions.txt")
  
#####################################################################  
# estimations
  
  XX = crossprod(X)
  descomposition = eigen(XX)
  Landa = diag(descomposition[[1]])
  Gamma = descomposition[[2]]
  alfa = t(X)%*%Y
  delta = t(Gamma)%*%alfa
  xi = t(Gamma)%*%beta
  
  ite = 5000 # bootstrap
  
  # kas
  
    k_OLS = 0 # OLS
    k_HKB = as.double((p*sigma2)/crossprod(beta))  # Hoerl, Kennard and Baldwing's k value
    k_HK = sigma2/max(xi^2) # Hoerl and Kennard's k value
    kls = kl_calculation(sigma2,xi,p,Gamma,Landa,K_mins,delta,Y) # function in "8.2_functions.txt"
      kls # results of Table 1
      pos_kl = which.min(kls[,2]) # k_l for MSE min
      #pos_kl = which.min(kls[kls[,3]==TRUE,][,2]) # k_l for MSE min and cond==TRUE
      k_l = kls[pos_kl,1]
    k_minimum = kmin_calculation(leap=0.00001,p,Landa,K,sigma2,xi)
      
    kas = c(k_OLS, k_HKB, k_HK, k_minimum, 0, k_l)
    GRR_kas = c(0,0,0,0,1, 2) # value 0 for regular ridge regression
  
    table0 = "k values"
    table1 = "beta 1"
    for (i in 2:p) table1 = c(table1, paste("beta",i))
    table2 = "MSE"
    table3 = "GoF"
    table_boots = c()
    j = 0
    for (k in kas){
      j = j + 1
      if (GRR_kas[j] == 0) K = diag(k,p)
      if (GRR_kas[j] == 1) K = diag(kls[,1])
      if (GRR_kas[j] == 2) {
        ks = array(0, p)
        ks[pos_kl] = k
        K = diag(ks) 
      }
      GRR_est = GRR(Gamma,Landa,K,delta,p,sigma2,xi,Y)
      table0 = cbind(table0, k)
      table1 = cbind(table1, GRR_est[[1]])
      table2 = cbind(table2, GRR_est[[2]])   
      table3 = cbind(table3, GRR_est[[3]])   
      boots = GRR_bootstrap(ite,n,p,beta,K,seed=0,estimation=1,option=2) # function in "8.2_functions.txt"
      table_boots = rbind(table_boots, "", boots)
    }  
    table = rbind(table0, table1, table2, table3)
    colnames(table) = c("", "k_OLS", "k_HKB", "k_HK", "k_min", "k_i", "k_l")
    table # results of Table 2
    table_boots
    
##################################################################### 
## comparison with existing R packages for the ridge estimator    
    
# genridge
    
    library("genridge") # install.packages("genridge")
    reg_ridge = ridge(Y, X[,-1], lambda = kas[GRR_kas==0])
    coef(reg_ridge) # does not coincide with the estimate given in line 95 (Table 2 of the paper)
    
    # GRR function with standardized data
    
      Xest = standardize(X[,-1])
      Yest = (y - mean(y))/sqrt((n-1)*var(y))
      
      reg_est = lm(Yest~Xest+0)
      beta_est = reg_est$coef
      sigma2est = as.double(summary(reg_est)[[6]]^2)
      
      XX = crossprod(Xest)
      descomposition = eigen(XX)
      Landa = diag(descomposition[[1]])
      Gamma = descomposition[[2]]
      alfa = t(Xest)%*%Yest
      delta = t(Gamma)%*%alfa
      xi = t(Gamma)%*%beta_est
      GRR(Gamma,Landa,diag(0,p-1),delta,p-1,sigma2est,xi,Yest) # does not coincide with the estimate given by 'genridge' in line 105 
    
# lrmest
      
    library("lrmest") # install.packages("lrmest")
    rid(Y~X+0, k = k_OLS) # coincides with the second column of Table 2
    rid(Y~X+0, k = k_HKB) # coincides with the third column of Table 2
    rid(Y~X+0, k = k_HK) # coincides with the fourth column of Table 2
    rid(Y~X+0, k = k_minimum) # coincides with the fifth column of Table 2
    
    rid(Yest~Xest+0, k = k_OLS) # coincides with the results of line 123
    
# lmridge
    
    library("lmridge") # install.packages("lmridge")
    data_frame = data.frame(Y, x)
    
    reg_lmridge = lmridge(Y~x, data = data_frame, K = kas[GRR_kas==0], scaling="non")    
    summary(reg_lmridge) # only coincides for k=0, for the rest it is similar
    kest(reg_lmridge) # differents estimations of k: the value given for k_HBK is different from the one calculated in line 57
    rstats1(reg_lmridge) # differents values for MSE
    
    reg_lmridge = lmridge(Y~x+0, data = data_frame, K = 0, scaling="scaled")    
    summary(reg_lmridge) # estimates similar to those given in row 105 by 'genridge' and different from 123/133
    kest(reg_lmridge) # differents estimations of k: the value given for k_HBK is different from the one calculated in line 57 and 142
    rstats1(reg_lmridge) # differents values for MSE from Table 2 and line 143

# ridge
    
    library("ridge") # install.packages("ridge")
    reg_linearRidge = linearRidge(Y~x, data = data_frame, lambda = k_OLS, scaling="none")
    summary(reg_linearRidge) # coincides with the second column of Table 2
    reg_linearRidge = linearRidge(Y~x, data = data_frame, lambda = k_HKB, scaling="none")
    summary(reg_linearRidge) # coincides with the estimation given in line 141 (similar to the third column of Table 2)
    reg_linearRidge = linearRidge(Y~x, data = data_frame, lambda = k_HK, scaling="none")
    summary(reg_linearRidge) # coincides with the estimation given in line 141 (similar to the fourth column of Table 2)
    reg_linearRidge = linearRidge(Y~x, data = data_frame, lambda = k_minimum, scaling="none")
    summary(reg_linearRidge) # coincides with the estimation given in line 141 (similar to the fifth column of Table 2)
    
    reg_linearRidge = linearRidge(Y~x, data = data_frame, lambda = k_OLS, scaling="corrForm") 
    summary(reg_linearRidge) # does not coincide with the estimate given in lines 107, 125/135, 146
    reg_linearRidge = linearRidge(Y~x, data = data_frame, lambda = k_OLS, scaling="scale") 
    summary(reg_linearRidge) # does not coincide with the estimate given in lines 105 and 123/133, coincide with the estimation in line 146
    
# glmnet    
    
    library("glmnet") # install.packages("glmnet")
    ridge_glmnet = glmnet(x, Y, alpha = 0, lambda = kas[GRR_kas==0], standardize = FALSE) # alpha = 0 denotes the ridge estimation
    rbind(ridge_glmnet$lambda, as.matrix(ridge_glmnet$beta))
    ridge_glmnet = glmnet(x, Y, alpha = 0, lambda = kas[GRR_kas==0]) # alpha = 0 denotes the ridge estimation
    rbind(ridge_glmnet$lambda, as.matrix(ridge_glmnet$beta))
    
##################################################################### 
## comparison with existing R packages for the lasso estimator 
    
    set.seed(2024)
    
    cv_output = cv.glmnet(x, Y, alpha=1, lambda = 10^seq(-3,3,0.1))
    best_landa = cv_output$lambda.min
    best_landa
    
    lasso_best = glmnet(x, Y, alpha=1, lambda = best_landa)
    coef(lasso_best)
    
    