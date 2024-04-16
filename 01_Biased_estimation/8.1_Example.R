rm(list=ls())
library(latex2exp) # for special characters in the graphs

#################################################################
## Gorman and Toman data

  library(ridge)
  data(Gorman)
  head(Gorman)

  y = Gorman[,1]
  Y = y - mean(y)
  n = length(y)
  x = Gorman[,-1]
  X = cbind(rep(1, n), x)
  p = dim(X)[2]
  
#################################################################
## OLS
  
  reg = lm(Y~x)
  summary(reg) # results of the second column of Table 2
  
  beta = as.double(reg$coef)
  sigma2 = as.double(summary(reg)[[6]]^2)
  GoF_OLS = as.double(summary(reg)[[8]])
  
#################################################################
## functions for the calculation of 
##      beta(K) from expressions (3) and (5), 
##      MSE from expressions (10)-(11) and (14)
##      GoF from expressions (16) and (21)
  
  source("8.1_functions.txt")

##################################################################################################################################
##################################################################################################################################
## decomposicion: information needed to use above the expressions 
  
  XX = crossprod(X)
  descomposicion = eigen(XX)
  Landa = diag(descomposicion[[1]])
  Gamma = descomposicion[[2]]
  alfa = t(X)%*%Y
  delta = t(Gamma)%*%alfa
  xi = t(Gamma)%*%beta
  
#################################################################
## results of the Table 2

  # a) OLS: obtaining the results for OLS from the GRR function  
  
    K = diag(0,p)
    OLS = GRR(Gamma,Landa,K,delta,p,sigma2,xi,Y)
    OLS # results of the second column of Table 2
  
  # b) RR for the Hoerl, Kennard and Baldwing's k value
  
    k_HKB = as.double((p*sigma2)/crossprod(beta)) 
    k_HKB
    K_HKB = diag(k_HKB,p)
    RR_HKB = GRR(Gamma,Landa,K_HKB,delta,p,sigma2,xi,Y) 
    RR_HKB # results of the third column of Table 2
  
  # b) RR threshold given in Hoerl and Kennard paper
    
    k_HK = sigma2/max(xi^2) 
    k_HK
    K_HK = diag(k_HK,p)
    RR_HK = GRR(Gamma,Landa,K_HK,delta,p,sigma2,xi,Y) 
    RR_HK # results of the fourth column of Table 2
  
  # b) RR for the k that minimizes the MSE
  
    j = 1 
    discr = seq(0,1,0.00001)
    mses = array(, length(discr))
    for (k in discr){ # Algorithm 1: Obtention of the k, k_min, that minimizes MSE
      K = diag(k,p)
      Omega = solve(Landa+K)
      Psi = Omega%*%Landa%*%Omega
      I = diag(1,p)
      Teta = (Omega%*%Landa-I)%*%(Omega%*%Landa-I)
      mses[j] = MSEbis(Psi,sigma2,xi,Teta)
      if (j > 1){
        if(mses[j]>mses[j-1]){
          k_minimum = j-1
          break
        }
      }
      j = j + 1
    }
    k_minimum = discr[k_minimum]
    k_minimum < sigma2/max(xi^2) # It should come out TRUE (if it doesn't come out, is it due to a bad estimate of sigma2?)
    K_minimum = diag(k_minimum,p)
    RR_minimum = GRR(Gamma,Landa,K_minimum,delta,p,sigma2,xi,Y)
    RR_minimum # results of the fifth column of Table 2
  
  # c) GRR for the k_i of Hoerl and Kennard
    
    ks_GHK = sigma2/xi^2
    K_GHK = diag(ks_GHK)
    GRR_GHK = GRR(Gamma,Landa,K_GHK,delta,p,sigma2,xi,Y)  
    GRR_GHK # results of the sixth column of Table 2
    
  # d) GRR for the k values that minimize MSE
  
    k_mins = sigma2/xi^2
    k_mins # results of the first column of Table 1

    MSEs_min = array(,p)
    var = 0
    for (k_l in k_mins){
      var = var + 1
      ks = array(0, p)
      ks[var] = k_l
      K_mins = diag(ks)
      GRRs = GRR(Gamma,Landa,K_mins,delta,p,sigma2,xi,Y)
      MSEs_min[var] = GRRs[[2]]
    }
    MSEs_min # results of the second column of Table 1

    cond = xi^2 - sigma2/diag(Landa) < 0 
    cond # results of the third column of Table 1
    
    data.frame(k_mins, MSEs_min, cond) # table 1
  
    # what is the k_l that gives the minimum MSE
  
      var = which.min(MSEs_min)
      k_min = k_mins[var]
      k_min
      ks = array(0, p)
      ks[var] = k_min
      K_min = diag(ks)
      GRR_kl = GRR(Gamma,Landa,K_min,delta,p,sigma2,xi,Y)  
      GRR_kl # results of the seventh column of Table 2
  
#################################################################
## results of the Table 2
    
  data.frame(0, k_HKB, k_HK, k_minimum, k_min) # first row of Table 2
  data.frame(OLS[[1]], RR_HKB[[1]], RR_HK[[1]], RR_minimum[[1]], GRR_GHK[[1]], GRR_kl[[1]]) # estimation: rows 2 to 12 of Table 2
  data.frame(OLS[[2]], RR_HKB[[2]], RR_HK[[2]], RR_minimum[[2]], GRR_GHK[[2]], GRR_kl[[2]]) # MSE: row 14 of Table 2
  data.frame(OLS[[3]], RR_HKB[[3]], RR_HK[[3]], RR_minimum[[3]], GRR_GHK[[3]], GRR_kl[[3]]) # GoF: row 15 of Table 2

  # row 13 of Table 2
  
    est = data.frame(OLS[[1]], RR_HKB[[1]], RR_HK[[1]], RR_minimum[[1]], GRR_GHK[[1]], GRR_kl[[1]]) 
    data.frame(crossprod(est[[1]]-est[[2]]), crossprod(est[[1]]-est[[3]]), crossprod(est[[1]]-est[[4]]), crossprod(est[[1]]-est[[5]]),   crossprod(est[[1]]-est[[6]]))
    
##################################################################################################################################
##################################################################################################################################
## Figures
  
  inicio = 0
  salto = 0.001
  tope = 1
  discretizacion = seq(inicio, tope, salto) 
  
  BetasK_RR = matrix(,p,length(discretizacion))
  BetasK_GRR = matrix(,p,length(discretizacion))
  norma_RR = array(,length(discretizacion))
  norma_GRR = array(,length(discretizacion))
  MSE_RR = array(,length(discretizacion))
  MSE_GRR = array(,length(discretizacion))
  GoF_RR = array(,length(discretizacion))
  GoF_GRR = array(,length(discretizacion))
  
  ks = array(0, p)
  
  j = 0
  for (k in discretizacion){
    j = j + 1
    ###
    K_RR = diag(k,p)
    RR = GRR(Gamma,Landa,K_RR,delta,p,sigma2,xi,Y)
    #
    BetasK_RR[,j] = RR[[1]]
    norma_RR[j] = crossprod(BetasK_RR[,j])
    MSE_RR[j] = RR[[2]] 
    GoF_RR[j] = RR[[3]] 
    ###
    ks[var] = k
    K_GRR = diag(ks)
    grr = GRR(Gamma,Landa,K_GRR,delta,p,sigma2,xi,Y) 
    #
    BetasK_GRR[,j] = grr[[1]]
    norma_GRR[j] = crossprod(BetasK_GRR[,j])
    MSE_GRR[j] = grr[[2]] 
    GoF_GRR[j] = grr[[3]]
  }
  
  # Figure 3 
  
    plot(discretizacion, MSE_GRR, type="l", col="blue", xlab=TeX('$k_l$'), ylab=TeX('$MSE(\\hat{\\beta(K)})$'), lwd = 2, cex.axis=1, cex.lab=1) 
    MSE_OLS = sigma2*sum(1/diag(Landa))
    asintota = MSE_OLS + xi[var]^2 - (sigma2)/(Landa[var,var])
    abline(h=asintota, lwd=2, col="black", lty=2)
    points(0, MSE_OLS, type="p", pch = 22, col="black", lwd=4)
    abline(v=k_min, lwd=2, col="black", lty=3)
    #savePlot("MSE_GRR", type="eps")
    #savePlot("MSE_GRR", type="jpg")  
  
  # Figure 4 top
  
    plot(ts(t(BetasK_RR), start=inicio, deltat=salto), plot.type="single", col=2:(p+1), xlab=TeX('$k$'), ylab=TeX('$\\hat{\\beta(K)}$'), lwd = 2, cex.axis=1.5, cex.lab=1.5)
    for (i in 1:p) points(0, BetasK_RR[i,1], type="p", pch = 22, col="black", lwd=2)
    #savePlot("betas_RR", type="eps")
    #savePlot("betas_RR", type="jpg")
  
  # Figure 4 bottom
  
    plot(ts(t(BetasK_GRR), start=inicio, deltat=salto), plot.type="single", col=2:(p+1), xlab=TeX('$k_l$'), ylab=TeX('$\\hat{\\beta(K)}$'), lwd = 2, cex.axis=1.5, cex.lab=1.5)
    for (i in 1:p) points(0, BetasK_RR[i,1], type="p", pch = 22, col="black", lwd=2)
    #savePlot("betas_GRR", type="eps")
    #savePlot("betas_GRR", type="jpg")
  
  # Figure 5
  
    plot(ts(cbind(norma_RR,norma_GRR), start=inicio, deltat=salto), plot.type="single", col=c("blue","red"), xlab=TeX('$k, k_l$'), ylab=TeX('$||\\hat{\\beta(K)}||$'), lwd = 2, cex.axis=1.5, cex.lab=1.5)
    points(0, crossprod(beta), type="p", pch = 22, col="black", lwd=4)
    #savePlot("normas", type="eps")
    #savePlot("normas", type="jpg")
  
  # Figure 6
    
    plot(ts(cbind(MSE_RR,MSE_GRR), start=inicio, deltat=salto), plot.type="single", col=c("blue","red"), xlab=TeX('$k, k_l$'), ylab=TeX('$MSE(\\hat{\\beta(K)})$'), lwd = 2, cex.axis=1.5, cex.lab=1.5)
    abline(h=MSE_OLS, lwd=2, col="black", lty=2)
    points(0, MSE_OLS, type="p", pch = 22, col="black", lwd=4)
    #savePlot("MSEs", type="eps")
    #savePlot("MSEs", type="jpg")
  
  # Figure 7
    
    plot(ts(cbind(GoF_RR,GoF_GRR), start=inicio, deltat=salto), plot.type="single", col=c("blue","red"), xlab=TeX('$k, k_l$'), ylab=TeX('$GoF(K)$'), lwd = 2, cex.axis=1.5, cex.lab=1.5)
    points(0, GoF_OLS, type="p", pch = 22, col="black", lwd=4)
    #savePlot("GoFs", type="eps")
    #savePlot("GoFs", type="jpg")
    
    