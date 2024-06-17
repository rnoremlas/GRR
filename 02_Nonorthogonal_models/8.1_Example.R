rm(list=ls())
library(latex2exp)
library(multiColl)

source("functions.txt")
source("GRR1.txt")
source("GRR2.txt")

#################################################################
## Example about 15 Spanish companies: number of employees, E (dependent variable), fixed assets, FA, operating income, OI, and sales, S.
##      see Salmerón, R., A. Rodríguez Sánchez, C. G. García, and J. García Pérez (2020). The vif and mse in raise regression. Mathematics 8(4), 605.

  E <- c(2637,15954,162503,162450,28389,132120,63387,26422,19661,8448,9819,7964,169766,15122,13881)
  FA <- c(44153,9389509,17374000,9723088,95980120,103667000,17588000,48848000,38901000,29109985,25529000,14502621,12639483,26787667,6681800)
  OI <- c(38903,4293386,23703000,23310532,29827663,52574000,22567000,35679000,23449000,14475554,18979000,3717342,32436917,4916125,4472900)
  S <- c(38867,4231043,23649000,23310532,29215382,52036000,22567000,34689000,23184000,14207980,18979000,3709581,31975212,4758244,4472900)
  y = E
  x = cbind(FA, OI, S)
  
  n = length(y)
  p = ncol(x)
  
  ## data standardisation
  
  y_s = (y-mean(y))/sqrt((n-1)*var(y))
  x_s = estandarizar(x) # function in 'functions.txt'
  xx_s = crossprod(x_s)
  
  descomposicion = eigen(xx_s)
  Landa = diag(descomposicion[[1]])
  Gamma = descomposicion[[2]]
  
#################################################################

## OLS estimation
  
  reg = lm(y_s~x_s+0)
  summary(reg) # results of the second column of Table 4
  
  beta = as.double(reg$coef)
  sigma2 = as.double(summary(reg)[[6]]^2)
  GoF_OLS = summary(reg)[[8]]
  
  alfa = t(x_s)%*%y_s
  delta = t(Gamma)%*%alfa
  xi = t(Gamma)%*%beta
  
  K = diag(0,p)
  OLS = GRR1(Gamma,Landa,K,delta,p,sigma2,xi,y_s) 
  OLS # results of the second column of Table 4 using function GRR1 of 'GRR1.txt' 
  MSE_OLS = OLS[[3]]
  
## multicollinearity detection
  
  R = cor(x_s) # the matrix of the correlation 
  R
  detR_OLS = det(R) # and its determinant
  detR_OLS
      cor12_OLS = R[1,2] # for Figure 6
      cor13_OLS = R[1,3]
      cor23_OLS = R[2,3]
  
  CVs(cbind(rep(1,n),x)) # firts row of Table 1 using function of multiColl package
  CVs(cbind(rep(1,n),x_s)) # second row of Table 1 using function of multiColl package 

  CNs(cbind(rep(1,n),x)) # firts row of Table 2 using function of multiColl package
  CNs(cbind(rep(1,n),x_s)) # second row of Table 2 using function of multiColl package

  VIF(cbind(rep(1,n),x)) # third row of Table 1 using function of multiColl package 
  VIF(cbind(rep(1,n),x_s)) # match the above values
  
#################################################################  
  
## k that minimise MSE
  
  k_mins = sigma2/xi^2 # firts column of Table 3
  
  MSEs_min = array(,p)
  var = 0
  for (k_l in k_mins){
    var = var + 1
    ks = array(0, p)
    ks[var] = k_l
    K_mins = diag(ks)
    GRR = GRR1(Gamma,Landa,K_mins,delta,p,sigma2,xi,y_s)
    MSEs_min[var] = GRR[[3]]
  }
  
  cond = xi^2 - sigma2/diag(Landa) < 0 # MSE always lower than OLS?
  
  data.frame(k_mins, MSEs_min, cond) # Table 3
  
#################################################################
  
## election of k
  
  landas = diag(Landa)
  landas
  NC = sqrt(max(landas)/min(landas)) 
  NC
  
  landas_ord = sort(landas)
  NC_min = sqrt(max(landas_ord)/landas_ord[2]) 
  NC_min # minimum CN
  
  knc_min_a = landas_ord[2]-min(landas)
  knc_min_b = max(landas)-min(landas)
  c(knc_min_a, knc_min_b) # interval for k_3 in which CN is equal to its minimum value (2.708444)
  
  #
  
  var = which.min(landas) # variable to be modified
  #kmse_mins = sigma2/xi^2
  kmse_min = k_mins[var] # kmse_mins[var]
  kmse_min
  
  # condition number for regular and generalized ridge regression (function in 'GRR2.txt') for differents value of k
  
    CN(landas, var, kmse_min) # for k = kmse_min
    CN(landas, var, 0.35) # it is observed that in the generalised case it is greater than the following
    CN(landas, var,  knc_min_a) # for k = knc_min_a (the generalised case must coincide with NC_min)
    CN(landas, var, 1) 
    CN(landas, var,  knc_min_b) # for k = knc_min_b (the generalised case must coincide with NC_min)
    CN(landas, var, 2.64) # it can be seen that in the generalised case it is higher than the previous ones
  
#################################################################
## some graphical representations: Figure 4 and 5
    
    #  k which leaves CN below 20 or 10 in the generalised (2) and regular (1) case.
    
    CNmin(20, 1) # function in 'GRR2.txt'
    kcn_menor_20 = CNmin(20, 2) # function in 'GRR2.txt'
    kcn_menor_20 
    
    CNmin(10, 1) # function in 'GRR2.txt'
    kcn_menor_10 = CNmin(10, 2) # function in 'GRR2.txt'
    kcn_menor_10
  
  ## local behaviour of the NC: top of Figure 4 
  
  inicio = 0
  final = 0.03
  salto = 0.0001
  discr = seq(inicio, final, salto)
  
  CN.r = array(NA, length(discr))
  CN.g = array(NA, length(discr))
  j = 1
  for (ki in discr){
    landas.r = landas + ki
    CN.r[j] = sqrt(max(landas.r)/min(landas.r))
    #
    landas.g = landas
    landas.g[var] = landas[var] + ki
    CN.g[j] = sqrt(max(landas.g)/min(landas.g))
    #
    j = j + 1
  }
  
  #CN.r - CN.g
  
  win.graph()
    plot(ts(cbind(CN.r,CN.g), start=inicio, deltat=salto), plot.type="single", col=c("blue","red"), xlab=TeX('$k, k_l$'), ylab=TeX('CN'), lwd = 2, cex.axis=1.5, cex.lab=1.5)
    abline(h=20, lwd=2, col="black", lty=2)
    abline(h=10, lwd=2, col="black", lty=2)
    points(0, NC, type="p", pch = 22, col="black", lwd=4)
    abline(v=kcn_menor_20, lwd=2, col="black", lty=3)
    abline(v=kcn_menor_10, lwd=2, col="black", lty=3)
    savePlot("Figure4top", type="eps")
    #savePlot("Figure4top", type="jpg")
    savePlot("Figure4top", type="pdf")
  dev.off()
  
  ## local behaviour of the MSE: bottom of Figure 4
  
  inicio = 0
  final = 0.0005
  salto = 0.00001
  discr = seq(inicio, final, salto)
  
  MSE_RR = array(,length(discr))
  MSE_GRR = array(,length(discr))
  
  ks = array(0, p)
  
  j = 0
  for (k in discr){
    j = j + 1
    ###
    K_RR = diag(k,p)
    RR = GRR1(Gamma,Landa,K_RR,delta,p,sigma2,xi,y_s)
    #
    MSE_RR[j] = RR[[3]] 
    ###
    ks[var] = k
    K_GRR = diag(ks)
    GRR = GRR1(Gamma,Landa,K_GRR,delta,p,sigma2,xi,y_s) 
    #
    MSE_GRR[j] = GRR[[3]] 
  } 
  
  #MSE_RR - MSE_GRR
  
  win.graph()
    plot(ts(cbind(MSE_RR,MSE_GRR), start=inicio, deltat=salto), plot.type="single", col=c("blue","red"), xlab=TeX('$k, k_l$'), ylab=TeX('MSE'), lwd = 2, cex.axis=1.5, cex.lab=1.5)
    abline(v=kmse_min, lwd=2, col="black", lty=3)
    abline(h=MSE_OLS, lwd=2, col="black", lty=2)
    points(0, MSE_OLS, type="p", pch = 22, col="black", lwd=4)
    abline(v=0.00011, lwd=2, col="black", lty=3)
    savePlot("Figure4bottom", type="eps")
    #savePlot("Figure4bottom", type="jpg")
    savePlot("Figure4bottom", type="pdf")
  dev.off()
  
  CN(landas, var, 0.00011)
  
  ## asymptotic behaviour of the NC 
  
  inicio = 0
  final = 100
  salto = 1
  discr = seq(inicio, final, salto)
  
  CN.r = array(NA, length(discr))
  CN.g = array(NA, length(discr))
  j = 1
  for (ki in discr){
    landas.r = landas + ki
    CN.r[j] = sqrt(max(landas.r)/min(landas.r))
    #
    landas.g = landas
    landas.g[var] = landas[var] + ki
    CN.g[j] = sqrt(max(landas.g)/min(landas.g))
    #
    j = j + 1
  }
  
  win.graph()
    plot(ts(cbind(CN.r,CN.g), start=inicio, deltat=salto), plot.type="single", col=c("blue","red"), xlab=TeX('$k, k_l$'), ylab=TeX('CN'), lwd = 2, cex.axis=1.5, cex.lab=1.5)
    points(0, NC, type="p", pch = 22, col="black", lwd=4)
    abline(v=kcn_menor_20, lwd=2, col="black", lty=3)
    abline(v=kcn_menor_10, lwd=2, col="black", lty=3)
    savePlot("Figure5", type="eps")
    #savePlot("Figure5", type="jpg")  
    savePlot("Figure5", type="pdf")
  dev.off()
  
  
###############################################################
  
## other k elections: results of Table 4
  
  # RR for the value of k of Hoerl, Kennard and Baldwing (results of the third column of Table 4)
  
    k_HKB = as.double((p*sigma2)/crossprod(beta))
    K_HKB = diag(k_HKB,p)
    RR_HKB = GRR1(Gamma,Landa,K_HKB,delta,p,sigma2,xi,y_s) 
    RR_HKB # results of the third column of Table 4
    CN(landas, var, k_HKB)[1]
  
  # RR for the value of k which minimises the MSE (result of the fourth column of Table 4)
    
    j = 1 
    discr = seq(0,1,0.00001)
    mses = array(, length(discr))
    for (k in discr){
      K = diag(k,p)
      Omega = solve(Landa+K)
      Psi = Omega%*%Landa%*%Omega
      I = diag(1,p)
      Teta = (Omega%*%Landa-I)%*%(Omega%*%Landa-I)
      mses[j] = MSEbis(Psi,sigma2,xi,Teta)
      if (j > 1){
        if(mses[j]>mses[j-1]){
          k_minimo = j-1
          break
        }
      }
      j = j + 1
    }
    k_minimo = discr[k_minimo]
    k_minimo < sigma2/max(xi^2) # deber?a salir TRUE (si no sale se debe a una mala estimaci?n de sigma2?)
    K_minimo = diag(k_minimo,p)
    RR_minimo = GRR1(Gamma,Landa,K_minimo,delta,p,sigma2,xi,y_s)
    RR_minimo # result of the fourth column of Table 4
    CN(landas, var, k_minimo)[1]
  
  # RR for the value of k of Hoerl and Kennard (result of the fifth column of Table 4)
  
    k_HK = sigma2/max(xi^2) 
    K_HK = diag(k_HK,p)
    RR_HK = GRR1(Gamma,Landa,K_HK,delta,p,sigma2,xi,y_s) 
    RR_HK # result of the fifth column of Table 4
    CN(landas, var, k_HK)[1]
  
  # GRR for the k_l giving the minimum MSE (coincides with the previous one)
  
    ks = array(0, p)
    ks[var] = kmse_min
    K_min = diag(ks)
    GRR_kl = GRR1(Gamma,Landa,K_min,delta,p,sigma2,xi,y_s)  
    GRR_kl # result of the fifth column of Table 4
    CN(landas, var, kmse_min)[2]
  
  # GRR for the k_i of Hoerl and Kennard
  
    #ks_GHK = sigma2/xi^2
    #K_GHK = diag(ks_GHK)
    #GRR_GHK = GRR1(Gamma,Landa,K_GHK,delta,p,sigma2,xi,y_s) 
    #GRR_GHK
    
  # RR para el k that mitigates multicollinearity (result of the sixth column of Table 4)
    
    K_minimo = diag(0.00652,p)
    RR_multicol = GRR1(Gamma,Landa,K_minimo,delta,p,sigma2,xi,y_s)
    RR_multicol # result of the sixth column of Table 4
    CN(landas, var, 0.00652)[1]
    
  # GRR para el k_l that mitigates multicollinearity (result of the seventh column of Table 4)
    
    ks = array(0, p)
    ks[var] = 0.00651
    K_min = diag(ks)
    GRR_multicol = GRR1(Gamma,Landa,K_min,delta,p,sigma2,xi,y_s) 
    GRR_multicol # result of the seventh column of Table 4
    CN(landas, var, 0.00651)[2]

#################################################################

## some other graphical representations
  
  # Figures 6 and 7
    
  inicio = 0
  final = 100
  salto = 0.1
  discr = seq(inicio, final, salto)
  
  detR.r = array(NA, length(discr))
  cor12.r = array(NA, length(discr))
  cor13.r = array(NA, length(discr))
  cor23.r = array(NA, length(discr))
  
  detR.g = array(NA, length(discr))
  cor12.g = array(NA, length(discr))
  cor13.g = array(NA, length(discr))
  cor23.g = array(NA, length(discr))
  
  j = 1
  for (ki in discr){
    K.r = diag(ki,p)
    Xa.r = as.matrix(rbind(x_s, sqrt(K.r)%*%t(Gamma)))
    R.r = cor(Xa.r)
    cor12.r[j] = R.r[1,2]
    cor13.r[j] = R.r[1,3]
    cor23.r[j] = R.r[2,3]
    detR.r[j] = det(R.r)
    #
    K.g = diag(0, p)
    K.g[var,var] = ki
    Xa.g = as.matrix(rbind(x_s, sqrt(K.g)%*%t(Gamma)))
    R.g = cor(Xa.g)
    cor12.g[j] = R.g[1,2]
    cor13.g[j] = R.g[1,3]
    cor23.g[j] = R.g[2,3]
    detR.g[j] = det(R.g)
    #
    j = j + 1
  }
  
  #cor12.r
  #cor13.r
  #cor23.r
  cor.r = cbind(cor12.r, cor13.r, cor23.r) # top of Figure 6
  win.graph()
    plot(ts(cor.r, start=inicio, deltat=salto), plot.type="single", col=2:(p+1), xlab=TeX('$k$'), ylab="R", lwd=2)
    points(0, cor12_OLS, type="p", pch = 22, col="black", lwd=4)
    points(0, cor13_OLS, type="p", pch = 22, col="black", lwd=4)
    points(0, cor23_OLS, type="p", pch = 22, col="black", lwd=4)
    savePlot("Figure6top", type="eps")
    #savePlot("Figure6top", type="jpg")
    savePlot("Figure6top", type="pdf")
  dev.off()
  
  #cor12.g
  #cor13.g
  #cor23.g
  cor.g = cbind(cor12.g, cor13.g, cor23.g) # bottom of Figure 6
  win.graph()
    plot(ts(cor.g, start=inicio, deltat=salto), plot.type="single", col=2:(p+1), xlab=TeX('$k_l$'), ylab="R", lwd=2)
    points(0, cor12_OLS, type="p", pch = 22, col="black", lwd=4)
    points(0, cor13_OLS, type="p", pch = 22, col="black", lwd=4)
    points(0, cor23_OLS, type="p", pch = 22, col="black", lwd=4)
    savePlot("Figure6bottom", type="eps")
    #savePlot("Figure6bottom", type="jpg")
    savePlot("Figure6bottom", type="pdf")
  dev.off()
  
  #detR.r
  win.graph() # top of Figure 7
    plot(discr, detR.r, type="l", col="blue", xlab=TeX('$k$'), ylab="det(R)", lwd=2)
    points(0, detR_OLS, type="p", pch = 22, col="black", lwd=4)
    savePlot("Figure7top", type="eps")
    #savePlot("Figure7top", type="jpg")
    savePlot("Figure7top", type="pdf")
  dev.off()
  
  #detR.g
  win.graph() # bottom of Figure 7
    plot(discr, detR.g, type="l", col="blue", xlab=TeX('$k_l$'), ylab="det(R)", lwd=2)
    points(0, detR_OLS, type="p", pch = 22, col="black", lwd=4)
    savePlot("Figure7bottom", type="eps")
    #savePlot("Figure7bottom", type="jpg")
    savePlot("Figure7bottom", type="pdf")
  dev.off()
  
#################################################################  

# asymptotic behaviour of the VIF: Tables 5 and 6
  
  inicio = 0
  final = 0.1 
  salto = 0.01 
  discr1 = seq(inicio, final, salto)
  inicio = 10
  final = 100
  salto = 10
  discr2 = seq(inicio, final, salto)
  discr = c(discr1, discr2)
  
  VIF.r = matrix(NA, nrow=length(discr), ncol=p)
  VIF.g = matrix(NA, nrow=length(discr), ncol=p)
  for (h in 1:p){
    j = 1
    for (ki in discr){
      K.r = diag(ki,p)
      Xa.r = as.matrix(rbind(x_s, sqrt(K.r)%*%t(Gamma)))
      R2.r = R2nc(Xa.r[,h],Xa.r[,-h])
      VIF.r[j,h] = 1/(1-R2.r)
      #
      K.g = diag(0, p)
      K.g[var,var] = ki
      Xa.g = as.matrix(rbind(x_s, sqrt(K.g)%*%t(Gamma)))
      R2.g = R2nc(Xa.g[,h],Xa.g[,-h])
      VIF.g[j,h] = 1/(1-R2.g)
      #
      j = j + 1
    }
  }
  
  data.frame(discr, VIF.r) # Table 5
  data.frame(discr, VIF.g) # Table 6
  