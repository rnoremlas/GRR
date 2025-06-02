# PRELIMINARIES
rm(list=ls())
library(multiColl)
library(glmnet)
library(ElemStatLearn)
source("8.1_functions.txt")

#################################################################
# DATA SET 1 (DIABETES)
data = read.table("diabetes.csv", header=T, sep=";")
data$SEX <- ifelse(data$SEX == 2, 1, 0)

y = data[,11]
n = length(y)
x = as.matrix(data[,-11])
cte=rep(1, n)
X = cbind(rep(1, n), x)
p = dim(X)[2]

# DIVIDED TRAIN AND TEST SAMPLE
num_iter <- 100
rss_results <- matrix(NA, nrow = num_iter, ncol = 4 * 2)  # 4 métodos x train y test

set.seed(123) 
for (i in 1:num_iter) {
  train_index <- sample(1:nrow(x), size = 0.7 * nrow(x))
  X_train <- x[train_index, ]
  y_train <- y[train_index]
  X_test <- x[-train_index, ]
  y_test <- y[-train_index]
  X_train_cte <- X[train_index, ]
  X_test_cte <- X[-train_index, ]
  
  # OLS
  reg <- lm(y_train ~ X_train)
  beta <- as.double(reg$coef)
  sigma2 <- as.double(summary(reg)[[6]]^2)

  # DESCOMPOSITION
  XX <- crossprod(X_train_cte)
  descomposicion <- eigen(XX)
  Landa <- diag(descomposicion[[1]])
  Gamma <- descomposicion[[2]]
  alfa <- t(X_train_cte) %*% y_train
  delta <- t(Gamma) %*% alfa
  xi <- t(Gamma) %*% beta
  
  k_mins <- sigma2 / xi^2
  var <- which.min(k_mins)
  k_min <- k_mins[var]
  ks <- array(0, p)
  ks[var] <- k_min
  K_min <- diag(ks)
  GRR_kl <- GRR(Gamma, Landa, k_min, delta, p, sigma2, xi, y)
  
  y_pred_GRR_train <- rowSums(sweep(X_train_cte, 2, GRR_kl$Estimation[, 1], "*"))
  rss_results[i, 1] <- sqrt(sum((y_train - y_pred_GRR_train)^2) / n)
  
  y_pred_GRR_test <- rowSums(sweep(X_test_cte, 2, GRR_kl$Estimation[, 1], "*"))
  rss_results[i, 2] <- sqrt(sum((y_test - y_pred_GRR_test)^2) / n)
  
  # Ridge
  ridge_model <- cv.glmnet(X_train, y_train, alpha = 0, grouped = FALSE)
  best_lambda_ridge <- ridge_model$lambda.min
  y_ridge_pred_train <- predict(ridge_model, X_train, s = best_lambda_ridge)
  rss_results[i, 3] <- sqrt(sum((y_train - y_ridge_pred_train)^2) / n)
  y_ridge_pred_test <- predict(ridge_model, X_test, s = best_lambda_ridge)
  rss_results[i, 4] <- sqrt(sum((y_test - y_ridge_pred_test)^2) / n)
  
  # Lasso
  lasso_model <- cv.glmnet(X_train, y_train, alpha = 1, grouped = FALSE)
  best_lambda_lasso <- lasso_model$lambda.min
  y_lasso_pred_train <- predict(lasso_model, X_train, s = best_lambda_lasso)
  rss_results[i, 5] <- sqrt(sum((y_train - y_lasso_pred_train)^2) / n)
  y_lasso_pred_test <- predict(lasso_model, X_test, s = best_lambda_lasso)
  rss_results[i, 6] <- sqrt(sum((y_test - y_lasso_pred_test)^2) / n)
  
  # Elastic Net
  elastic_model <- cv.glmnet(X_train, y_train, alpha = 0.5, grouped = FALSE)
  best_lambda_elastic <- elastic_model$lambda.min
  y_elastic_pred_train <- predict(elastic_model, X_train, s = best_lambda_elastic)
  rss_results[i, 7] <- sqrt(sum((y_train - y_elastic_pred_train)^2) / n)
  y_elastic_pred_test <- predict(elastic_model, X_test, s = best_lambda_elastic)
  rss_results[i, 8] <- sqrt(sum((y_test - y_elastic_pred_test)^2) / n)
}

# MEAN
rss_means <- colMeans(rss_results, na.rm = TRUE)
rss_sds <- apply(rss_results, 2, sd, na.rm = TRUE)
rss_table_diabetes <- data.frame(
  Method = c("GRR", "Ridge", "Lasso", "Elastic Net"),
  RSS_Train_Mean = rss_means[c(1, 3, 5, 7)],
  RSS_Train_SD = rss_sds[c(1, 3, 5, 7)],
  RSS_Test_Mean = rss_means[c(2, 4, 6, 8)],
  RSS_Test_SD = rss_sds[c(2, 4, 6, 8)]
)



#################################################################)
## DATA SET 2 (PROSTATE)
data=data(prostate)#SecarganlosdatosProstatedelalibreriaElemStatLearn
attach(prostate)
head(prostate)
summary(prostate)

y = lpsa
n = length(y)
x = as.matrix(cbind(prostate[,-((ncol(prostate)-1):ncol(prostate))]))
cte=rep(1, n)
X = cbind(rep(1, n), x)
p = dim(X)[2]

## DIVIDED TRAIN AND TEST SAMPLE
num_iter <- 100
rss_results <- matrix(NA, nrow = num_iter, ncol = 4 * 2)  # 4 métodos x train y test

set.seed(123) 
for (i in 1:num_iter) {
  train_index <- sample(1:nrow(x), size = 0.7 * nrow(x))
  X_train <- x[train_index, ]
  y_train <- y[train_index]
  X_test <- x[-train_index, ]
  y_test <- y[-train_index]
  X_train_cte <- X[train_index, ]
  X_test_cte <- X[-train_index, ]
  
  
  # OLS
  reg <- lm(y_train ~ X_train)
  beta <- as.double(reg$coef)
  sigma2 <- as.double(summary(reg)[[6]]^2)
  
  # DESCOMPOSITION
  XX <- crossprod(X_train_cte)
  descomposicion <- eigen(XX)
  Landa <- diag(descomposicion[[1]])
  Gamma <- descomposicion[[2]]
  alfa <- t(X_train_cte) %*% y_train
  delta <- t(Gamma) %*% alfa
  xi <- t(Gamma) %*% beta
  
  k_mins <- sigma2 / xi^2
  var <- which.min(k_mins)
  k_min <- k_mins[var]
  ks <- array(0, p)
  ks[var] <- k_min
  K_min <- diag(ks)
  GRR_kl <- GRR(Gamma, Landa, k_min, delta, p, sigma2, xi, y)
  
  y_pred_GRR_train <- rowSums(sweep(X_train_cte, 2, GRR_kl$Estimation[, 1], "*"))
  rss_results[i, 1] <- sqrt(sum((y_train - y_pred_GRR_train)^2) / n)
  
  y_pred_GRR_test <- rowSums(sweep(X_test_cte, 2, GRR_kl$Estimation[, 1], "*"))
  rss_results[i, 2] <- sqrt(sum((y_test - y_pred_GRR_test)^2) / n)
  
  # Ridge
  ridge_model <- cv.glmnet(X_train, y_train, alpha = 0, grouped = FALSE)
  best_lambda_ridge <- ridge_model$lambda.min
  y_ridge_pred_train <- predict(ridge_model, X_train, s = best_lambda_ridge)
  rss_results[i, 3] <- sqrt(sum((y_train - y_ridge_pred_train)^2) / n)
  y_ridge_pred_test <- predict(ridge_model, X_test, s = best_lambda_ridge)
  rss_results[i, 4] <- sqrt(sum((y_test - y_ridge_pred_test)^2) / n)
  
  # Lasso
  lasso_model <- cv.glmnet(X_train, y_train, alpha = 1, grouped = FALSE)
  best_lambda_lasso <- lasso_model$lambda.min
  y_lasso_pred_train <- predict(lasso_model, X_train, s = best_lambda_lasso)
  rss_results[i, 5] <- sqrt(sum((y_train - y_lasso_pred_train)^2) / n)
  y_lasso_pred_test <- predict(lasso_model, X_test, s = best_lambda_lasso)
  rss_results[i, 6] <- sqrt(sum((y_test - y_lasso_pred_test)^2) / n)
  
  # Elastic Net
  elastic_model <- cv.glmnet(X_train, y_train, alpha = 0.5, grouped = FALSE)
  best_lambda_elastic <- elastic_model$lambda.min
  y_elastic_pred_train <- predict(elastic_model, X_train, s = best_lambda_elastic)
  rss_results[i, 7] <- sqrt(sum((y_train - y_elastic_pred_train)^2) / n)
  y_elastic_pred_test <- predict(elastic_model, X_test, s = best_lambda_elastic)
  rss_results[i, 8] <- sqrt(sum((y_test - y_elastic_pred_test)^2) / n)
}

# MEAN
rss_means <- colMeans(rss_results, na.rm = TRUE)
rss_sds <- apply(rss_results, 2, sd, na.rm = TRUE)
rss_table_prostate <- data.frame(
  Method = c("GRR", "Ridge", "Lasso", "Elastic Net"),
  RSS_Train_Mean = rss_means[c(1, 3, 5, 7)],
  RSS_Train_SD = rss_sds[c(1, 3, 5, 7)],
  RSS_Test_Mean = rss_means[c(2, 4, 6, 8)],
  RSS_Test_SD = rss_sds[c(2, 4, 6, 8)]
)

#######################################################################
## DATA SET 3 (ABALONE)

url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/abalone/abalone.data"
abalone <- read.csv(url, header = FALSE)
colnames(abalone) <- c("Sex", "Length", "Diameter", "Height", 
                       "WholeWeight", "ShuckedWeight", "VisceraWeight", 
                       "ShellWeight", "Rings")
abalone$Sex <- as.factor(abalone$Sex)
sex_dummies <- model.matrix(~ Sex - 1, data = abalone)
abalone <- cbind(abalone, sex_dummies)
attach(abalone)

y = Rings
n = length(y)
x = as.matrix(cbind(Length, Diameter, Height, WholeWeight, ShuckedWeight, VisceraWeight, ShellWeight, SexF, SexI))
cte=rep(1, n)
X = cbind(rep(1, n), x)
p = dim(X)[2]

## DIVIDED TRAIN AND TEST SAMPLE
num_iter <- 100
rss_results <- matrix(NA, nrow = num_iter, ncol = 4 * 2)  # 4 métodos x train y test
set.seed(123)
for (i in 1:num_iter) {
  train_index <- sample(1:nrow(x), size = 0.7 * nrow(x))
  X_train <- x[train_index, ]
  y_train <- y[train_index]
  X_test <- x[-train_index, ]
  y_test <- y[-train_index]
  X_train_cte <- X[train_index, ]
  X_test_cte <- X[-train_index, ]
  
  
  # OLS
  reg <- lm(y_train ~ X_train)
  beta <- as.double(reg$coef)
  sigma2 <- as.double(summary(reg)[[6]]^2)
  cn<-CN(X_train)
  
  # DESCOMPOSITION
  XX <- crossprod(X_train_cte)
  descomposicion <- eigen(XX)
  Landa <- diag(descomposicion[[1]])
  Gamma <- descomposicion[[2]]
  alfa <- t(X_train_cte) %*% y_train
  delta <- t(Gamma) %*% alfa
  xi <- t(Gamma) %*% beta
  
  k_mins <- sigma2 / xi^2
  var <- which.min(k_mins)
  k_min <- k_mins[var]
  ks <- array(0, p)
  ks[var] <- k_min
  K_min <- diag(ks)
  GRR_kl <- GRR(Gamma, Landa, k_min, delta, p, sigma2, xi, y)
  
  y_pred_GRR_train <- rowSums(sweep(X_train_cte, 2, GRR_kl$Estimation[, 1], "*"))
  rss_results[i, 1] <- sqrt(sum((y_train - y_pred_GRR_train)^2) / n)
  
  y_pred_GRR_test <- rowSums(sweep(X_test_cte, 2, GRR_kl$Estimation[, 1], "*"))
  rss_results[i, 2] <- sqrt(sum((y_test - y_pred_GRR_test)^2) / n)
  
  # Ridge
  ridge_model <- cv.glmnet(X_train, y_train, alpha = 0, grouped = FALSE)
  best_lambda_ridge <- ridge_model$lambda.min
  y_ridge_pred_train <- predict(ridge_model, X_train, s = best_lambda_ridge)
  rss_results[i, 3] <- sqrt(sum((y_train - y_ridge_pred_train)^2) / n)
  y_ridge_pred_test <- predict(ridge_model, X_test, s = best_lambda_ridge)
  rss_results[i, 4] <- sqrt(sum((y_test - y_ridge_pred_test)^2) / n)
  
  # Lasso
  lasso_model <- cv.glmnet(X_train, y_train, alpha = 1, grouped = FALSE)
  best_lambda_lasso <- lasso_model$lambda.min
  y_lasso_pred_train <- predict(lasso_model, X_train, s = best_lambda_lasso)
  rss_results[i, 5] <- sqrt(sum((y_train - y_lasso_pred_train)^2) / n)
  y_lasso_pred_test <- predict(lasso_model, X_test, s = best_lambda_lasso)
  rss_results[i, 6] <- sqrt(sum((y_test - y_lasso_pred_test)^2) / n)
  
  # Elastic Net
  elastic_model <- cv.glmnet(X_train, y_train, alpha = 0.5, grouped = FALSE)
  best_lambda_elastic <- elastic_model$lambda.min
  y_elastic_pred_train <- predict(elastic_model, X_train, s = best_lambda_elastic)
  rss_results[i, 7] <- sqrt(sum((y_train - y_elastic_pred_train)^2) / n)
  y_elastic_pred_test <- predict(elastic_model, X_test, s = best_lambda_elastic)
  rss_results[i, 8] <- sqrt(sum((y_test - y_elastic_pred_test)^2) / n)
}

# Mean
rss_means <- colMeans(rss_results, na.rm = TRUE)
rss_sds <- apply(rss_results, 2, sd, na.rm = TRUE)
rss_table_abalone <- data.frame(
  Method = c("GRR", "Ridge", "Lasso", "Elastic Net"),
  RSS_Train_Mean = rss_means[c(1, 3, 5, 7)],
  RSS_Train_SD = rss_sds[c(1, 3, 5, 7)],
  RSS_Test_Mean = rss_means[c(2, 4, 6, 8)],
  RSS_Test_SD = rss_sds[c(2, 4, 6, 8)]
)

########################################################
## DATA SET 4 (BODY FAT)
bodyfat = read.table("bodyfat.csv", header=T, sep=",")
attach(bodyfat)

y = BodyFat
n = length(y)
x = as.matrix(cbind(Density, Age, Weight, Height, Neck, Chest, Abdomen, Hip, Thigh, Knee, Ankle, Biceps, Forearm, Wrist))
cte=rep(1, n)
X = cbind(rep(1, n), x)
p = dim(X)[2]

## DIVIDED TRAIN AND TEST SAMPLE
num_iter <- 100
rss_results <- matrix(NA, nrow = num_iter, ncol = 4 * 2)  # 4 métodos x train y test

set.seed(123)
for (i in 1:num_iter) {
  train_index <- sample(1:nrow(x), size = 0.7 * nrow(x))
  X_train <- x[train_index, ]
  y_train <- y[train_index]
  X_test <- x[-train_index, ]
  y_test <- y[-train_index]
  X_train_cte <- X[train_index, ]
  X_test_cte <- X[-train_index, ]
  
  
  # OLS
  reg <- lm(y_train ~ X_train)
  beta <- as.double(reg$coef)
  sigma2 <- as.double(summary(reg)[[6]]^2)
  
  # DESCOMPOSITION
  XX <- crossprod(X_train_cte)
  descomposicion <- eigen(XX)
  Landa <- diag(descomposicion[[1]])
  Gamma <- descomposicion[[2]]
  alfa <- t(X_train_cte) %*% y_train
  delta <- t(Gamma) %*% alfa
  xi <- t(Gamma) %*% beta
  
  k_mins <- sigma2 / xi^2
  var <- which.min(k_mins)
  k_min <- k_mins[var]
  ks <- array(0, p)
  ks[var] <- k_min
  K_min <- diag(ks)
  GRR_kl <- GRR(Gamma, Landa, k_min, delta, p, sigma2, xi, y)
  
  y_pred_GRR_train <- rowSums(sweep(X_train_cte, 2, GRR_kl$Estimation[, 1], "*"))
  rss_results[i, 1] <- sqrt(sum((y_train - y_pred_GRR_train)^2) / n)
  
  y_pred_GRR_test <- rowSums(sweep(X_test_cte, 2, GRR_kl$Estimation[, 1], "*"))
  rss_results[i, 2] <- sqrt(sum((y_test - y_pred_GRR_test)^2) / n)
  
  # Ridge
  ridge_model <- cv.glmnet(X_train, y_train, alpha = 0, grouped = FALSE)
  best_lambda_ridge <- ridge_model$lambda.min
  y_ridge_pred_train <- predict(ridge_model, X_train, s = best_lambda_ridge)
  rss_results[i, 3] <- sqrt(sum((y_train - y_ridge_pred_train)^2) / n)
  y_ridge_pred_test <- predict(ridge_model, X_test, s = best_lambda_ridge)
  rss_results[i, 4] <- sqrt(sum((y_test - y_ridge_pred_test)^2) / n)
  
  # Lasso
  lasso_model <- cv.glmnet(X_train, y_train, alpha = 1, grouped = FALSE)
  best_lambda_lasso <- lasso_model$lambda.min
  y_lasso_pred_train <- predict(lasso_model, X_train, s = best_lambda_lasso)
  rss_results[i, 5] <- sqrt(sum((y_train - y_lasso_pred_train)^2) / n)
  y_lasso_pred_test <- predict(lasso_model, X_test, s = best_lambda_lasso)
  rss_results[i, 6] <- sqrt(sum((y_test - y_lasso_pred_test)^2) / n)
  
  # Elastic Net
  elastic_model <- cv.glmnet(X_train, y_train, alpha = 0.5, grouped = FALSE)
  best_lambda_elastic <- elastic_model$lambda.min
  y_elastic_pred_train <- predict(elastic_model, X_train, s = best_lambda_elastic)
  rss_results[i, 7] <- sqrt(sum((y_train - y_elastic_pred_train)^2) / n)
  y_elastic_pred_test <- predict(elastic_model, X_test, s = best_lambda_elastic)
  rss_results[i, 8] <- sqrt(sum((y_test - y_elastic_pred_test)^2) / n)
}

# MEAN
rss_means <- colMeans(rss_results, na.rm = TRUE)
rss_sds <- apply(rss_results, 2, sd, na.rm = TRUE)
rss_table_bodyfat <- data.frame(
  Method = c("GRR", "Ridge", "Lasso", "Elastic Net"),
  RSS_Train_Mean = rss_means[c(1, 3, 5, 7)],
  RSS_Train_SD = rss_sds[c(1, 3, 5, 7)],
  RSS_Test_Mean = rss_means[c(2, 4, 6, 8)],
  RSS_Test_SD = rss_sds[c(2, 4, 6, 8)]
)

########################################################
## RESULTS
print(rss_table_diabetes)
print(rss_table_prostate)
print(rss_table_abalone)
print(rss_table_bodyfat)
