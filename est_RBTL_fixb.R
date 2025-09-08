rm(list=ls())
set.seed(2)
setwd("~/Desktop")
source("~/Desktop/stats_2/est_fun_RBTL_set1.R")
source("~/Desktop/stats_2/est_fun_RBTL_set2.R")
source("~/Desktop/stats_2/est_fun_RBTL_set3.R")
source("~/Desktop/stats_2/est_fun_RBTL_set4.R")
library(maxLik)
library(stats4)
library(numDeriv)
library(miscTools)

MLE_fun_param1B <- function(params, x){
  sigma <- params[1]
  a <- params[2]
  c <- params[3]
  b <- 0.4  # Assuming sigma is fixed
  
  f_x <- dgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  term1 <- 4 * b * a * c *(x^(c-1)) * (x^c)^(a-1) *(1 - (x^c/(2+x^c))^a) /(2+x^c)^(a+1)
  term2 <- 1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2
  
  f_x <- f_x * term1/term2
  
  # # Step 1: Calculate the internal expression inside pgamma
  # term <- (1 - (1 - ((1 - (1 + x^c)^(-1)) / (1 + (1 + x^c)^(-1)))^a)^2)^b
  # z <- -log(term)
  # 
  # # Step 2: Calculate the derivative of z with respect to x
  # dz_dx <- b * 2 * (1 - term)^(b-1) * term * a * (1 - ((1 - (1 + x^c)^(-1)) / (1 + (1 + x^c)^(-1)))^(a-1)) * c * x^(c-1) / ((1 + (1 + x^c)^(-1))^2)
  # 
  # # Step 3: Calculate the PDF
  # f_x <- dgamma(z, shape = sigma) * dz_dx
  
  # term1 <- 4*a*b/gamma(sigma)
  # term2 <- 1 - ((1-(1+x^c)^(-1)) / (1+(1-(1+x^c)^(-1))))^a
  # term3 <- (-log((1 - ( 1 - ((1-(1+x^c)^(-1)) / (1+(1-(1+x^c)^(-1))))^a )^2)^b))^(sigma -1)
  # term4 <- (1 - (1 - ((1-(1+x^c)^(-1)) / (1+(1-(1+x^c)^(-1))))^a)^2)^(b-1)
  # term5 <- c * (x^(c-1)) * ((1+x^c)^(-2)) * ((1- (1+x^c)^(-1))^(a-1))
  # term6 <- ((1- (1+x^c)^(-1))^(a+1))
  # 
  # f_x <- term1 * term2 * term3 * term4 * term5 / term6
  
  # Negative log-likelihood function
  LnL <- - sum(log(f_x))
  return(LnL)
}

# Quantile Function
quantile_function <- function(p, sigma, b, a, c) {
  # gamma_inv <- qgamma(gamma(sigma) * (1 - p), shape = sigma) # this is wrong
  
  gamma_inv <- qgamma((1 - p), shape = sigma) 
  
  term <- 1 - exp(-gamma_inv / b)
  G <- 2 * (1 - term^(1/2))^(1/a) / (1 + (1 - term^(1/2))^(1/a))
  
  # Compute x from G(x)
  x <- ((1 / (1 - G)) - 1)^(1 / c)
  
  return(x)
}

### True parameters
sigma1 = 0.4; b1 = 0.4; a1 = 0.9; c1 = 2.5;
init_cond <- c(0.45, 0.95, 2.55)  # b is fixed

# sigma1 = 0.4; b1 = 0.9; a1 = 1.2; c1 = 1.0;
# init_cond <- c(0.45, 1.21, 1.1)  # b is fixed

# sigma1 = 0.6; b1 = 0.6; a1 = 1.2; c1 = 2.1;
# init_cond <- c(0.65, 1.24, 2.4)  # b is fixed

# sigma1 = 1.1; b1 = 1.5; a1 = 0.4; c1 = 2.5;
# init_cond <- c(1.15, 0.45, 2.55)  # b is fixed


# init_cond <- c(0.45, 0.45, 0.95)  # c is fixed
# init_cond <- c(0.45, 0.45, 2.55)  # a is fixed
# init_cond <- c(0.45, 0.95, 2.55)  # b is fixed

ns <- c(25, 50, 100, 200, 400, 800)

sigmas <- c()
bs <- c()
as <- c()
cs <- c()
error1 <- c(); error2 <- c(); error3 <- c();
mse1 <- c();   mse2 <- c();   mse3 <- c();

NN <- 900
for (j in 1:length(ns)) {
  n <- ns[j]
  columns <- c('param1', 'param2', 'param3')
  error <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(error) <- columns
  
  for (k in 1:NN) {
    Fx <- runif(n)
    x <- sapply(Fx, function(p) quantile_function(p, sigma1, b1, a1, c1))  # Generate random samples in one step
    ins <- (2 * (1:n) - 1) / (2 * n)
    ls_ins <- (1:n) / (n + 1)
    wls_ins <- (n + 1)^2 * (n + 2) / ((1:n) * (n - (1:n) + 1))
    iMPS <- 1/(n+1)
    i21 <- 2 * (1:n) - 1
    
    x <- sort(x)
    x_rev <- sort(x, decreasing = TRUE)
    
    #est_nlminb <- nlminb(init_cond, MLE_fun_param1B, lower = 0.001, x = x)
    
    est_nlminb <- nlminb(init_cond, LS_fun_param1B, lower = 0.001, x = x, ls_ins = ls_ins)
    
    # est_nlminb <- nlminb(init_cond, WLS_fun_param1B, lower = 0.001, x = x, ls_ins = ls_ins, wls_ins = wls_ins)
    
    # est_nlminb <- nlminb(init_cond, MPS_fun_param1B, lower = 0.001, x = x, n = n, iMPS = iMPS)
    
    # est_nlminb <- nlminb(init_cond, CVM_fun_param1B, lower = 0.001, x = x, n = n, ins = ins)
    
    # est_nlminb <- nlminb(init_cond, ADE_fun_param1B, lower = 0.001, x = x, x_rev = x_rev, i21 = i21)
    
    error[k, ] <- est_nlminb$par
    
    if (est_nlminb$convergence == 0) {  # Check for successful convergence
      error[k, ] <- est_nlminb$par
    } else {
      warning("Optimization did not converge on iteration ", k)
      error[k, ] <- NA
    }
  }
  
  # # fix c
  # error1 <- c(error1, mean(error$param1) - sigma1)
  # error2 <- c(error2, mean(error$param2) - b1)
  # error3 <- c(error3, mean(error$param3) - a1)
  # mse1 <- c(mse1, sum((error$param1 - sigma1)^2) / NN)
  # mse2 <- c(mse2, sum((error$param2 - b1)^2) / NN)
  # mse3 <- c(mse3, sum((error$param3 - a1)^2) / NN)
  
  # #fix a
  # error1 <- c(error1, mean(error$param1) - sigma1)
  # error2 <- c(error2, mean(error$param2) - b1)
  # error3 <- c(error3, mean(error$param3) - c1)
  # mse1 <- c(mse1, sum((error$param1 - sigma1)^2) / NN)
  # mse2 <- c(mse2, sum((error$param2 - b1)^2) / NN)
  # mse3 <- c(mse3, sum((error$param3 - c1)^2) / NN)
  
  # fix b
  error1 <- c(error1, mean(error$param1, na.rm = TRUE) - sigma1)
  error2 <- c(error2, mean(error$param2, na.rm = TRUE) - a1)
  error3 <- c(error3, mean(error$param3, na.rm = TRUE) - c1)
  
  mse1 <- c(mse1, mean((error$param1 - sigma1)^2, na.rm = TRUE))
  mse2 <- c(mse2, mean((error$param2 - a1)^2, na.rm = TRUE))
  mse3 <- c(mse3, mean((error$param3 - c1)^2, na.rm = TRUE))
  
  # # fix signma
  # error1 <- c(error1, mean(error$param1) - b1)
  # error2 <- c(error2, mean(error$param2) - a1)
  # error3 <- c(error3, mean(error$param3) - c1)
  # mse1 <- c(mse1, sum((error$param1 - b1)^2) / NN)
  # mse2 <- c(mse2, sum((error$param2 - a1)^2) / NN)
  # mse3 <- c(mse3, sum((error$param3 - c1)^2) / NN)

  rmse1 <- sqrt(mse1)
  rmse2 <- sqrt(mse2)
  rmse3 <- sqrt(mse3)
  print(j)
  print(k)
}

error1
error2
error3

rmse1
rmse2
rmse3
