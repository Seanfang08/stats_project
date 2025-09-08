################ List of parameters with fixed values
# b <- params[1]
# a <- params[2]
# c <- params[3]
# sigma <- 0.6  # Assuming sigma is fixed
# 
# sigma <- params[1]
# a <- params[2]
# c <- params[3]
# b <- 0.6  # Assuming b is fixed
# 
# sigma <- params[1]
# b <- params[2]
# c <- params[3]
# a <- 1.2  # Assuming a is fixed
# 
# sigma <- params[1]
# b <- params[2]
# a <- params[3]
# c <- 2.1  # Assuming c is fixed

############### ADE fix sigma SSSSSSSSSSSSSSSS ###########################
ADE_fun_param1 <- function(params, x, x_rev, i21) {
  b <- params[1]
  a <- params[2]
  c <- params[3]
  sigma <- 0.6  # Assuming sigma is fixed
  
  F_x <- 1 - pgamma(-log((1 - (1 - ((1 - (1 + x^c)^(-1)) / (1 + (1 + x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  # ADE <- -(sum(i21 * (log(F_x(x)) + log(1 - F_x(x_rev)))) WRONG
  
  # ADE <- -(sum(i21 * (log(F_x) + log(1 - F_x[length(x):1]))))
  
  # ADE <- -(sum(i21 * (log(F_x) + log(1 - F_x[order(x_rev)]))))
  
  ADE <- - sum(i21 * (log(F_x) + log(1 - rev(F_x))))
  return(ADE)
}

############### ADE fix b. BBBBBBBBBBBB########################
ADE_fun_param1B <- function(params, x, x_rev, i21) {
  sigma <- params[1]
  a <- params[2]
  c <- params[3]
  b <- 0.6  # Assuming b is fixed
  
  F_x <- 1 - pgamma(-log((1 - (1 - ((1 - (1 + x^c)^(-1)) / (1 + (1 + x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  # ADE <- -(sum(i21 * (log(F_x(x)) + log(1 - F_x(x_rev)))) WRONG
  
  # ADE <- -(sum(i21 * (log(F_x) + log(1 - F_x[length(x):1]))))
  
  # ADE <- -(sum(i21 * (log(F_x) + log(1 - F_x[order(x_rev)]))))
  
  ADE <- -sum(i21 * (log(F_x) + log(1 - rev(F_x))))
  return(ADE)
}

############### ADE fix a AAAAAAAAAAAAAAAAAA##########################
ADE_fun_param1A <- function(params, x, x_rev, i21) {
  sigma <- params[1]
  b <- params[2]
  c <- params[3]
  a <- 1.2  # Assuming a is fixed
  
  F_x <- 1 - pgamma(-log((1 - (1 - ((1 - (1 + x^c)^(-1)) / (1 + (1 + x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  # ADE <- -(sum(i21 * (log(F_x(x)) + log(1 - F_x(x_rev)))) WRONG
  
  # ADE <- -(sum(i21 * (log(F_x) + log(1 - F_x[length(x):1]))))
  
  # ADE <- -(sum(i21 * (log(F_x) + log(1 - F_x[order(x_rev)]))))
  
  ADE <- -sum(i21 * (log(F_x) + log(1 - rev(F_x))))
  return(ADE)
}

############### ADE fix c CCCCCCCCCCCCCCCCC#####################
ADE_fun_param1C <- function(params, x, x_rev, i21) {
  sigma <- params[1]
  b <- params[2]
  a <- params[3]
  c <- 2.1  # Assuming c is fixed
  
  F_x <- 1 - pgamma(-log((1 - (1 - ((1 - (1 + x^c)^(-1)) / (1 + (1 + x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  # ADE <- -(sum(i21 * (log(F_x(x)) + log(1 - F_x(x_rev)))) WRONG
  
  # ADE <- -(sum(i21 * (log(F_x) + log(1 - F_x[length(x):1]))))
  
  # ADE <- -(sum(i21 * (log(F_x) + log(1 - F_x[order(x_rev)]))))
  
  ADE <- -sum(i21 * (log(F_x) + log(1 - rev(F_x))))
  return(ADE)
}

############### MPS fix sigma SSSSSSSSSSSSSSSSSSSS #######################
MPS_fun_param1 <- function(params, x, n, iMPS) {
  b <- params[1]
  a <- params[2]
  c <- params[3]
  sigma <- 0.6  # Assuming sigma is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  x1 <- F_x[1:(n-1)]
  x2 <- F_x[2:n]
  
  MPS <- iMPS * (log(F_x[1]) + log(1 - F_x[n]) + sum(log(x2 - x1)))
  MPS <- -MPS
  return(MPS)
}

############### MPS fix b BBBBBBBBBBBBBBB ################################
MPS_fun_param1B <- function(params, x, n, iMPS) {
  sigma <- params[1]
  a <- params[2]
  c <- params[3]
  b <- 0.6  # Assuming b is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  x1 <- F_x[1:(n-1)]
  x2 <- F_x[2:n]
  
  MPS <- iMPS * (log(F_x[1]) + log(1 - F_x[n]) + sum(log(x2 - x1)))
  MPS <- -MPS
  return(MPS)
}

############### MPS fix a AAAAAAAAAAA###################
MPS_fun_param1A <- function(params, x, n, iMPS) {
  sigma <- params[1]
  b <- params[2]
  c <- params[3]
  a <- 1.2  # Assuming a is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  x1 <- F_x[1:(n-1)]
  x2 <- F_x[2:n]
  
  MPS <- iMPS * (log(F_x[1]) + log(1 - F_x[n]) + sum(log(x2 - x1)))
  MPS <- -MPS
  return(MPS)
}

############### MPS fix c CCCCCCCCCCCCCCCCCCCC##############################
MPS_fun_param1C <- function(params, x, n, iMPS) {
  sigma <- params[1]
  b <- params[2]
  a <- params[3]
  c <- 2.1  # Assuming c is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  x1 <- F_x[1:(n-1)]
  x2 <- F_x[2:n]
  
  MPS <- iMPS * (log(F_x[1]) + log(1 - F_x[n]) + sum(log(x2 - x1)))
  MPS <- -MPS
  return(MPS)
}


############### WLS fix sigma SSSSSSSSSSS################################
WLS_fun_param1 <- function(params, x, ls_ins, wls_ins) {
  b <- params[1]
  a <- params[2]
  c <- params[3]
  sigma <- 0.6  # Assuming sigma is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  WLS <- sum(wls_ins*(F_x - ls_ins)**2)
  
  return(WLS)
}

############### WLS fix b BBBBBBBBBBBB########################
WLS_fun_param1B <- function(params, x, ls_ins, wls_ins) {
  sigma <- params[1]
  a <- params[2]
  c <- params[3]
  b <- 0.6  # Assuming b is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  WLS <- sum(wls_ins*(F_x - ls_ins)**2)
  
  return(WLS)
}

############### WLS fix aAAAAAAAAAAA#######################
WLS_fun_param1A <- function(params, x, ls_ins, wls_ins) {
  sigma <- params[1]
  b <- params[2]
  c <- params[3]
  a <- 1.2  # Assuming a is fixed 
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  WLS <- sum(wls_ins*(F_x - ls_ins)**2)
  
  return(WLS)
}

############### WLS fix c CC##############################
WLS_fun_param1C <- function(params, x, ls_ins, wls_ins) {
  sigma <- params[1]
  b <- params[2]
  a <- params[3]
  c <- 2.1  # Assuming c is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  WLS <- sum(wls_ins*(F_x - ls_ins)**2)
  
  return(WLS)
}
############### LS fix sigma SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
LS_fun_param1 <- function(params, x, ls_ins) {
  b <- params[1]
  a <- params[2]
  c <- params[3]
  sigma <- 0.6  # Assuming sigma is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  LS <- sum((F_x-ls_ins)**2)
  
  return(LS)
}

############### LS fix b BBBBBBBBBBBBBBBBBB
LS_fun_param1B <- function(params, x, ls_ins) {
  sigma <- params[1]
  a <- params[2]
  c <- params[3]
  b <- 0.6  # Assuming b is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  LS <- sum((F_x-ls_ins)**2)
  
  return(LS)
}

############### LS fix a AAAAAAAAAAAAAAAAAAA
LS_fun_param1A <- function(params, x, ls_ins) {
  sigma <- params[1]
  b <- params[2]
  c <- params[3]
  a <- 1.2  # Assuming a is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  LS <- sum((F_x - ls_ins)^2)
  
  return(LS)
}

############### LS fix c CCCCCCCCCCCCCCCCCCCCCC
LS_fun_param1C <- function(params, x, ls_ins) {
  sigma <- params[1]
  b <- params[2]
  a <- params[3]
  c <- 2.1  # Assuming c is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  LS <- sum((F_x-ls_ins)**2)
  
  return(LS)
}

################## CVM sigma SSSSSSSSSSSSSSSSSS
CVM_fun_param1 <- function(params, x, n, ins) { 
  b <- params[1]
  a <- params[2]
  c <- params[3]
  sigma <- 0.6  # Assuming sigma is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  CVM <- (1/n) * sum((F_x - ins)^2)
  
  return(CVM)
}

################## CVM b BBBBBBBBBBBBBBBB
CVM_fun_param1B <- function(params, x, n, ins) { 
  sigma <- params[1]
  a <- params[2]
  c <- params[3]
  b <- 0.6  # Assuming b is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  CVM <- (1/n) * sum((F_x - ins)^2)
  
  return(CVM)
}

################## CVM AAAAAAAAAAAAAAAAAAAA
CVM_fun_param1A <- function(params, x, n, ins) { 
  sigma <- params[1]
  b <- params[2]
  c <- params[3]
  a <- 1.2  # Assuming a is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  CVM <- (1/n) * sum((F_x - ins)^2)
  
  return(CVM)
}

################## CVM c CCCCCCCCCCCCCCCCCCC
CVM_fun_param1C <- function(params, x, n, ins) { 
  sigma <- params[1]
  b <- params[2]
  a <- params[3]
  c <- 2.1  # Assuming c is fixed
  
  F_x <- 1 - pgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  
  CVM <- (1/n) * sum((F_x - ins)^2)
  
  return(CVM)
}
################### MLE Functionfix sigma SSSSSSSSSSSSSSS
MLE_fun_param1 <- function(params, x){
  b <- params[1]
  a <- params[2]
  c <- params[3]
  sigma <- 0.6  # Assuming sigma is fixed
  
  f_x <- dgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  term1 <- 4 * b * a * c *(x^(c-1)) * (x^c)^(a-1) *(1 - (x^c/(2+x^c))^a) /(2+x^c)^(a+1)
  term2 <- 1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2
  
  f_x <- f_x * term1/term2
  
  # Negative log-likelihood function
  LnL <- - sum(log(f_x))
  return(LnL)
}

################### MLE Function fix b BBBBBBBBBBBBBBBBBBB
MLE_fun_param1B <- function(params, x){
  sigma <- params[1]
  a <- params[2]
  c <- params[3]
  b <- 0.6  # Assuming b is fixed
  
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
  
  # Negative log-likelihood function
  LnL <- -sum(log(f_x))
  return(LnL)
}

################## MLE Function fix a AAAAAAAAAAAAAAAAAAAAA
MLE_fun_param1A <- function(params, x){
  sigma <- params[1]
  b <- params[2]
  c <- params[3]
  a <- 1.2  # Assuming a is fixed
  
  f_x <- dgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  term1 <- 4 * b * a * c *(x^(c-1)) * (x^c)^(a-1) *(1 - (x^c/(2+x^c))^a) /(2+x^c)^(a+1)
  term2 <- 1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2
  
  f_x <- f_x * term1/term2
  
  # Negative log-likelihood function
  LnL <- -sum(log(f_x))
  return(LnL)
}

################## MLE Function fix c CCCCCCCCCCCCCCCCCCC
MLE_fun_param1C <- function(params, x){
  sigma <- params[1]
  b <- params[2]
  a <- params[3]
  c <- 2.1  # Assuming c is fixed
  
  f_x <- dgamma( -log((1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2)^b), shape = sigma)
  term1 <- 4 * b * a * c *(x^(c-1)) * (x^c)^(a-1) *(1 - (x^c/(2+x^c))^a) /(2+x^c)^(a+1)
  term2 <- 1 - (1 - ((1 - (1+x^c)^(-1))/(1 + (1+ x^c)^(-1)))^a)^2
  
  f_x <- f_x * term1/term2
  
  # Negative log-likelihood function
  LnL <- -sum(log(f_x))
  return(LnL)
}
