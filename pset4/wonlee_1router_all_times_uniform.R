#!/usr/bin/env Rscript
# Problem Set 4
# Won Lee

# Command line args
args = commandArgs(trailingOnly = TRUE)
t = as.numeric(args[1])

# Load data/MCMC
source("wonlee_mcmc.R")
library(coda)
data = read.csv("1router_allcount.dat", header = TRUE)
A = matrix(c(rep(1,4),rep(0,12),rep(0,4),rep(1,4),rep(0,8),rep(0,8),rep(1,4),rep(0,4),rep(0,12),rep(1,4),rep(c(1,0,0,0),4),rep(c(0,1,0,0),4),rep(c(0,0,1,0),4)), ncol = 16, byrow = TRUE)
Y = data[(25*(t-1)+17):(25*(t-1)+23),"value"]
r = qr(A)$rank
c = ncol(A)
idx1 = qr(A)$pivot[0:r]
A1 = A[,idx1]
idx2 = qr(A)$pivot[(r+1):c]
A2 = A[,idx2]

# Run MCMC (10 chains)
prior = "uniform"
iter = 1.2e5
burnin = 2e4
X1 = matrix(rep(0, 10*r*(iter-burnin)), ncol = iter-burnin)
X2 = matrix(rep(0, 10*(c-r)*(iter-burnin)), ncol = iter-burnin)
lambdas = matrix(rep(0, 10*c*(iter-burnin)), ncol = iter-burnin)
for (i in 1:10) {
  mcmc_results = network_mcmc(Y, A, prior, iter, burnin, TRUE, t)
  X1[(r*(i-1)+1):(r*i),] = mcmc_results$X1
  X2[((c-r)*(i-1)+1):((c-r)*i),] = mcmc_results$X2
  lambdas[(c*(i-1)+1):(c*i),] = mcmc_results$lambda
}

# Make quantiles

# Save results
write.table(X1, file = "wonlee_1router_t5_uniform_X1.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(X2, file = "wonlee_1router_t5_uniform_X2.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(lambdas, file = "wonlee_1router_t5_uniform_lambdas.csv", sep = ",", row.names = FALSE, col.names = FALSE)



