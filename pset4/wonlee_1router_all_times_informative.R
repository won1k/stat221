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
prior = "informative"
iter = 1.2e5
burnin = 2e4
X1 = matrix(rep(0, 10*r*(iter-burnin)), ncol = 10*(iter-burnin))
X2 = matrix(rep(0, 10*(c-r)*(iter-burnin)), ncol = 10*(iter-burnin))
lambdas = matrix(rep(0, 10*c*(iter-burnin)), ncol = 10*(iter-burnin))
for (i in 1:10) {
  mcmc_results = network_mcmc(Y, A, prior, iter, burnin, TRUE, t)
  X1[,((i-1)*(iter-burnin)+1):(i*(iter-burnin))] = mcmc_results$X1
  X2[,((i-1)*(iter-burnin)+1):(i*(iter-burnin))] = mcmc_results$X2
  lambdas[,((i-1)*(iter-burnin)+1):(i*(iter-burnin))] = mcmc_results$lambda
}

# Make quantiles
q1 = data.frame(summary(mcmc(t(X1)))$quantiles[,1], summary(mcmc(t(X1)))$quantiles[,5])
colnames(q1) = c("2.5","97.5")
q2 = data.frame(summary(mcmc(t(X2)))$quantiles[,1], summary(mcmc(t(X2)))$quantiles[,5])
colnames(q2) = c("2.5","97.5")
l  = data.frame(summary(mcmc(t(lambdas)))$quantiles[,1], summary(mcmc(t(lambdas)))$quantiles[,5])
colnames(l) = c("2.5","97.5")


# Save results
write.table(q1, file = paste("wonlee_1router_t", t, "_informative_X1_quantiles.csv", sep=""), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(q2, file = paste("wonlee_1router_t", t, "_informative_X2_quantiles.csv", sep=""), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(l, file = paste("wonlee_1router_t", t, "_informative_lambda_quantiles.csv", sep=""), sep = ",", row.names = FALSE, col.names = FALSE)




