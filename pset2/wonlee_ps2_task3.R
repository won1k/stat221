#!/usr/bin/env Rscript

# Problem Set 1
# Task 3
# Won Lee

library(mvtnorm)
library(matrixStats)
source("wonlee_ps2_functions.R")
source("possonLogN_MCMC.R")

args = commandArgs(trailingOnly=TRUE)

# Hyperparameters
mus = c(1.6, 2.5, 5.2, 4.9)
sigmas = c(0.7, 1.3, 1.3, 1.6)
J = 1000
N = 2
w = rep(1, J)
ndat = 25 # Number of data sets per simulation (1 theta per simulation)
ndraws = 1e4

# Simulate data
curr_sim = args[1]
mu = mus[curr_sim]
sigma = sigmas[curr_sim]
log.theta = rnorm(J, mu, sigma)

# Run MCMC to obtain posterior
post.means = matrix(rep(0, ndat*J), nrow = J)
post.std = matrix(rep(0, ndat*J), nrow = J)
lb.95 = matrix(rep(0, ndat*J), nrow = J)
ub.95 = matrix(rep(0, ndat*J), nrow = J)
start = proc.time()
for (i in 1:ndat) {
  Y = matrix(c(rpois(J, exp(log.theta)), rpois(J, exp(log.theta))), ncol = N)
  results = poisson.log.mcmc(Y, w)$logTheta
  post.means[,i] = rowMeans(results)
  post.std[,i] = rowSds(results)
  lb.95[,i] = apply(results, 2, function(x) return(sort(x)[floor(0.025*ndraws)]))
  ub.95[,i] = apply(results, 2, function(x) return(sort(x)[floor(0.975*ndraws)]))
}
time = (proc.time() - start)[3]
print(time)



