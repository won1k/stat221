#!/usr/bin/env Rscript

# Problem Set 1
# Task 3
# Won Lee

library(mvtnorm)
library(matrixStats)
source("wonlee_ps2_functions.R")
source("poissonLogN_MCMC.R")

args = commandArgs(trailingOnly=TRUE)

# Hyperparameters
mus = c(1.6, 2.5, 5.2, 4.9)
sigmas = c(0.7, 1.3, 1.3, 1.6)
J = 1000
N = 2
w = rep(1, J)
nsims = 11 # Number of theta draws per simulation
ndat = 25 # Number of data sets per theta
ndraws = 1e4
curr_sim = strtoi(args[1])
mu = mus[curr_sim]
sigma = sigmas[curr_sim]

# Store simulation values
log.thetas = matrix(rep(0, J*nsims), nrow = J)
coverage.probs = matrix(rep(0, J*nsims), nrow = J)
runtimes = rep(0, nsims)

for (t in 1:nsims) {
  # Simulate theta
  log.thetas[,t] = rnorm(J, mu, sigma)
  # Run MCMC to obtain posterior
  lb.95 = matrix(rep(0, ndat*J), nrow = J)
  ub.95 = matrix(rep(0, ndat*J), nrow = J)
  cov = rep(0, J)
  start = proc.time()
  for (i in 1:ndat) {
    print(i)
    Y = matrix(c(rpois(J, exp(log.thetas[,t])), rpois(J, exp(log.thetas[,t]))), ncol = N)
    results = poisson.logn.mcmc(Y, w)$logTheta
    lb.95[,i] = apply(results, 1, function(x) return(sort(x)[floor(0.025*ndraws)]))
    ub.95[,i] = apply(results, 1, function(x) return(sort(x)[floor(0.975*ndraws)]))
    for (j in 1:J) {
      if (log.thetas[j,t] > lb.95[j,i]) {
        if (log.thetas[j,t] < ub.95[j,i]) {
          cov[j] = cov[j] + 1
        }
      }
    }
  }
  coverage.probs[,t] = cov/ndat
  runtimes[t] = (proc.time() - start)[3]
}
# Vectorize log.thetas, coverage.probs
output = data.frame(as.vector(log.thetas), as.vector(coverage.probs))
write.table(output, file = paste("wonlee_ps2_task3_par", curr_sim, "_theta.dat", sep = ""),
            row.names = FALSE, sep="\t", quote = FALSE)
write.table(runtimes, file =  paste("wonlee_ps2_task3_par", curr_sim, "_runtimes.dat", sep = ""),
            row.names = FALSE)

