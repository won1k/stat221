#!/usr/bin/env Rscript
# Problem Set 4
# Won Lee

# Load data/MCMC
source("wonlee_mcmc.R")
library(coda)
data = read.csv("1router_allcount.dat", header = TRUE)
A = matrix(c(rep(1,4),rep(0,12),rep(0,4),rep(1,4),rep(0,8),rep(0,8),rep(1,4),rep(0,4),rep(0,12),rep(1,4),rep(c(1,0,0,0),4),rep(c(0,1,0,0),4),rep(c(0,0,1,0),4)), ncol = 16, byrow = TRUE)
Y = data[117:123,"value"]
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
X1 = matrix(rep(0, 10*r*(iter-burnin)), ncol = iter-burnin)
X2 = matrix(rep(0, 10*(c-r)*(iter-burnin)), ncol = iter-burnin)
lambdas = matrix(rep(0, 10*c*(iter-burnin)), ncol = iter-burnin)
for (i in 1:10) {
  mcmc_results = network_mcmc(Y, A, prior, iter, burnin, TRUE, 5)
  X1[(r*(i-1)+1):(r*i),] = mcmc_results$X1
  X2[((c-r)*(i-1)+1):((c-r)*i),] = mcmc_results$X2
  lambdas[(c*(i-1)+1):(c*i),] = mcmc_results$lambda
}

# Plots (Diagnostics)
# ACF
for (i in 1:(c-r)){
  png(paste("wonlee_1router_t5_informative_acf_",i,".png",sep=""))
  x = acf(X2[i,], plot = FALSE)$lag
  mcmc_acf = acf(X2[i,], plot = FALSE)$acf
  plot(x, mcmc_acf, main = paste("X2[",i,"]", sep = ""), type = "l", xlab = "Lag", ylab = "ACF",
       ylim = c(0.8,1))
  for (j in 2:10) {
    mcmc_acf = acf(X2[(c-r)*(j-1)+i,], plot = FALSE)$acf
    lines(x, mcmc_acf, col = j)
  }
  dev.off()
}

# ESS
#png(paste("wonlee_1router_t5_informative_ess.png",sep=""))
#heights = rep(0,c)
#for (i in 1:c) {
#  if (i <= r) {
#    x = X1[i,]
#    for (j in 2:10) {
#      x = c(x, X1[(j-1)*r+i,])
#    }
#  } else {
#    x = X2[i-r,]
#    for (j in 2:10) {
#      x = c(x, X2[(j-1)*(c-r)+i-r,])
#    }
#  }
#  heights[i] = effectiveSize(x)
#}
#barplot(heights, names.arg = 1:c)
#dev.off()

# Plots (Marginals)
png(paste("wonlee_1router_t5_informative_marginals.png",sep=""))
par(mfrow = c(4,4))
for (i in 1:c) {
  if (i <= r) {
    x = X1[i,]
    for (j in 2:10) {
      x = c(x, X1[(j-1)*r+i,])
    }
  } else {
    x = X2[i-r,]
    for (j in 2:10) {
      x = c(x, X2[(j-1)*(c-r)+i-r,])
    }
  }
  hist(x, main = paste("X",i,sep=""))
}
dev.off()

png(paste("wonlee_1router_t5_informative_densities.png",sep=""))
par(mfrow = c(4,4))
for (i in 1:c) {
  lamb = lambdas[i,]
  for (j in 2:10) {
    lamb = c(x, lambdas[(j-1)*c+i,])
  }
  hist(lamb, main = "", ylab = "", xlab = paste("Lambda", i))
}
dev.off()

# Save data
write.table(X1, file = "wonlee_1router_t5_informative_X1.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(X2, file = "wonlee_1router_t5_informative_X2.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(lambdas, file = "wonlee_1router_t5_informative_lambdas.csv", sep = ",", row.names = FALSE, col.names = FALSE)


