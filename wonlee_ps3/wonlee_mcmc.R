library(MASS)

accept = function(next.N, next.theta, N, theta, y) {
  if (next.N < max(y)) {
    return(0)
  } else{
    num = N * prod(choose(next.N, y) * next.theta^y * (1-next.theta)^(next.N-y))
    den = next.N * prod(choose(N, y) * theta^y * (1-theta)^(N-y))
    return(min(1, num/den))
  }
}

mcmc = function(y, burnin, iterations, chains) {
  T = iterations - burnin
  samples.N = matrix(0, nrow = chains, ncol = T)
  samples.theta = matrix(0, nrow = chains, ncol = T)
  for (chain in 1:chains) {
    chain.samples.N = rep(0, iterations)
    chain.samples.N[1] = max(y) + rpois(1, 10)
    chain.samples.theta = rep(0, iterations)
    chain.samples.theta[1] = runif(1,0,1)
    for (i in 2:iterations) {
      next.N = max(0, chain.samples.N[i-1] + (2*rbinom(1,1,0.5)-1) * rpois(1,1))
      next.theta = min(1,max(0,rnorm(1,chain.samples.theta[i-1], 0.05)))
      alpha = accept(next.N, next.theta, chain.samples.N[i-1], chain.samples.theta[i-1], y)
      if (runif(1,0,1) < alpha) {
        chain.samples.N[i] = next.N
        chain.samples.theta[i] = next.theta
      } else {
        chain.samples.N[i] = chain.samples.N[i-1]
        chain.samples.theta[i] = chain.samples.theta[i-1]
      }
    }
    samples.N[chain,] = chain.samples.N[burnin+1:T]
    samples.theta[chain,] = chain.samples.theta[burnin+1:T]
  }
  return(list(N = samples.N, theta = samples.theta))
}

# Load data
waterbuck = read.csv("waterbuck.txt")[,1]
impala = read.csv("impala.txt")[,1]

# parameters
chains = 10
burnin = 5000
iterations = 20000

### Impala
# mcmc
samples.impala = mcmc(impala, burnin, iterations, chains)
# diagnostics
uniqs = rep(0, chains)
for (i in 1:chains) {
  uniqs[i] = length(unique(samples.impala$N[i,]))
}
#mean(uniqs) # ~1000 per 5000 samples, so around 0.2
#plot(1:5000, x$N[2,], type = "l")
#plot(1:5000, x$theta[1,], type = "l")

# scatterplots
for (i in 1:chains) {
  png(paste("wonlee_mcmc_impala_", i, ".png", sep = ""))
  plot(samples.impala$N[i,], samples.impala$theta[i,], xlab = "N", ylab = "theta", pch  = 20, 
       cex = 0.1, main = paste("Impala Data, Chain", i, sep = " "), col = rgb(0,0,0,0.2))
  den = kde2d(samples.impala$N[i,], samples.impala$theta[i,])
  contour(den$x, den$y, den$z, add = TRUE, col = 2, lty = 1, lwd = 1.2) 
  dev.off()
}


### Waterbuck
# mcmc
samples.waterbuck = mcmc(waterbuck, burnin, iterations, chains)
# scatterplots
for (i in 1:chains) {
  png(paste("wonlee_mcmc_waterbuck_", i, ".png", sep = ""))
  plot(samples.waterbuck$N[i,], samples.waterbuck$theta[i,], xlab = "N", ylab = "theta", pch  = 20, 
       cex = 0.1, main = paste("Waterbuck Data, Chain", i, sep = " "), col = rgb(0,0,0,0.2))
  den = kde2d(samples.waterbuck$N[i,], samples.waterbuck$theta[i,])
  contour(den$x, den$y, den$z, add = TRUE, col = 2, lty = 1, lwd = 1.2) 
  dev.off()
}


### Prob N > 100
# Impala
probs.impala = rep(0, chains)
for (i in 1:chains) {
  probs.impala[i] = sum(samples.impala$N[i,] > 100) / length(samples.impala$N[i,])
}
mean(probs.impala) # 0.169
sqrt(var(probs.impala)) # 0.305
# Waterbuck
probs.waterbuck = rep(0, chains)
for (i in 1:chains) {
  probs.waterbuck[i] = sum(samples.waterbuck$N[i,] > 100) / length(samples.waterbuck$N[i,])
}
mean(probs.waterbuck) # 0.939
sqrt(var(probs.waterbuck)) # 0.102




