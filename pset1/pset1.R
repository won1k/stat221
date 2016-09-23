### Problem Set 1
### STAT 221
### Won Lee

library(mvtnorm)

# Problem 1

# 1.4
dhypercube = function(u, mu, alpha, beta) {
  d = length(u)
  Sigma = (alpha+beta) * diag(d) - matrix(beta, ncol = d, nrow = d)
  x = log(u/(1-u))
  return(dmvnorm(x, mean = mu, sigma = Sigma)/(prod(u)*prod(1-u)))
}

rhypercube = function(n, mu, alpha, beta) {
  d = length(mu)
  Sigma = (alpha+beta) * diag(d) - matrix(beta, ncol = d, nrow = d)
  x = rmvnorm(n, mean = mu, sigma = Sigma)
  return(1/(1+exp(-x)))
}

hypercube.mean = function(mu, alpha, beta) {
  mu.t = 1/(1+exp(-mu))
  return(mu.t + 0.5*alpha*mu.t*(1-mu.t)*(1-2*mu.t))
}

hypercube.var = function(mu, alpha, beta) {
  d = length(mu)
  mu.t = 1/(1+exp(-mu))
  Sigma = matrix(0, ncol = d, nrow = d)
  for (i in 1:d) {
    for (j in 1:d) {
      if (i == j) {
        Sigma[i,i] = mu.t[i]^2*(1-mu.t[i])^2*(alpha+0.25*alpha^2*(1-2*mu.t[i])^2)
      } else {
        Sigma[i,j] = mu.t[i]*mu.t[j]*(0.5*alpha*(1-mu.t[j])*(1-2*mu.t[j])) + 
          0.5*alpha*(1-mu.t[i])*(1-2*mu.t[i]) - 0.25*alpha^2*(1-mu.t[i])*(1-mu.t[j])*(1-2*mu.t[i])*(1-2*mu.t[j]) -
          beta*(1-mu.t[i])*(1-mu.t[j])
      }
    }
  }
  return(Sigma)
}

hypercube.mle = function(U) {
  N = nrow(U)
  d = ncol(U)
  X = log(U/(1-U))
  mu.hat = colSums(X) / N
  Sigma.hat = 0
  for (n in 1:N) {
    Sigma.hat = Sigma.hat + (unlist(X[n,]) - mu.hat) %o% (unlist(X[n,]) - mu.hat)
  }
  Sigma.hat = Sigma.hat / N
  alpha.hat = 0
  beta.hat = 0
  for (i in 1:d) {
    alpha.hat = alpha.hat + Sigma.hat[i,i]
    for (j in 1:d) {
      if (i != j) {
        beta.hat = beta.hat - Sigma.hat[i,j]
      }
    }
  }
  alpha.hat = alpha.hat / d
  beta.hat = beta.hat / (d*(d-1))
  return(list("mu.hat" = mu.hat, "alpha.hat" = alpha.hat, "beta.hat" = beta.hat))
}

# 1.5
n = 250
mu1 = c(5,5)
alpha1 = 10
beta1 = 1
U1 = rhypercube(n, mu1, alpha1, beta1)
hypercube.mle(U1)

mu2 = c(1,1,10)
alpha2 = 1
beta2 = 0.1
U2 = rhypercube(n, mu2, alpha2, beta2)
hypercube.mle(U2) # (0.92, 1.03, 9.96) & 0.924 & 0.064

mu3 = c(5,-5,10,-10)
alpha3 = 50
beta3 = 5
U3 = rhypercube(n, mu3, alpha3, beta3)
hypercube.mle(U3) # (4.52, -5.34, 9.92, -9.54) & 48.568 & 3.163

# 1.6
data = read.table("dataHypercube3D.txt", header = TRUE)
mle.list = hypercube.mle(data) # (0.379, 0.366) & 2.915 & 0.75

# 1.7 (Simulations)
nsims = 1000
N = nrow(data)
mu.hats = data.frame(matrix(rep(0, nsims*d), ncol = d))
alpha.hats = rep(0, nsims)
beta.hats = rep(0, nsims)
for (t in 1:nsims) {
  idx = sample(1:N, N, replace = TRUE)
  U = data[idx,]
  mle = hypercube.mle(U)
  mu.hats[t,] = mle$mu.hat
  alpha.hats[t] = mle$alpha.hat
  beta.hats[t] = mle$beta.hat
}
par(mfrow = c(2,2))
hist(mu.hats[,1], xlab = "mu_1", main = "", breaks = 20)
hist(mu.hats[,2], xlab = "mu_2", main = "", breaks = 20)
hist(alpha.hats, xlab = "alpha", main = "", breaks = 20)
hist(beta.hats, xlab = "beta", main = "", breaks = 20)
dev.off()

# 1.8

beta.mom = function(U) {
  d = ncol(U)
  U.vec = c(U[,1])
  for (i in 2:d) {
    U.vec = c(U.vec, U[,i])
  }
  x.bar = mean(U.vec)
  s.2 = var(U.vec)
  a.hat = x.bar*(x.bar*(1-x.bar)/s.2 - 1)
  b.hat = (1 - x.bar)*(x.bar*(1-x.bar)/s.2 - 1)
  return(list("a.hat" = a.hat, "b.hat"= b.hat))
}

mom.list = beta.mom(data) # 1.085 & 0.849
beta.samples = matrix(rbeta(200, mom.list$a.hat, mom.list$b.hat), ncol = 2)
hyper.samples = matrix(rhypercube(200, mle.list$mu, mle.list$alpha, mle.list$beta), ncol = 2)

par(mfrow = c(2,1))
hist(data[,1], breaks = 10, main = "", xlab = "u_1")
x.points = 1:10000 / 10000
lines(x.points, dbeta(x.points, mom.list$a.hat, mom.list$b.hat) * 10)
hist(data[,2], breaks = 10, main = "", xlab = "u_2")
x.points = 1:10000 / 10000
lines(x.points, dbeta(x.points, mom.list$a.hat, mom.list$b.hat) * 10)
dev.off()

library(goftest)
ks.test(data[,1], beta.samples) # 0.06 & 0.97
ks.test(data[,2], beta.samples) # 0.075 & 0.8475
ks.test(data[,1], hyper.samples[,1]) # 0.085 & 0.721
ks.test(data[,2], hyper.samples[,2]) # 0.055 & 0.988
