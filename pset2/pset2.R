### Problem Set 1
### STAT 221
### Won Lee

library(mvtnorm)
source("wonlee_ps2_functions.R")

# Task 1
n = 10
y.bar = 2
mu = 1
sigma.2 = 5
w = 0.01
log.theta.points = c(-10000:10000/1000)

cond.post = function(logtheta) {
  exp.term = exp(-0.5*(1/sigma.2)*(logtheta-mu)^2-n*w*exp(logtheta))
  return(exp.term*(w*exp(logtheta))^(n*y.bar))
}

post.points = cond.post(log.theta.points)
plot(log.theta.points, post.points, type= "l")

# Task 2

# Task 3
poisson.logn.mcmc()



