stirling = function(n) {
  return(sqrt(2*pi*n)*(n/exp(1))^n)
}

log.stirling = function(n) {
  return(n*log(n) - n + 0.5*log(2*pi*n))
}

unnormalized.waterbuck = function(N) {
  y = waterbuck
  n = length(y)
  if (N < max(y)) {
    return(0)
  } else {
    if (N < 250) {
      num = prod(choose(N,y))
      den = choose(n*N+1, sum(y))
      dens = num/den
      return(dens/(N*n*(N-mean(y)+1)))
    } else {
      dens = 1/N
      for (i in 1:n) {
        dens = dens * stirling(N)/(stirling(N-y[i]))
      }
      dens = dens * stirling(n*(N-mean(y)))
      return(dens/stirling(n*N+1))
    }
  }
}
  
unnormalized.impala = function(N) {
  y = impala
  n = length(y)
  if (N < max(y)) {
    return(0)
  } else {
    num = prod(choose(N,y))
    den = choose(n*N+1, sum(y))
    dens = num/den
    return(dens/(N*n*(N-mean(y)+1)))
  }
}

log.unnormalized = function(N) {
  y = waterbuck
  n = length(y)
  if (N <= max(y)) {
    return(-Inf)
  } else {
    log.dens = log(1/N)
    for (i in 1:n) {
      log.dens = log.dens + log.stirling(N) - log.stirling(N-y[i])
    }
    log.dens = log.dens + log.stirling(n*(N-mean(y))) - log.stirling(n*N+1)
    return(log.dens)
  }
}
 
# Data
waterbuck = read.csv("waterbuck.txt")[,1]
impala = read.csv("impala.txt")[,1]

# Diag (Impala)

x = 1:2000
dens.impala = rep(0, length(x))
for (i in 1:length(x)) { dens.impala[i] = unnormalized.impala(x[i]) }


# Diag (Waterbuck)
dens.waterbuck = rep(0, length(x))
for (i in 1:length(x)) { dens.waterbuck[i] = exp(log.unnormalized(x[i]))*exp(log.stirling(sum(y)) - sum(log.stirling(y))) }

# Normalizing
norm.impala = sum(dens.impala) # 
norm.waterbuck = sum(dens.waterbuck)
# Plots
png("wonlee_marginal_impala.png")
plot(x, dens.impala/norm.impala, xlab = "N", ylab = "Prob", type = "l")
dev.off()
png("wonlee_marginal_waterbuck.png")
plot(x, dens.waterbuck/norm.waterbuck, xlab = "N", ylab = "Prob", type = "l")
dev.off()

# Prob(N > 100)
sum(dens.impala[101:2000]/norm.impala) # 0.3279
sum(dens.waterbuck[101:2000]/norm.waterbuck) # Waterbuck: 0.9546
