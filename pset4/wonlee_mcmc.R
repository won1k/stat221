library(MASS)

MAX = 1e10
alpha = 3

log.stirling = function(n) {
  if (n > 0) {
    return(0.5*log(2*pi*n) + n*log(n) - n)
  } else {
    return(0)
  }
}

accept_prob = function(Y, A1, A2, currX2, propX2, lambdas, i) {
  currX1 = solve(A1) %*% (Y - A2 %*% currX2)
  propX1 = solve(A1) %*% (Y - A2 %*% propX2)
  r = ncol(A1)
  log.accept = dpois(round(currX2[i-r]), lambdas[i], log = TRUE) - dpois(propX2[i-r], lambdas[i], log = TRUE)
  #log.accept = 0
  log.accept = log.accept + log.stirling(currX2[i-r]) - log.stirling(propX2[i-r])
  log.accept = log.accept + (propX2[i-r] - currX2[i-r]) * log(lambdas[i])
  for (j in 1:r) {
    log.accept = log.accept + (propX1[j] - currX1[j]) * log(lambdas[j]) +
      log.stirling(currX1[j]) - log.stirling(propX1[j])
  }
  return(exp(log.accept))
}

network_mcmc = function(Y, A, prior, iter = 1.2e5, burnin = 2e4, verbose = FALSE) {
  # QR to find indices
  if (verbose) {
    print("Performing QR to find indices")
  }
  r = qr(A)$rank
  c = ncol(A)
  idx1 = qr(A)$pivot[0:r]
  A1 = A[,idx1]
  idx2 = qr(A)$pivot[(r+1):c]
  A2 = A[,idx2]
  
  # Initialize
  if (verbose) {
    print("Setting up matrices, initializing using prior")
  }
  lambs = matrix(rep(0, c*iter), ncol = iter)
  X1 = matrix(rep(0, r*iter), ncol = iter)
  X2 = matrix(rep(0, (c-r)*iter), ncol = iter)
  # Setup if informative
  if (prior == "informative") {
    a = 0.02
    X = read.csv("1router_allcount.dat", header = TRUE)[76:91,"value"]
    X = X[qr(A)$pivot]
    lambs[1:r,1] = rgamma(r, a*X[0:r], scale = a)
    lambs[(r+1):c,1] = rgamma(c-r, a*X[(r+1):c], scale = a)
  } else {
    lambs[,1] = runif(c, 0, max(Y)/16)
  }
  X2[,1] = read.csv("1router_allcount.dat", header = TRUE)[100+idx2,"value"] + 0.1
  X1[,1] = read.csv("1router_allcount.dat", header = TRUE)[100+idx1,"value"] + 1e-5
  
  # MCMC
  if (verbose) {
    print("Done! Starting MCMC iterations...")
    time = proc.time()
  }
  for (t in 2:iter) {
    if (verbose) {
      print(paste("MCMC Iteration: ", t))
    }
    # Draw lambdas
    if (prior == "uniform") {
      lambs[1:r,t] = rgamma(r, X1[,t-1]+1)
      lambs[(r+1):c,t] = rgamma(c-r, X2[,t-1]+1)
    } else {
      lambs[1:r,t] = rgamma(r, a*X[1:r] + X1[,t-1], scale = a)
      lambs[(r+1):c,t] = rgamma(c-r, a*X[(r+1):c] + X2[,t-1], scale = a)
    }
    #print(lambs[,t])
    
    # Draw X2 and compute X1
    X2[,t] = X2[,t-1]
    X1[,t] = X1[,t-1]
    
    for (i in (r+1):c) {
      #print(i)
      
      # Poisson proposal
      prop = rpois(1, lambs[i,t])
      #prop = max(0, runif(1, X2[i-r,t] - alpha, X2[i-r,t] + alpha))
      propX2 = X2[,t]
      propX2[i-r] = prop
      # Check feasibility
      propX1 = solve(A1) %*% (Y - A2 %*% propX2)
      #print(sum(propX1 >= 0))
      if (sum(propX1 >= 0) < length(X1[,t])) {
        X2[i-r,t] = X2[i-r,t-1]
      } else {
        # Accept-reject (i.e. reject if prob. too low)
        #print(accept_prob(Y, A1, A2, X2[,t], propX2, lambs[,t], i))
        #print(propX2)
        #print(X2[,t])
        #print(lambs[,t])
        if (runif(1) < accept_prob(Y, A1, A2, X2[,t], propX2, lambs[,t], i)) {
          X2[i-r,t] = prop
          X1[,t] = solve(A1) %*% (Y - A2 %*% X2[,t])
        }
      }
    }
  }
  if (verbose) {
    print("MCMC Complete!")
    print(paste("Running time: ", (proc.time()-time)[3]))
  }
  return(list("X1" = X1[,(burnin+1):iter], "X2" = X2[,(burnin+1):iter], 
              "lambda" = lambs[,(burnin+1):iter], "r" = r, "idx1" = idx1, "idx2" = idx2))
}