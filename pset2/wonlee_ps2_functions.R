M = 1e5

simYgivenTheta = function(theta, w, N) {
  rpois(N, theta*w)
}