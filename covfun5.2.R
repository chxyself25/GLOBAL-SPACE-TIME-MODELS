# covariance function for section 5.2: single latitude
# Sat Oct 27 11:23:01 2018 ------------------------------
library(Matrix)
library(doParallel)
registerDoParallel(cores = 8)

# covariance function for a single latitude band
# N should be fixed at 96 all through the calculation
# l is longitudinal lag, ranges from 0 to 2*pi*(N-2)/N, sysmmetric about pi
KL_l <- function(phi, alpha, nu, l, N = 96) {
  # spectral density
  fL_c <- sapply(0:(N-1), function(c) {
    phi/((alpha^2 + 4*(sin(c/N*pi))^2)^(nu+1/2))
  })
  mean(fL_c*cos(l*(0:(N-1))))
}
# check the sysmetry initially 
l <- 2*pi*seq(0, 95)/96
KLs <- sapply(l, function(x) {KL_l(phi = 1, alpha = 0.1, nu = 0.1, l = x)})
plot(0:95, KLs, pch = 20)

# function for getting Sigmas in temporal structure, assuming axially symmetry
Sigmas <- function(phi, alpha, nu, N = 96) {
  entries <- combn(0:(N-1), m = 2)
  covs <- apply(entries, 2, function(x) {KL_l(phi = phi, alpha = alpha, nu = nu, l = abs(x[1]-x[2])*2*pi/N)})
  mat <- matrix(0, ncol = N, nrow = N)
  mat[lower.tri(mat)] <- covs
  mat <- mat + t(mat)
  diag(mat) <- KL_l(phi = phi, alpha = alpha, nu = nu, l = 0)
  return(mat)
}

# function for calculating Sigma for a single latitude band
# landfrac: a vector of length N including all land fractions in the latitude band
# phi, alpha, and nu are parameters in spectral density
# psis = c(psi0, psi1) (water and land coefficients)
Sigma <- function(phi, alpha, nu, psis, landfrac, T.len = 15, N = 96) {
  sigmas <- Sigmas(phi = phi, alpha = alpha, nu = nu)
  Psi <- diag(sapply(landfrac, function(x) {psis[x+1]}))
  # considering speed, calculate these only once, but may have memory 
  temps <- array(NA, dim = c(N, N, T.len))
  temps[,,1] <- sigmas
  for (i in 2:T.len) {
    temps[,,i] <- (Psi^(i-1)) %*% sigmas %*% (Psi^(i-1)) 
  }
  mat <- foreach(i = 1:T.len, .combine = 'rbind') %dopar% {
    res <- matrix(0, ncol = T.len*N, nrow = N)
    if (i != T.len) {
      for (j in (i+1):T.len) { # j > i
        res[, (N*(j-1)+1):(N*j)] <- apply(temps[,,(j-i+1):j], 1:2, sum)
      }
    }
    res
  }
  mat <- mat+t(mat)
  mat <- mat + as.matrix(bdiag(foreach(i = 1:T.len) %dopar% {
    apply(temps[,,1:i], 1:2, sum)
  }))
  return(mat)
}


################################################################
###################old version: slow############################
################################################################
Sigma <- function(phi, alpha, nu, psis, landfrac, T.len = 15, N = 96) {
  sigmas <- Sigmas(phi = phi, alpha = alpha, nu = nu)
  Psi <- diag(sapply(landfrac, function(x) {psis[x+1]}))
  # considering speed, calculate these only once, but may have memory issue
  temps <- array(NA, dim = c(N, N, T.len))
  temps[,,1] <- sigmas
  for (i in 2:T.len) {
    temps[,,i] <- (Psi^(i-1)) %*% sigmas %*% (Psi^(i-1)) 
  }
  mat <- matrix(0, ncol = T.len*N, nrow = T.len*N)
  for (i in 1:T.len) {
    for (j in i:T.len) {
      if(i == j) {
        mat[(N*(i-1)+1):(N*i), (N*(j-1)+1):(N*j)] <- apply(temps[,,1:i], 1:2, sum)
      } else { # j > i
        mat[(N*(i-1)+1):(N*i), (N*(j-1)+1):(N*j)] <- apply(temps[,,(j-i+1):j], 1:2, sum) 
        mat[(N*(j-1)+1):(N*j), (N*(i-1)+1):(N*i)] <- apply(temps[,,(j-i+1):j], 1:2, sum)
      }
    }
  }
  return(mat)
}


