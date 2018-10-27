# covariance function for section 5.2: single latitude
# Sat Oct 27 11:23:01 2018 ------------------------------
library(dplyr)
library(ggplot2)

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
Sigma <- function(phi, alpha, nu, landfrac, T.len = 15, N = 96) {
  sigmas <- Sigmas(phi = phi, alpha = alpha, nu = nu)
  Psi <- diag(landfrac)
  temp.mat <- matrix(0, ncol = T.len, nrow = T.len)
  temp.mat[lower.tri(temp.mat)] <- unlist(lapply(1:(T.len-1), function(x) {rep(x, T.len-x)}))
  temp.mat <- temp.mat + t(temp.mat)
  diag(temp.mat) <- 1:T.len
  mat <- kronecker(temp.mat, Psi%*%sigmas%*%Psi)
  # the first block is Sigmas
  mat[1:N, 1:N] <- sigmas
  return(mat)
}


