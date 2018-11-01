############### compute Sigmas by FFT #################################
# covariance function for a single latitude band
# N should be fixed at 96 all through the calculation
# l is longitudinal lag, ranges from 0 to 2*pi*(N-1)/N, sysmmetric about pi
KL_l <- function(phi, alpha, nu, l, N = 96) {
  # spectral density
  fL_c <- sapply(0:(N-1), function(c) {
    phi/((alpha^2 + 4*(sin(c/N*pi))^2)^(nu+1/2))
  })
  return( mean(fL_c*cos(l*(0:(N-1)))) )
}




## DFT matrix   (eigenvector of a circulant matrix)
DFT_mat <- function(N){
  M <- outer(0:(N-1), 0:(N-1), function(i, j) exp(-(2*pi*1i*i*j)/N) ) / sqrt(N)
  return(M)
}




# function for getting Sigmas in temporal structure, assuming axially symmetry
# covariance matrix for one band at one time
Sigmas_fft <- function(phi, alpha, nu, N = 96) {
  col1 <- sapply(0:(N-1), function(x) {KL_l(phi = phi, alpha = alpha, nu = nu, l =  (2*pi*x)/N )})
  eig_val <- Re(fft(col1))
  eig_vec <- DFT_mat(N)
  Sigmas <- Re( eig_vec%*%diag(eig_val, ncol = N)%*%Conj(t(eig_vec)) )
  return(list(Sigmas = Sigmas, sig_val = eig_val, eig_vec = eig_vec))
}

# this function return a list that contains: 
# 1. Sigmas, 2. eigenvalue, 3. complex eigenvectors.









