# covariance function for section 5.2: single latitude
# Sat Oct 27 11:23:01 2018 ------------------------------
library(Matrix)
library(foreach)
library(doParallel)
registerDoParallel(cores=16)


## F operator
DFT_mat <- function(N){
  M <- outer(0:(N-1), 0:(N-1), function(i, j) exp(-(2*pi*1i*i*j)/N) ) / sqrt(N)
  return(M)
}

# Foutier transform to Z
FFT <- function(D_mat, T.len = 15, N=96) {
  #FF <- outer(0:(N-1), 0:(N-1), function(i, j) exp(-(2*pi*1i)/N*i*j))
  # res <- apply(D_mat, 2, function(x) {
  #   unlist(lapply(1:T.len, function(y) FF%*%x[(N*(y-1)+1):(y*N)]))
  # })
  res <- apply(D_mat, 2, function(x) {
    unlist(lapply(1:T.len, function(y) fft(x[(N*(y-1)+1):(y*N)])/sqrt(N)))
  })
  return((res))
}

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


# function for getting Sigmas in temporal structure, assuming axially symmetry
# covariance matrix for one band at one time
Sigmas_fft <- function(phi, alpha, nu, N = 96) {
  col1 <- sapply(0:(N-1), function(x) {KL_l(phi = phi, alpha = alpha, nu = nu, l =  (2*pi*x)/N )})
  eig_val <- Re(fft(col1))
  #eig_vec <- DFT_mat(N)
  #Sigmas <- Re( eig_vec%*%diag(eig_val, ncol = N)%*%Conj(t(eig_vec)) )
  return(diag(eig_val))
  #return(list(Sigmas = Sigmas, sig_val = eig_val, eig_vec = eig_vec))
}

# this function return a list that contains: 
# 1. Sigmas, 2. eigenvalue, 3. complex eigenvectors.


# function for calculating Sigma for a single latitude band
# landfrac: a vector of length N including all land fractions in the latitude band
# phi, alpha, and nu are parameters in spectral density
# psis = c(psi0, psi1) (water and land coefficients)
Sigma_fft <- function(phi, alpha, nu, psis, landfrac, T.len = 15, N = 96) {
  #sigmas_out <- Sigmas_fft(phi = phi, alpha = alpha, nu = nu)
  #sigmas <- sigmas_out$Sigmas
  sigmas <- Sigmas_fft(phi = phi, alpha = alpha, nu = nu)
  psi <- rep(0,N)
  psi[which(landfrac==0)] <- psis[1]
  psi[which(landfrac==1)] <- psis[2]
  
  temps <- array(NA, dim = c(N, N, T.len))
  temps[,,1] <- sigmas
  for (i in 2:T.len) {
    MX <- sweep(sigmas,1,psi^(i-1),"*")
    MXM <- sweep(MX,2,psi^(i-1),"*")
    temps[,,i] <- MXM
  }
  mat <- foreach(i = 1:T.len, .combine = 'rbind') %dopar% {
    res <- matrix(0, ncol = T.len*N, nrow = N)
    if (i != T.len) {
      for (j in (i+1):T.len) { # j > i
        res[, (N*(j-1)+1):(N*j)] <- sweep( (apply(temps[,,1:i], 1:2, sum)), 2, psi^(j-i), "*" )
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


################ reml negative likelihood function for single band #######################
library(Rcpp)
library(RcppArmadillo)


code_inv_sympd <- 'arma::mat inv_sympd(arma::mat & X){
arma::mat Y;
Y = arma::inv_sympd(X);
return( Y );
}'

cppFunction(code=code_inv_sympd, depends="RcppArmadillo")


reml_neglik_single <- function(pars,landfrac1=landfrac[20,], 
                               T.len = 15, N = 96, R = 5, D = D_mat){
  # Foutier transform
  Z_mat <- FFT(D, T.len = T.len, N = N)
  phi1=pars[1]
  alpha1=pars[2]
  nu1=pars[3]
  psis1=pars[4:5]
  cat("paras = ", pars, ", \t")
  #cl <- makeCluster(8)
  #registerDoParallel(cl)
  Sigma_mat<-Sigma_fft(phi=phi1, alpha=alpha1, nu=nu1, psis=psis1, 
                       landfrac=landfrac1, T.len = T.len, N = N)
  # rearrange columns and rows
  arg.idx <- c(sapply(1:N, function(x) {seq(x, T.len*N, by = N)}))
  sigma_mat <- Sigma_mat[arg.idx, arg.idx]
  z_mat <- Z_mat[arg.idx,]
  sigma_mat_inv <- as.matrix(bdiag(foreach(n = 1:N) %dopar% {
    solve(sigma_mat[(T.len*(n-1)+1):(T.len*n), (T.len*(n-1)+1):(T.len*n)])
  }))
  #stopCluster(cl)
  eigen_sigma <- foreach(n=1:N, .combine = 'c') %dopar% {
    eigen(sigma_mat[(T.len*(n-1)+1):(T.len*n), (T.len*(n-1)+1):(T.len*n)], symmetric = TRUE)$values
  }
  A <- -0.5 * (R-1) * sum( log( eigen_sigma ) )
  B <- 0
  #inv_Sigma <- inv_sympd(Sigma_mat)
  for (i in 1:R){
    #B0 <- sum( z_mat[,i] * c( sigma_mat_inv%*%z_mat[,i] ) )
    B0 <- Conj(z_mat[,i]) %*% sigma_mat_inv %*% z_mat[,i] 
    B <- B + B0
  }
  loglike<- Mod(A - 0.5 * B) #- 0.5*T.len*N*(R-1)*log(2*pi) - 0.5*T.len*N*log(R)
  cat("loglike = ", loglike, "\n\n")
  return(-loglike)
}


