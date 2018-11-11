# covariance function for section 5.2: single latitude
library(Matrix)

## F operator
DFT_mat <- function(N){
  M <- outer(0:(N-1), 0:(N-1), function(i, j) exp(-(2*pi*1i*i*j)/N) ) / sqrt(N)
  return(M)
}

# Foutier transform to Z
FFT <- function(D_mat, T.len = 15, N=96) {
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
  return(eig_val)
}


# function for calculating Sigma for a single latitude band
# landfrac: a vector of length N including all land fractions in the latitude band
# phi, alpha, and nu are parameters in spectral density
# psis = c(psi0, psi1) (water and land coefficients)
Sigma_fft <- function(phi, alpha, nu, psis, T.len = 15, N = 96) {
  sigmas <- Sigmas_fft(phi = phi, alpha = alpha, nu = nu)
  psi <- rep(1,N)*psis
  
  temps <- matrix(NA, ncol = N, nrow = T.len)
  temps[1,] <- sigmas
  for (i in 2:T.len) {
    temps[i,] <- sigmas*psi^(2*i-2)  
  }
  temps1 <- apply(temps,2,cumsum)
  mat <- array(NA, dim = c(T.len, T.len, N))
  for (i in 1:(T.len-1)){
    for (j in (i+1):T.len){ # j>i
      mat[i,j,] <- temps1[i,]*psi^(j-i)
      mat[j,i,] <- mat[i,j,]
    }
  }
  for (i in 1:T.len){
    mat[i,i,] <- temps1[i,]
  }
  return(mat)
}


## rearrange Sigma and Z_mat
rearr_Sigma <- function(Sigma){
  N <- dim(Sigma)[3]
  diag_block <- vector("list",N)
  for (i in 1:N){
    diag_block[[i]] <- Sigma[,,i]
  }
  return(diag_block)
}


rearr_Data <- function(Z_mat, T.len = 15, N = 96){
  Z_block <- vector("list",N)
  for (i in 1:N){
    Z_block[[i]] <- Z_mat[seq(i, T.len*N, by = N),]
  }
  return(Z_block)
}




################ reml negative likelihood function for single band #######################



reml_neglik_reduced <- function(pars, T.len = 15, N = 96, R = 5, D = D_mat){
  # Foutier transform
  Z_mat <- FFT(D, T.len = T.len, N = N)
  phi1=pars[1]
  alpha1=pars[2]
  nu1=pars[3]
  psis1=pars[4]
  cat("paras = ", pars, ", \t")
  Sigma_mat<-Sigma_fft(phi=phi1, alpha=alpha1, nu=nu1, psis=psis1, T.len = T.len, N = N)
  # rearrange columns and rows

  Sigma_block <- rearr_Sigma(Sigma_mat)
  Z_block <- rearr_Data(Z_mat, T.len, N)
  Sigma_block_inv <- lapply(Sigma_block, solve)
  eigen_sigma_block <- lapply(Sigma_block, function(x){
    return(eigen(x,symmetric = TRUE)$values)
    })
  eigen_sigma <- unlist(eigen_sigma_block)
  A <- -0.5 * (R-1) * sum( log( eigen_sigma ) )
  B <- 0
  for (i in 1:R){
    for (j in 1:N){
      B0 <- Re( Conj(t( Z_block[[j]][,i] )) %*% Sigma_block_inv[[j]] %*% Z_block[[j]][,i] )
      B <- B + B0
    }
  }
  loglike<- A - 0.5 * B #- 0.5*T.len*N*(R-1)*log(2*pi) - 0.5*T.len*N*log(R)
  cat("loglike = ", loglike, "\n\n")
  return(-loglike)
}














##### run
R01 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/2000_5_2100/2000_5_2100_R01.rds")
R02 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/2000_5_2100/2000_5_2100_R02.rds")
R03 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/2000_5_2100/2000_5_2100_R03.rds")
R04 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/2000_5_2100/2000_5_2100_R04.rds")
R05 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/2000_5_2100/2000_5_2100_R05.rds")

lat_id <- 32
R01=c(R01[lat_id,,])
R02=c(R02[lat_id,,])
R03=c(R03[lat_id,,])
R04=c(R04[lat_id,,])
R05=c(R05[lat_id,,])
R <- cbind(R01,R02,R03,R04,R05)
R_bar <- apply(R,1,mean)

D1=R01-R_bar
D2=R02-R_bar
D3=R03-R_bar
D4=R04-R_bar
D5=R05-R_bar
# T=21;N=96;M=1;R=5
D_mat = cbind(D1,D2,D3,D4,D5)


################# constrOptim
UI <- diag(4)
UI
CI <- c(0,0,0,0)
CI

result <- constrOptim(theta = c(0.5, 0.001, 0.004, 0.005), f = reml_neglik_reduced, grad = NULL,
                      ui = UI, ci = CI)

