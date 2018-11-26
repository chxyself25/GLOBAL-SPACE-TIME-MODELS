# covariance function for section 5.2 and 5.3: single latitude and multiple latitude
library(Matrix)
library(doParallel)
registerDoParallel(cores = 8)
library(data.table)

## F operator
DFT_mat <- function(N){
  M <- outer(0:(N-1), 0:(N-1), function(i, j) exp(-(2*pi*1i*i*j)/N) ) / sqrt(N)
  return(M)
}

# Foutier transform to Z: transform for each latitude and year
# ordered by longitude, year, latitude, but will be rearranged later
# input is original data, all 5 realizations
FFT_multiple <- function(R1, R2, R3, R4, R5, psi.same = FALSE) {
  R_bar <- (R1+R2+R3+R4+R5)/5
  D <- list(R1-R_bar, R2-R_bar, R3-R_bar, R4-R_bar, R5-R_bar)
  T.len <- dim(R_bar)[3]; N <- dim(R_bar)[2]; M <- dim(R_bar)[1]
  res <- rep(list(NULL), 5)
  for (r in 1:5) {
    Dr <- D[[r]]
    ## fft: ordered by year, longitude, and latitude
    Zr <- c(sapply(1:T.len, function(t) {c(t(apply(Dr[,,t], 1, function(x) fft(x)/sqrt(N))))}))
    if (psi.same) {
      ## rearrange to order of longitude, year and latitude
      res[[r]] <- lapply(1:N, function(n) {Zr[c(outer(1:M, 1:T.len, function(m, t) {M*N*(t-1)+(n-1)*M+m}))]})
    } else {
      res[[r]] <- Zr
    }
  }
  return(res)
}

## Foutier transform to Z for single band, input is organized D_mat
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

## inverse of Sigma, which can be expressed in terms of Sigmas
Sigmas_fft_inv <- function(phi, alpha, nu, psis, landfrac) {
  N <- 96
  sigmas <- Sigmas_fft(phi, alpha, nu)
  psi <- rep(0, N)
  psi[which(landfrac==0)] <- psis[1]
  psi[which(landfrac==1)] <- psis[2]
  # inverse of Sigmas
  sigmas_inv <- 1/sigmas
  # log determinant
  log_det <- sum(log(sigmas))
  # intermediate results in Sigma inverse
  A.mat <- sigmas_inv + psi*sigmas_inv*psi
  B.mat <- -psi*sigmas_inv
  
  return(list('A' = A.mat, 'B' = B.mat, 'inv' = sigmas_inv, 'log_det' = log_det))
}

## Sigmas for multiple bands, each block is diagonal, should be diagonal block matrix after rearrange
## phis are all phi estimation from 42 latitudes, same with alphas and nus
MSigmas_fft <- function(phis, alphas, nus, xi, tau) {
  M <- 42; N <- 96
  # spectral densities for 42 latitude bands
  fL_c <- foreach(m = 1:M, .combine = "rbind") %do% {
    sapply(0:(N-1), function(c) {
      phis[m]/((alphas[m]^2 + 4*(sin(c/N*pi))^2)^(nus[m]+1/2))
    })
  }
  # cross spectral densities
  lagm <- 3.71*pi/180
  pho_c <- foreach(d = 1:(M-1), .combine = "rbind") %do% {
    sapply(0:(N-1), function(c) {
      (xi/(1 + 4*(sin(c/N*pi))^2)^tau)^(d*lagm)
    })
  }
  # cosine functions
  cos.mat <- outer(0:(N-1), 0:(N-1), function(n, c) {
    cos((2*pi*n/N) * c)
  })
  res <- array(NA, dim = c(M, M, N))
  for (i in 1:(M-1)) {
    for (j in (i+1):M) { # j > i
      fLL_c <- sqrt(fL_c[i,]*fL_c[j,])*pho_c[j-i,]
      res[i,j,] <- Re(fft(cos.mat %*% fLL_c / N))
      res[j,i,] <- res[i,j,]
    }
  }
  for (i in 1:M) {
    res[i,i,] <- Re(fft(cos.mat %*% fL_c[i,] / N))
  }
  # rearrange the matrix such that it is diagonal
  rearr_Sigma(res)
}

## inverse of Sigma, which can be expressed by Sigmas, psis is a vector of length 2
MSigmas_fft_inv <- function(phi, alpha, nu, psis, landfrac, xi, tau) {
  sigmas <- MSigmas_fft(phis = phi, alphas = alpha, nus = nu, xi = xi, tau = tau)
  M = 42; N = 96
  landfrac[which(landfrac==0)] <- psis[1]
  landfrac[which(landfrac==1)] <- psis[2]
  psi <- lapply(1:N, function(n) {landfrac[,n]})
  # log of determinant of sigmas
  eigen_values <- lapply(sigmas, function(x) {eigen(x)$values})
  log_det <- sum(log(unlist(eigen_values)))
  # inverse of sigmas
  sigmas_invlist <- lapply(sigmas, solve)
  # intermediate results in Sigma inverse
  A.mat <- lapply(1:N, function(n) {
    sigmas_invlist[[n]] + sweep(sweep(sigmas_invlist[[n]], 1, psi[[n]], "*"), 2, psi[[n]], "*")
    #sigmas_invlist[[n]] + diag(psi[[n]]) %*% sigmas_invlist[[n]] %*% diag(psi[[n]])
  })
  B.mat <- lapply(1:N, function(n) {
    -1 * sweep(sigmas_invlist[[n]], 1, psi[[n]], "*")
  })
  # construct inverse of Sigma
  #band1 <- rbind(A.mat, t(B.mat))
  #band2 <- rbind(B.mat, A.mat, t(B.mat))
  #band3 <- rbind(B.mat, sigmas_inv)
  #return(list('band1' = band1, 'band2' = band2, 'band3' = band3, 'A' = A.mat, 'B' = B.mat, 'inv' = sigmas_inv))
  return(list('invlist' = sigmas_invlist, 'log_det' = log_det, 
              'A' = A.mat, 'B' = B.mat))
}

## multiplication of D and inverse of Sigma which can be expressed in terms of Sigmas
Dxmat_list <- function(mat.list, Dryear) {
  N <- 96; M <- 42
  res <- lapply(1:N, function(n) {
    t(Dryear[((n-1)*M+1):(n*M)]) %*% mat.list[[n]]
  })
  return(unlist(res))
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

reml_neglik_single <- function(pars, T.len = 15, D, landfrac){
  R <- 5; N <- 96
  # Foutier transform
  Z_mat <- FFT(D, T.len = T.len, N = N)
  phi1=pars[1]
  alpha1=pars[2]
  nu1=pars[3]
  psis1=pars[4:5]
  cat("paras = ", pars, ", \t")
  Sigmas_res<-Sigmas_fft_inv(phi=phi1, alpha=alpha1, nu=nu1, psis=psis1, landfrac)
  
  ## the first relative term in reml
  log_det <- T.len * Sigmas_res$log_det
  v1 <- -0.5 * (5-1) * log_det
  
  ## the second relative term in reml
  A.mat <- Sigmas_res$A
  B.mat <- Sigmas_res$B # B transpose is the same as B
  sigmas_inv <- Sigmas_res$inv
  
  # Sigma_inv <- matrix(0, ncol = 96*21, nrow = 96*21)
  # for (t in 2:20) {
  #   Sigma_inv[((t-2)*N+1):((t+1)*N),((t-1)*N+1):(t*N)] <- rbind(diag(B.mat), diag(A.mat), diag(B.mat))
  # }
  # Sigma_inv[1:(2*N),1:N] <- rbind(diag(A.mat), diag(B.mat))
  # Sigma_inv[((T.len-2)*N+1):(T.len*N), ((T.len-1)*N+1):(T.len*N)] <- rbind(diag(B.mat), diag(sigmas_inv))
  
  
  v2 <- 0
  for (r in 1:5) {
    Dr <- Z_mat[,r]
    res1 <- Dr[1:N] * A.mat + Dr[(N+1):(2*N)] * B.mat
    res2 <- lapply(2:(T.len-1), function(t) {
      Dr[((t-2)*N+1):((t-1)*N)] * B.mat + Dr[((t-1)*N+1):(t*N)] * A.mat + Dr[(t*N+1):((t+1)*N)] * B.mat
    })
    res3 <- Dr[((T.len-2)*N+1):((T.len-1)*N)] * B.mat + Dr[((T.len-1)*N+1):(T.len*N)] * sigmas_inv
    v2 <- v2 + sum(c(res1, unlist(res2), res3) * Conj(Dr))
  }
  loglike <- -0.5*Re(v2) + v1
  cat("loglike = ", loglike, "\n")
  return(-loglike)
}

## loglikelihood functions for multiple bands, psi0 and psi1 are not the same
reml_neglik_multiple <- function(pars, T.len = 21, fix.pars, R, landfrac) {
  R1 <- R[[1]]; R2 <- R[[2]]; R3 <- R[[3]]; R4 <- R[[4]]; R5 <- R[[5]]
  M <- 42; N <- 96
  # Foutier transform: ordered by year, longitude and latitude
  Z_block <- FFT_multiple(R1, R2, R3, R4, R5, psi.same = FALSE)
  psis=pars[1:2] # psi0 and psi1
  #psis <- c(0.0002632907, 0.005489152)
  xi <- pars[3]
  tau <- log(pars[4], 5) # input is 5^tau
  # psis <- c(0.0001,0.0001)
  # xi <- 0.9; tau <- 0.2
  cat("paras = ", pars, ", \t")
  
  # Sigmas calculation results
  Sigmas_res <- MSigmas_fft_inv(phi = fix.pars$phi, alpha = fix.pars$alpha, nu = fix.pars$nu, 
                                psis = psis, landfrac = landfrac, xi = xi, tau = tau)
  ## the first relative term in reml
  log_det <- T.len * Sigmas_res$log_det
  v1 <- -0.5 * (5-1) * log_det
  
  ## the second relative term in reml
  A.list <- Sigmas_res$A
  B.list <- Sigmas_res$B; tB.list <- lapply(B.list, t)
  sigmas_invlist <- Sigmas_res$invlist
  v2 <- 0
  for (r in 1:5) {
    Dr <- Z_block[[r]]
    res1 <- Dxmat_list(A.list, Dr[1:(M*N)]) + Dxmat_list(tB.list, Dr[(M*N+1):(2*M*N)])
    res2 <- lapply(2:(T.len-1), function(t) {
      Dxmat_list(B.list, Dr[((t-2)*M*N+1):((t-1)*M*N)]) + Dxmat_list(A.list, Dr[((t-1)*M*N+1):(t*M*N)]) + Dxmat_list(tB.list, Dr[(t*M*N+1):((t+1)*M*N)])
    })
    res3 <- Dxmat_list(B.list, Dr[((T.len-2)*M*N+1):((T.len-1)*M*N)]) + Dxmat_list(sigmas_invlist, Dr[((T.len-1)*M*N+1):(T.len*M*N)])
    v2 <- v2 + sum(c(res1, unlist(res2), res3) * Conj(Dr))
  }
  loglike <- -0.5*Re(v2) + v1
  cat("loglike = ", loglike, "\n")
  return(-loglike)
}




########################################################################
#####################useless functions##################################
########################################################################
## loglikelihood functions for single laititude bands: reduced version: psi0 and psi1 are the same
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
## loglikelihood functions for multiple laititude bands: reduced version: psi0 and psi1 are the same
reml_neglik_mreduced <- function(pars, T.len = 15, fix.pars, R){
  R1 <- R[[1]]; R2 <- R[[2]]; R3 <- R[[3]]; R4 <- R[[4]]; R5 <- R[[5]]
  
  # Foutier transform: ordered by longitude, year and latitude
  Z_block <- FFT_multiple(R1, R2, R3, R4, R5, psi.same = TRUE)
  psis=pars[1]
  xi <- pars[2]
  tau <- pars[3]
  cat("paras = ", pars, ", \t")
  
  Sigma_block<-Sigma_fft(phi=fix.pars$phi, alpha=fix.pars$alpha, nu=fix.pars$nu, psis=psis, xi=xi, tau=tau, T.len = T.len)
  
  # Sigma_block_inv <- lapply(Sigma_block, solve)
  # eigen_sigma_block <- lapply(Sigma_block, function(x){
  #   return(eigen(x,symmetric = TRUE)$values)
  # })
  # eigen_sigma <- unlist(eigen_sigma_block)
  eigen_blocks <- lapply(Sigma_block, function(x) {eigen(x, symmetric = TRUE)})
  Sigma_block_inv <- lapply(eigen_blocks, function(e) {(e$vectors %*% diag(1/e$values) %*% t(e$vectors))})
  eigen_sigma <- unlist(lapply(eigen_blocks, function(e) {e$values}))
  A <- -0.5 * (5-1) * sum( log( eigen_sigma ) )
  B <- 0
  for (r in 1:5){
    B0 = sapply(1:N, function(n) {
      Re( Conj(t(Z_block[[r]][[n]])) %*% Sigma_block_inv[[n]] %*% Z_block[[r]][[n]] )
    })
    B = B + sum(B0)
  }
  loglike<- A - 0.5 * B #- 0.5*T.len*N*(R-1)*log(2*pi) - 0.5*T.len*N*log(R)
  cat("loglike = ", loglike, "\n\n")
  return(-loglike)
}

# cross covariance function, may not be useful
KLL_l <- function(phis, alphas, nus, xi, tau, L1, L2, l) {
  N <- 96
  # spectral density
  f1 <- sapply(0:(N-1), function(c) {
    phis[1]/((alphas[1]^2 + 4*(sin(c/N*pi))^2)^(nus[1]+1/2))
  })
  f2 <- sapply(0:(N-1), function(c) {
    phis[2]/((alphas[2]^2 + 4*(sin(c/N*pi))^2)^(nus[2]+1/2))
  })
  # cross spectral density
  pho_c <- sapply(0:(N-1), function(c) {
    (xi/(1 + 4*(sin(c/N*pi))^2)^tau)^(abs(L1-L2))
  })
  sumres <- sqrt(f1*f2) * pho_c * cos(l*(0:(N-1)))
  mean(sumres)
}

# function for calculating Sigma for a single latitude band or multiple latitude bands
# landfrac: a vector of length N including all land fractions in the latitude band
# phi, alpha, and nu are parameters in spectral density
# psis = psi0 = psi1
Sigma_fft <- function(phi, alpha, nu, psis, xi=NULL, tau=NULL, T.len = 15) {
  if (any(is.null(c(xi, tau)))) {
    sigmas <- Sigmas_fft(phi, alpha, nu)
    psi <- rep(1,96)*psis
    temps <- matrix(NA, ncol = 96, nrow = T.len)
    temps[1,] <- sigmas
    for (i in 2:T.len) {
      temps[i,] <- sigmas*psi^(2*i-2)  
    }
    temps1 <- apply(temps,2,cumsum)
    mat <- array(NA, dim = c(T.len, T.len, 96))
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
  } else {
    sigmas <- MSigmas_fft(phis = phi, alphas = alpha, nus = nu, xi = xi, tau = tau)
    M = 42; N = 96
    psi <- rep(1,N)*psis
    temps <- array(NA, dim = c(M, M, N, T.len))
    temps[,,,1] <- sigmas
    for (i in 2:T.len) {
      temps[,,,i] <- sigmas*(psi^(2*i-2))
    }
    temps1 <- apply(temps, 1:3, cumsum) # dimension becomes T.len, M, M, N
    mat <- array(NA, dim = c(T.len, T.len, M, M, N))
    for (i in 1:(T.len-1)) {
      for (j in (i+1):T.len) {# j > i
        mat[i,j,,,] <- temps1[i,,,] * (psi^(j-i))
        mat[j,i,,,] <- mat[i,j,,,]
      } 
    }
    for (i in 1:T.len) {
      mat[i,i,,,] <- temps1[i,,,]
    }
    #dimnames(mat) <- list("t1", "t2", "lat1", "lat2", "long")
    ## mat: single latitude: T.len, T.len, N; multiple latitudes: T.len, T.len, M, M, N
    # rearrange matrix: NxN blocks: each is TMxTM matrix ordered by year and latitude
    lapply(1:N, function(x) {
      matx <- mat[,,,,x]
      do.call("rbind", lapply(1:T.len, function(t1) {
        do.call("cbind", lapply(1:T.len, function(t2) {matx[t1,t2,,]}))
      }))
    })
  }
}