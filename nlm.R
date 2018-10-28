# covariance function for section 5.2: single latitude
# Sat Oct 27 11:23:01 2018 ------------------------------
library(Matrix)
library(dplyr)
library(foreach)
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

landfrac <- readRDS("~/Documents/Iowa State/2018Fall/606/Project/Data/land_fraction_2010-01.rds")
R01 <- readRDS("~/Documents/Iowa State/2018Fall/606/Project/Data/1870_20_2150_R01.rds")
R02 <- readRDS("~/Documents/Iowa State/2018Fall/606/Project/Data/1870_20_2150_R02.rds")
R03 <- readRDS("~/Documents/Iowa State/2018Fall/606/Project/Data/1870_20_2150_R03.rds")
R04 <- readRDS("~/Documents/Iowa State/2018Fall/606/Project/Data/1870_20_2150_R04.rds")
R05 <- readRDS("~/Documents/Iowa State/2018Fall/606/Project/Data/1870_20_2150_R05.rds")
# R_bar=mean(c(R01,R02,R03,R04,R05))
# dim(R05)
# D1=R01-R_bar;D2=R02-R_bar;D3=R03-R_bar;D4=R04-R_bar;D5=R05-R_bar
# D1=matrix(D1,1,prod(c(dim(D1)[1:3])))
# D2=matrix(D1,1,prod(c(dim(D2)[1:3])))
# D3=matrix(D1,1,prod(c(dim(D3)[1:3])))
# D4=matrix(D1,1,prod(c(dim(D4)[1:3])))
# D5=matrix(D1,1,prod(c(dim(D5)[1:3])))
R01=R01[20,,];R02=R02[20,,];R03=R03[20,,];R04=R04[20,,];R05=R05[20,,]
R_bar=mean(c(R01,R02,R03,R04,R05))
D1=R01-R_bar;D2=R02-R_bar;D3=R03-R_bar;D4=R04-R_bar;D5=R05-R_bar
D1=matrix(D1,1,prod(c(dim(D1)[1:2])))
D2=matrix(D1,1,prod(c(dim(D2)[1:2])))
D3=matrix(D1,1,prod(c(dim(D3)[1:2])))
D4=matrix(D1,1,prod(c(dim(D4)[1:2])))
D5=matrix(D1,1,prod(c(dim(D5)[1:2])))

# T=15;N=96;M=1;R=5
D=matrix(c(D1,D2,D3,D4,D5),1,5*dim(D1)[2])
# loglike<--0.5*T*N*M*(R-1)*log(2*pi)
#         -0.5*(R-1)*log(det(Sigma(phi=ph_2, alpha=alpha_2, nu=nu_2, psis=psis_2, landfrac=landfrac_2, T.len = 15, N = 96)))
#         -0.5*T*N*M*log(R)
#         -0.5*t(D)%*%solve(kronecker(diag(R),Sigma(phi=ph_2, alpha=alpha_2, nu=nu_2, psis=psis_2, landfrac=landfrac_2, T.len = 15, N = 96)))%*%D
whole_func<-function(pars,landfrac_2=landfrac[20,],T.len = 15, N = 96,D=D){
  phi_2=pars[1]
  alpha_2=pars[2]
  nu_2=pars[3]
  psis_2=pars[4]
  T=15;N=96;M=1;R=5
  Sigma_mat<-Sigma(phi=phi_2, alpha=alpha_2, nu=nu_2, psis=psis_2, landfrac=landfrac_2, T.len = 15, N = 96)
  loglike<--0.5*T*N*M*(R-1)*log(2*pi)-0.5*(R-1)*log(det(Sigma_mat))-0.5*T*N*M*log(R)-0.5*D%*%solve(kronecker(diag(R),Sigma_mat))%*%t(D)
  return(loglike)
}


# try<-Sigma(phi=0.01, alpha=0.5, nu=2, psis=c(0.1,0.1), landfrac=landfrac[20,], T.len = 15, N = 96)
# try2<-solve(kronecker(diag(5),try))

result<-nlm(whole_func,p=c(0.01,0.5,2,c(0.1,0.1)))

################################################################
###################old version: slow############################
################################################################
Sigma <- function(phi_2, alpha_2, nu_2, psis_2, landfrac, T.len = 15, N = 96) {
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
