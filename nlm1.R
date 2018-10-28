source(file = "/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/covfun5.2.R")
landfrac <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/land_fraction_2010-01.rds")
R01 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/1870_20_2150_R01.rds")
R02 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/1870_20_2150_R02.rds")
R03 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/1870_20_2150_R03.rds")
R04 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/1870_20_2150_R04.rds")
R05 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/1870_20_2150_R05.rds")
# R_bar=mean(c(R01,R02,R03,R04,R05))
# dim(R05)
# D1=R01-R_bar;D2=R02-R_bar;D3=R03-R_bar;D4=R04-R_bar;D5=R05-R_bar
# D1=matrix(D1,1,prod(c(dim(D1)[1:3])))
# D2=matrix(D1,1,prod(c(dim(D2)[1:3])))
# D3=matrix(D1,1,prod(c(dim(D3)[1:3])))
# D4=matrix(D1,1,prod(c(dim(D4)[1:3])))
# D5=matrix(D1,1,prod(c(dim(D5)[1:3])))
R01=c(R01[20,,])
R02=c(R02[20,,])
R03=c(R03[20,,])
R04=c(R04[20,,])
R05=c(R05[20,,])
R <- cbind(R01,R02,R03,R04,R05)
R_bar <- apply(R,1,mean)

D1=R01-R_bar
D2=R02-R_bar
D3=R03-R_bar
D4=R04-R_bar
D5=R05-R_bar


# T=15;N=96;M=1;R=5
D_mat = cbind(D1,D2,D3,D4,D5)
# loglike<--0.5*T*N*M*(R-1)*log(2*pi)
#         -0.5*(R-1)*log(det(Sigma(phi=ph_2, alpha=alpha_2, nu=nu_2, psis=psis_2, landfrac=landfrac_2, T.len = 15, N = 96)))
#         -0.5*T*N*M*log(R)
#         -0.5*t(D)%*%solve(kronecker(diag(R),Sigma(phi=ph_2, alpha=alpha_2, nu=nu_2, psis=psis_2, landfrac=landfrac_2, T.len = 15, N = 96)))%*%D
whole_func<-function(pars,landfrac_2=landfrac[20,], T.len = 15, N = 96, D = D_mat){
  phi_2=pars[1]
  alpha_2=pars[2]
  nu_2=pars[3]
  psis_2=pars[4:5]
  M=1
  R=5
  Sigma_mat<-Sigma(phi=phi_2, alpha=alpha_2, nu=nu_2, psis=psis_2, landfrac=landfrac_2, T.len = T.len, N = N)
  eigen_sigma <- svd(Sigma_mat)$d
  A <- -0.5 * (R-1) * sum( log( eigen_sigma ) )
  B <- 0
  for (i in 1:R){
    B0 <- sum( D[,i] * c(solve(Sigma_mat)%*%D[,i]) )
    B <- B + B0
  }
  loglike<- A - 0.5 * B
  return(loglike)
}


# try<-Sigma(phi=0.01, alpha=0.1, nu=0.4, psis=c(0.05,0.05), landfrac=landfrac[20,], T.len = 15, N = 96)
# try2<-solve(kronecker(diag(5),try))

result<-nlm(whole_func, p=c(0.01, 0.1, 0.4, 0.05, 0.05) )




### test
Sigma_mat_test <- Sigma(phi=0.01, alpha=0.01, nu=0.01, psis=c(0.05,0.05), 
                        landfrac=landfrac[20,], T.len = 15, N = 96)

svd_sig <- svd(Sigma_mat_test)
svd_sig$d

sum( log(svd_sig$d))
