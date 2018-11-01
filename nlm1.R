source(file = "/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/covfun5.2.R")
landfrac <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/land_fraction_2010-01.rds")
R01 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/1870_20_2150_R01.rds")
R02 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/1870_20_2150_R02.rds")
R03 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/1870_20_2150_R03.rds")
R04 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/1870_20_2150_R04.rds")
R05 <- readRDS("/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/1870_20_2150_R05.rds")
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
  return((-1)*loglike)
}

nlmeControl(msTol = 1e-04) # 1e-07 by default
result<-nlm(whole_func, p=c(0.01, 0.1, 0.4, 0.05, 0.05), print.level = 2, gradtol = 1e-04)
## optim does not work better than nlm


####################################################################
source(file = "./covfun5.2.R")
landfrac <- readRDS("./Data/land_fraction_2010-01.rds")
R01 <- readRDS("./Data/1870_20_2150_R01.rds")
R02 <- readRDS("./Data/1870_20_2150_R02.rds")
R03 <- readRDS("./Data/1870_20_2150_R03.rds")
R04 <- readRDS("./Data/1870_20_2150_R04.rds")
R05 <- readRDS("./Data/1870_20_2150_R05.rds")
# estimate multiple bands
lat.bds <- seq(2, 42, by = 5)
initials <- list(b1 = c(0.01, 0.1, 1.8, 0.1, 0.1), 
                 b2 = c(0.01, 0.5, 3.7, 0.2, 0.2),
                 b3 = c(exp(-3), 0.7, 3.2, 0.1, 0.1),
                 b4 = c(exp(-2), 0.9, 2.5, 0.2, 0.1),
                 b5 = c(0.02, 0.1, 0.5, 1e-04, 1e-04),
                 b6 = c(0.02, 0.6, 1.7, 0.2, 0.1),
                 b7 = c(exp(-2), 0.6, 2.6, 0.1, 0.1),
                 b8 = c(exp(-3), 0.3, 2.2, 0.1, 0.1),
                 b9 = c(exp(-4), 0.2, 1.6, 0.3, 0.3))
pars.bds <- foreach(i = 1:length(lat.bds), .combine = 'rbind') %dopar% {
  R1=c(R01[lat.bds[i],,])
  R2=c(R02[lat.bds[i],,])
  R3=c(R03[lat.bds[i],,])
  R4=c(R04[lat.bds[i],,])
  R5=c(R05[lat.bds[i],,])
  R <- cbind(R1,R2,R3,R4,R5)
  R_bar <- apply(R,1,mean)
  
  D1=R1-R_bar
  D2=R2-R_bar
  D3=R3-R_bar
  D4=R4-R_bar
  D5=R5-R_bar
  
  D_mat = cbind(D1,D2,D3,D4,D5)
  nlmeControl(msTol = 1e-04)
  #result<-nlm(whole_func, p=initials[[paste0("b",i)]], D = D_mat, landfrac = landfrac[lat.bds[i],], 
              #print.level = 2, gradtol = 1e-04, iterlim = 50)
  ui.mat <- matrix(0, ncol = 5, nrow = 3) # constrain optim, first three parameters are positive
  ui.mat[1:3, 1:3] <- diag(1, 3)
  result <- constrOptim(theta = initials[[paste0("b",i)]], whole_func, D = D_mat, landfrac = landfrac[lat.bds[i],],
                        grad = NULL, ui = ui.mat, ci = rep(0, 3), outer.eps = 1e-04)
  result$estimate
}
