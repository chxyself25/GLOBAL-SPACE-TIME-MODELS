############################# input data ##########################
source(file = "./fft_funs_v2.R")
landfraction <- readRDS("~/RDCEPData/land_fraction_2010-01.rds")
lats <- readRDS("~/RDCEPData/T31grids_latlong.rds")$lat
R01 <- readRDS("~/RDCEPData/2000_5_2100_R01.rds")
R02 <- readRDS("~/RDCEPData/2000_5_2100_R02.rds")
R03 <- readRDS("~/RDCEPData/2000_5_2100_R03.rds")
R04 <- readRDS("~/RDCEPData/2000_5_2100_R04.rds")
R05 <- readRDS("~/RDCEPData/2000_5_2100_R05.rds")


library(readxl)
#fix.pars <- as.data.frame(na.omit(read_excel("./single_band_results.xlsx")))
fix.pars <- read.csv("./single_band_results_new.csv")

R <- list(R01, R02, R03, R04, R05)

ui.mat <- matrix(0, ncol = 4, nrow = 5)
ui.mat[1:3, 1:3] <- diag(1,3)
ui.mat[4,3] <- -1
ui.mat[5,3:4] <- c(-1,1)
result <- constrOptim(theta = c(0.1,0.1, 0.9,1.4), reml_neglik_multiple, 
                      T.len = 21, fix.pars = fix.pars, R = R, landfrac = landfraction,
                      method = "Nelder-Mead", ui = ui.mat, ci = c(0,0,0, -1,0), outer.eps = 1e-04, hessian = FALSE)

## maybe later for standard deviation: hessian matrix
par.est <- result$par
library(numDeriv)
xc <- hessian(reml_neglik_multiple, x = par.est, method="Richardson",
        T.len = 21, fix.pars = fix.pars, R = R, landfrac = landfrac)






##########################################################
###############single band################################
##########################################################
single.res <- NULL
for(i in 1:42) {
  R1=c(R01[i,,])
  R2=c(R02[i,,])
  R3=c(R03[i,,])
  R4=c(R04[i,,])
  R5=c(R05[i,,])
  R <- cbind(R1,R2,R3,R4,R5)
  R_bar <- apply(R,1,mean)
  
  D1=R1-R_bar
  D2=R2-R_bar
  D3=R3-R_bar
  D4=R4-R_bar
  D5=R5-R_bar
  
  D_mat = cbind(D1,D2,D3,D4,D5)
  
  #nlmeControl(msTol = 1e-04)
  #result<-nlm(whole_func, p=initials[[paste0("b",i)]], D = D_mat, landfrac = landfrac[lat.bds[i],], 
  #print.level = 2, gradtol = 1e-04, iterlim = 50)
  ui.mat <- matrix(0, ncol = 5, nrow = 3)
  ui.mat[1:3, 1:3] <- diag(1, 3)
  result <- constrOptim(theta = c(0.2, 0.5, 1, 0.1, 0.1), reml_neglik_single, D = D_mat, T.len = 21, landfrac = landfrac[i,],
                        grad = NULL, ui = ui.mat, ci = rep(0,3), outer.eps = 1e-04)
  #control = list(reltol = 1e-04, ndeps = 1e-04))
  single.res <- rbind(single.res, data.frame(id = i, lat = lats[i], phi = result$par[1], alpha = result$par[2], nu = result$par[3],
                                       psi0 = result$par[4], psi1 = result$par[5], loglike = -result$value))
}

