# data preprocessing, aggregation, organization, etc.
# Wed Oct 24 19:48:22 2018 ------------------------------
library(ncdf4)
library(dplyr)
library(doParallel)
registerDoParallel(cores=8)

# dir: directory containing all the nc files
# year is the vector of years(or range) you want to extract
ccsm3_data <- function(dir, years = NULL) {
  files <- list.files(dir)
  yrsin <- sapply(strsplit(files, split = "\\."), function(x){substr(x[4], 1,4)}) %>% unique
  years <- years[years %in% yrsin]
  mdoys <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  mdoys0 <- c(31,29,31,30,31,30,31,31,30,31,30,31)
  acomb3 <- function(...) abind::abind(..., along=3)
  res <- foreach(yr = years, .combine = 'acomb3', .multicombine = TRUE) %dopar% {
    fs <- sort(files[grepl(as.character(yr), files)])
    if (length(fs)==0) {return(NULL)}
    temps <- array(NA, dim = c(42, 96, 12))
    for (i in 1:12) {
      nc1 <- nc_open(paste0(dir, fs[i]))
      # remove the northernmost and southernmost latitude bands
      temps[,,i] <- t(ncvar_get(nc1, "TS")[,4:45]) # surface temprature
    }
    apply(temps, 1:2, function(x, y) {sum(x*y)/sum(y)}, y = ifelse(rep(yr%%4==0, 12), mdoys0, mdoys))
  }
  dimnames(res) <- list(paste("lat", 1:42, sep = ""), 
                        paste("long", 1:96, sep = ""), years)
  return(res)
}

# run 01, 2010-2020
dir <- "~/RDCEPData/700_from_historical_gradual_200years_640years/R05/"
years <- seq(2000, 2100, by = 5)
res <- ccsm3_data(dir, years)
saveRDS(res, file = "~/RDCEPData/700_from_historical_gradual_200years_640years/2000_5_2100_R05.rds")


# land fraction information, assuming consistent through years
nc1 <- nc_open("~/RDCEPData/700_from_historical_gradual_200years_640years/R01/2010_2050_med_R01.cam2.h0.2010-01.nc")
landfrac <- ncvar_get(nc1, "LANDFRAC")
landfrac[landfrac <= 0.5] <- 0
landfrac[landfrac > 0.5] <- 1
saveRDS(t(landfrac[,4:45]), file = "~/RDCEPData/700_from_historical_gradual_200years_640years/land_fraction_2010-01.rds")

long <- ncvar_get(nc1, "lon") # equally spaced
lat <- ncvar_get(nc1, "lat") # nearly equally spaced: 0.680158 to 3.711143
saveRDS(list(lat = lat[4:45], long = long), file = "T31grids_latlong.rds")
