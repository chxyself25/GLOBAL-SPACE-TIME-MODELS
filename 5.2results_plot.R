library(readxl)
A <- read_excel(path = "/Users/apple/Desktop/ISU 2018 fall/STAT 606/final_project/data/single_band_results.xlsx")
A <- as.data.frame(A)[1:42,]


fit_phi <- smooth.spline(x = A$lat, y = log(A$phi))
fit_alpha <- smooth.spline(x = A$lat, y = A$alpha)
fit_nu <- smooth.spline(x = A$lat, y = A$nu)

par(mfrow = c(2,2))
plot(A$lat, log(A$phi), pch = 16, 
     xlab = "latitude", ylab = "log(phi)", col = "blue")
lines(A$lat, fit_phi$y, col = "red")
plot(A$lat, A$alpha, pch = 16, 
     xlab = "latitude", ylab = "alpha", col = "blue")
lines(A$lat, fit_alpha$y, col = "red")
plot(A$lat, A$nu, pch = 16, 
     xlab = "latitude", ylab = "nu", col = "blue")
lines(A$lat, fit_nu$y, col = "red")
par(mfrow = c(1,1))