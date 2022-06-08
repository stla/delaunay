library(uniformly)
pts <- cbind(runif_in_sphere(12, d = 2), rgamma(12, 10, 10))
cdt:::delXY_cpp(t(pts))
