library(uniformly)
pts <- runif_in_sphere(12, d = 2)
cdt:::del2D_cpp(t(pts))
