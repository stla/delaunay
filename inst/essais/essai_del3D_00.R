library(uniformly)
pts <- runif_in_sphere(12, d = 3)
cdt:::del3D_cpp(t(pts))
