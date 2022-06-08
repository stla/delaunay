library(microbenchmark)
library(delaunay)
library(geometry)
library(uniformly)

pts <- runif_in_sphere(8000L, d = 3)

microbenchmark(
  delaunay = delaunay(pts, quick3d = TRUE),
  geometry = delaunayn(pts),
  times = 5
)