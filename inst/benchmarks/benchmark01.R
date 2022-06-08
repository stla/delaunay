library(cdt)
library(RCDT)
library(microbenchmark)

nsides <- 1440L
angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
points <- cbind(cos(angles), sin(angles))
points <- rbind(points, points/1.5)
# constraint edges
indices <- 1L:nsides
edges_outer <- cbind(
  indices, c(indices[-1L], indices[1L])
)
edges_inner <- edges_outer + nsides
edges <- rbind(edges_outer, edges_inner)

microbenchmark(
  RCDT = delaunay(points, edges),
  cdt = cdel(points, edges),
  times = 5
)