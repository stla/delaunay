library(delaunay)
library(cgalMeshes)

nsides <- 6L
angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
outer_points <- cbind(cos(angles), sin(angles))
inner_points <- outer_points / 2
middle_points <- outer_points / 1.5
points <- rbind(outer_points, inner_points, middle_points)
# constraint edges
indices <- 1L:nsides
edges <- cbind(
  indices, c(indices[-1L], indices[1L])
)
edges <- rbind(edges, edges + nsides)
# constrained Delaunay triangulation
d <- delaunay(points, constraints = edges)

#####



#####
points <- uniformly::runif_in_annulus(100, c(0, 0), 2, 1)

nsides <- 6L
angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
opoints <- 2*cbind(cos(angles), sin(angles))
ipoints <- cbind(cos(angles), sin(angles))
points <- rbind(opoints, ipoints, points)
# constraint edges
indices <- 1L:nsides
edges <- cbind(
  indices, c(indices[-1L], indices[1L])
)
constraints <- rbind(edges, edges + nsides)


d <- delaunay(points, constraints = constraints)
vertices <- cbind(points, 0)
mesh <- cgalMesh$new(vertices = vertices, faces = d[["faces"]])
Edges <- mesh$getEdges()
Faces <- mesh$getFacesInfo()
Circumcenters <- Faces[, c("ccx", "ccy")]
Vertices <- points

voronoiEdgeFromDelaunayEdge <- function(edgeIndex) {
  faces <- unlist(Edges[edgeIndex, c("f1", "f2")])
  face1 <- faces[1L]
  face2 <- faces[2L]
  c1 <- Circumcenters[face1, ]
  if(is.na(face2)){
    return(NULL)
  }
  c2 <- Circumcenters[face2, ]
  if(isTRUE(all.equal(c1, c2))){
    return(NULL)
  }
  cbind(c1, c2)
}

vertexNeighborEdges <- function(vertexId) {
  i1 <- Edges[["i1"]]
  i2 <- Edges[["i2"]]
  which((i1 == vertexId) | (i2 == vertexId))
}

voronoiCell <- function(facetsQuotienter, edgeTransformer){
  function(tessellation, vertexId){
    delaunayEdges <- facetsQuotienter(
      unname(vertexNeighborEdges(vertexId))
    )
    voronoiEdges <- Filter(Negate(is.null), lapply(delaunayEdges, function(dedge){
      voronoiEdgeFromDelaunayEdge(dedge)
    }))
    cell <- edgeTransformer(voronoiEdges)
    attr(cell, "bounded") <- TRUE
    cell
  }
}

voronoi0 <- function(cellGetter, tessellation) {
  bounded <- logical(nrow(Vertices))
  L <- lapply(1L:nrow(Vertices), function(i) {
    cell <- cellGetter(tessellation, i)
    bounded[i] <<- attr(cell, "bounded")
    cell
  })
  out <- tessellation:::zip(Vertices, L)
  attr(out, "nbounded") <- sum(bounded)
  out
}

v <- voronoi0(voronoiCell(identity, identity), NULL)

plot(Vertices, type = "n")
for(i in seq_along(v)) {
  vi <- v[[i]]
  site <- vi$site
  points(rbind(site), pch = 19)
  cell <- vi$cell
  for(j in seq_along(cell)) {
    edge <- cell[[j]]
    A <- edge[, 1L]
    B <- edge[, 2L]
    segments(A[1L], A[2L], B[1L], B[2L], col = "black")
  }
}
