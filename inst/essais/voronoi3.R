library(delaunay)
library(cgalMeshes)

nsides <- 12L
angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
outer_points <- cbind(cos(angles), sin(angles))
inner_points <- outer_points / 4

nsides <- 36L
angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
middle_points <- cbind(cos(angles), sin(angles)) / 2
points <- rbind(outer_points, inner_points, middle_points)

angles <- angles + pi/36
middle_points <- cbind(cos(angles), sin(angles)) / 3
points <- rbind(points, middle_points)
middle_points <- cbind(cos(angles), sin(angles)) / 1.5
points <- rbind(points, middle_points)


# constraint edges
indices <- 1L:12L
edges <- cbind(
  indices, c(indices[-1L], indices[1L])
)
edges <- rbind(edges, edges + 12L)
# constrained Delaunay triangulation
d <- delaunay(points, constraints = edges)

#####

# polygon(outer_points, lwd = 6, border = "black")
# polygon(inner_points, lwd = 6, border = "black")

#####
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
    return(cbind(c1, NA))
    # if(Edges[edgeIndex, "length"] < 0.5) {
    #   return(cbind(c1, c(0,0))) # NA  
    # } else {
    #   M <- (Vertices[Edges[edgeIndex, "i1"], ] + Vertices[Edges[edgeIndex, "i2"], ])/2
    #   return(cbind(c1, 2*M))
    # }
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
    nedges <- length(voronoiEdges)
    bounded <- nedges > 0L
    while(nedges > 0L){
      bounded <- bounded && !anyNA(edges[[nedges]])
      nedges <- nedges - 1L
    }
    cell <- edgeTransformer(voronoiEdges)
    attr(cell, "bounded") <- bounded
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


library(randomcoloR)
opar <- par(mar = c(0,0,0,0))
plot(Vertices*1, type = "n", asp = 1, axes=FALSE, xlab=NA, ylab=NA)
#plotDelaunay2D(d, asp = 1)
for(i in seq_along(v)) {
  vi <- v[[i]]
  site <- vi$site
  cell <- vi$cell
  vertices <- unique(t(do.call(cbind, cell)))
  if(anyNA(vertices)) {
    next
  }
  vectors <- sweep(vertices, 2L, colMeans(vertices), check.margin = FALSE)
  vector1 <- vectors[1L, ]
  a <- atan2(vector1[2L], vector1[1L])
  vectors <- vectors[-1L, ]
  angles <- c(0, apply(vectors, 1L, function(v) atan2(v[2L], v[1L]) - a))
  vertices <- vertices[order(angles %% (2*pi)), ]
  polygon(vertices, col = randomColor(hue="random", luminosity = "dark"))
  points(rbind(site), pch = 19)
  # for(j in seq_along(cell)) {
  #   edge <- cell[[j]]
  #   A <- edge[, 1L]
  #   B <- edge[, 2L]
  #   segments(A[1L], A[2L], B[1L], B[2L], col = "black", lwd=6)
  # }
}
