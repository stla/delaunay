voronoiEdgeFromDelaunayEdge <- function(edgeIndex, Edges, Circumcenters) {
  faces <- unlist(Edges[edgeIndex, c("f1", "f2")])
  face1 <- faces[1L]
  face2 <- faces[2L]
  c1 <- Circumcenters[face1, ]
  if(is.na(face2)){
    return(cbind(c1, NA_real_))
  }
  c2 <- Circumcenters[face2, ]
  if(isTRUE(all.equal(c1, c2))) {
    return(NULL)
  }
  cbind(c1, c2)
}

vertexNeighborEdges <- function(vertexId, Edges) {
  i1 <- Edges[["i1"]]
  i2 <- Edges[["i2"]]
  which((i1 == vertexId) | (i2 == vertexId))
}

VoronoiCell <- function(mesh, facetsQuotienter, edgeTransformer) {
  Edges <- mesh[["edges"]]
  Faces <- mesh[["faces"]]
  Circumcenters <- Faces[, c("ccx", "ccy")]
  function(vertexId) {
    delaunayEdges <- facetsQuotienter(
      unname(vertexNeighborEdges(vertexId, Edges))
    )
    voronoiEdges <- Filter(
      Negate(is.null), lapply(delaunayEdges, function(dedge) {
        voronoiEdgeFromDelaunayEdge(dedge, Edges, Circumcenters)
      })
    )
    nedges <- length(voronoiEdges)
    bounded <- nedges > 0L
    while(nedges > 0L) {
      bounded <- bounded && !anyNA(voronoiEdges[[nedges]])
      nedges <- nedges - 1L
    }
    cell <- edgeTransformer(voronoiEdges)
    attr(cell, "bounded") <- bounded
    cell
  }
}

#' @importFrom sets pair
#' @noRd
zip <- function(matrix, list) {
  lapply(seq_len(nrow(matrix)), function(i) {
    pair(site = matrix[i, ], cell = list[[i]])
  })
}

Voronoi0 <- function(cellGetter, mesh) {
  vertices <- mesh[["vertices"]]
  nvertices <- nrow(vertices)
  bounded <- logical(nvertices)
  L <- lapply(1L:nvertices, function(i) {
    cell <- cellGetter(i)
    bounded[i] <<- attr(cell, "bounded")
    cell
  })
  out <- zip(vertices, L)
  attr(out, "nbounded") <- sum(bounded)
  out
}

#' @title Voronoï diagram
#' @description Returns the Voronoï tessellation corresponding to a 2D 
#'   Delaunay triangulation. If the Delaunay triangulation is constrained, 
#'   the output is not a true constrained Voronoï tessellation.
#' 
#' @param triangulation a 2D Delaunay triangulation obtained with the 
#'   \code{\link{delaunay}} function
#'
#' @return The Voronoï diagram given as list of pairs; each pair is made of 
#'   one site (a vertex of the Delaunay triangulation) and one Voronoï cell, 
#'   a polygon given as numeric matrix with two columns (that can be plotted 
#'   with \code{\link[graphics:polygon]{polygon}}). 
#' @export
#' @seealso \code{\link{plotVoronoi}}
#'
#' @examples
#' library(delaunay)
#' # make the vertices
#' nsides <- 6L
#' angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
#' outer_points <- cbind(cos(angles), sin(angles))
#' inner_points <- outer_points / 4
#' nsides <- 12L
#' angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
#' middle_points <- cbind(cos(angles), sin(angles)) / 2
#' points <- rbind(outer_points, inner_points, middle_points)
#' angles <- angles + pi/24
#' middle_points <- cbind(cos(angles), sin(angles)) / 3
#' points <- rbind(points, middle_points)
#' middle_points <- cbind(cos(angles), sin(angles)) / 1.5
#' points <- rbind(points, middle_points)
#' # constraint edges
#' indices <- 1L:6L
#' edges <- cbind(
#'   indices, c(indices[-1L], indices[1L])
#' )
#' edges <- rbind(edges, edges + 6L)
#' ## | constrained Delaunay triangulation 
#' del <- delaunay(points, constraints = edges)
#' opar <- par(mar = c(0,0,0,0))
#' plotDelaunay2D(
#'   del, type = "n", xlab = NA, ylab = NA, axes = FALSE, asp = 1,
#'   fillcolor = "random", luminosity = "dark", 
#'   col_borders = "black", lwd_borders = 3
#' )
#' par(opar)
#' ## | corresponding Voronoï diagram
#' vor <- Voronoi(del)
#' opar <- par(mar = c(0,0,0,0))
#' plot(
#'   NULL, asp = 1, axes = FALSE, xlab = NA, ylab = NA, 
#'   xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5) 
#' )
#' plotVoronoi(vor, luminosity = "dark")
#' points(points, pch = 19)
#' par(opar)
Voronoi <- function(triangulation) {
  if(!inherits(triangulation, "delaunay")) {
    stop(
      "The argument `triangulation` must be an output of the `delaunay` function.",
      call. = TRUE
    )
  }
  dimension <- attr(triangulation, "dimension")
  if(dimension != 2) {
    stop(
      sprintf("Invalid dimension (%s instead of 2).", dimension),
      call. = TRUE
    )
  }
  mesh <- triangulation[["mesh"]]
  Voronoi0(VoronoiCell(mesh, identity, identity), mesh)  
}

#' @title Plot Voronoï diagram
#' @description Plot a Voronoï tessellation.
#'
#' @param tessellation an output of \code{\link{Voronoi}}
#' @param colors this can be \code{"random"} to use random colors for the cells
#'   (with \code{\link[randomcoloR]{randomColor}}), \code{"distinct"} to use
#'   distinct colors with the help of
#'   \code{\link[randomcoloR]{distinctColorPalette}}, \code{NA} for no colors, 
#'   or this can be a vector of colors; the length of this vector
#'   of colors must match the number of Voronoï cells, that you can get by
#'   typing \code{length(tessellation)}
#' @param hue,luminosity if \code{colors = "random"}, these arguments are passed
#'   to \code{\link[randomcoloR]{randomColor}}
#' @param alpha opacity, a number between 0 and 1
#'   (used when \code{colors} is not \code{NA})
#' @param ... arguments passed to \code{\link[graphics]{polygon}} to plot the 
#'   cells
#'
#' @return No returned value.
#' @export
#' @importFrom graphics polygon
#' @importFrom randomcoloR randomColor distinctColorPalette
#' @importFrom scales alpha
#'
#' @examples 
#' library(delaunay)
#' # make vertices
#' nsides <- 12L
#' angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
#' outer_points <- cbind(cos(angles), sin(angles))
#' inner_points <- outer_points / 4
#' nsides <- 36L
#' angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
#' middle_points <- cbind(cos(angles), sin(angles)) / 2
#' vertices <- rbind(outer_points, inner_points, middle_points)
#' angles <- angles + pi/36
#' middle_points <- cbind(cos(angles), sin(angles)) / 3
#' vertices <- rbind(vertices, middle_points, 2*middle_points)
#' # constraint edges
#' indices <- 1L:12L
#' edges <- cbind(
#'   indices, c(indices[-1L], indices[1L])
#' )
#' edges <- rbind(edges, edges + 12L)
#' ## | constrained Delaunay triangulation 
#' del <- delaunay(vertices, constraints = edges)
#' opar <- par(mar = c(0,0,0,0))
#' plotDelaunay2D(
#'   del, type = "n", xlab = NA, ylab = NA, axes = FALSE, asp = 1,
#'   fillcolor = "random", luminosity = "dark", 
#'   col_borders = "black", lwd_borders = 3
#' )
#' par(opar)
#' ## | corresponding Voronoï diagram
#' vor <- Voronoi(del)
#' opar <- par(mar = c(0,0,0,0))
#' plot(
#'   vertices, type = "n", asp = 1, axes = FALSE, xlab = NA, ylab = NA 
#' )
#' plotVoronoi(vor, luminosity = "dark")
#' points(vertices, pch = 19)
#' par(opar)
plotVoronoi <- function(
    tessellation, colors = "random", hue = "random", luminosity = "light", 
    alpha = 1, ...
){
  ncells <- length(tessellation)
  if(identical(colors, "random")) {
    colors <- scales::alpha(
      randomColor(ncells, hue = hue, luminosity = luminosity), alpha
    )
  } else if(identical(colors, "distinct")) {
    colors <- scales::alpha(distinctColorPalette(ncells), alpha)
  } else if(identical(colors, NA)) {
    colors <- rep(NA, ncells)
  } else {
    if(length(colors) != ncells) {
      stop(
        sprintf(
          "There are %d Vorono\u00ef cells and you supplied %d colors.",
          ncells, length(colors)
        )
      )
    }
    colors <- scales::alpha(colors, alpha)
  }
  for(i in 1L:ncells) {
    cell <- tessellation[[i]][["cell"]]
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
    polygon(vertices, col = colors[i], ...)
  }
}
