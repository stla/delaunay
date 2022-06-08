#' Title
#'
#' @param points xx
#' @param edges xx
#'
#' @return xx
#' @export
#'
#' @importFrom Rvcg vcgGetEdge
#' 
#' @examples
#' nsides <- 12L
#' angles <- seq(0, 2*pi, length.out = nsides+1L)[-1L]
#' points <- cbind(cos(angles), sin(angles))
#' points <- rbind(points, points/1.5)
#' # constraint edges
#' indices <- 1L:nsides
#' edges_outer <- cbind(
#'   indices, c(indices[-1L], indices[1L])
#' )
#' edges_inner <- edges_outer + nsides
#' edges <- rbind(edges_outer, edges_inner)
#' cdel(points, edges)
cdel <- function(points, edges){
  if(!is.matrix(points) || !is.numeric(points)){
    stop(
      "The `points` argument must be a numeric matrix with two or three columns.", 
      call. = TRUE
    )
  }
  dimension <- ncol(points)
  if(!is.element(dimension, c(2L, 3L))){
    stop(
      "The `points` argument must be a numeric matrix with two or three columns.", 
      call. = TRUE
    )
  }
  if(any(is.na(points))){
    stop("Missing values are not allowed.", call. = TRUE)
  }
  if(anyDuplicated(points)){
    stop("There are some duplicated points.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  tpoints <- t(points)
  if(!is.matrix(edges) || !is.numeric(edges) || ncol(edges) != 2L){
    stop(
      "The `edges` argument must be an integer matrix with two columns.", 
      call. = TRUE
    )
  }
  if(any(is.na(edges))){
    stop("Missing values in `edges` are not allowed.", call. = TRUE)
  }
  storage.mode(edges) <- "integer"
  stopifnot(all(edges >= 1L))
  stopifnot(all(edges <= nrow(points)))
  edges <- t(apply(edges, 1L, sort))
  if(anyDuplicated(edges)){
    stop("There are some duplicated constraint edges.", call. = TRUE)
  }
  if(any(edges[, 1L] == edges[, 2L])){
    stop("There are some invalid constraint edges.", call. = TRUE)
  }
  del <- del2d_constrained_cpp(tpoints, t(edges))
  rglfake <- list(vb = rbind(tpoints, 1), it = del)
  class(rglfake) <- "mesh3d"
  edges <- `colnames<-`(
    as.matrix(vcgGetEdge(rglfake))[, -3L], c("v1", "v2", "border")
  )
  edges
}