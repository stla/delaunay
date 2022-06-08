distance <- function(A, B){
  sqrt(c(crossprod(A-B)))
}

triangleArea <- function(A, B, C){
  a <- distance(B, C)
  b <- distance(A, C)
  c <- distance(A, B)
  s <- (a + b + c) / 2
  sqrt(s*(s-a)*(s-b)*(s-c))
}

delaunayArea <- function(vertices, triangles){
  ntriangles <- nrow(triangles)
  areas <- numeric(ntriangles)
  for(i in 1L:ntriangles){
    points <- vertices[triangles[i, ], ]
    areas[i] <- triangleArea(points[1L, ], points[2L, ], points[3L, ])
  }
  sum(areas)
}
