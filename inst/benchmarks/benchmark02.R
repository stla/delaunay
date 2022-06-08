library(microbenchmark)

n <- 600L
M <- matrix(rpois(2L*n, 10), nrow = n, ncol = 2L)

f1 <- function(){
  apply(M, 1L, sort)
}

f2 <- function(){
  col1 <- M[, 1L]
  col2 <- M[, 2L]
  cbind(pmin(col1, col2), pmax(col1, col2))
}

f3 <- function(){
  col1 <- M[, 1L]
  col2 <- M[, 2L]
  mx <- pmax(col1, col2)
  sm <- rowSums(M)
  cbind(sm - mx, mx)
}

microbenchmark(
  f1 = f1(),
  f2 = f2(),
  f3 = f3(),
  times = 50
)
