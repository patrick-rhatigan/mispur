#Simple helper functions that will be used throughout the package


#' Positive Elements Only
#'
#' Keeps only the positive elements from the vector and replaces
#'  all negative elements with zeroes.
#'
#' @param x A vector of numerics
#'
#' @return A new vector where the negative values in x are
#'  replaced by zeros
#' @seealso [neg()]
#'
#'
#' @examples pos(-2:2) == c(0,0,0,1,2)
pos <- function(x) { #[ ]_+
  (x > 0) * x
}


#' neg
#'
#' @param x vector of numerics
#'
#' @return vector x with positive values replaced by zeros
#'
#'
#' @examples neg(-2:2) == c(-2,-1,0,0,0)
neg <- function(x) { #[ ]_-
  x[x == Inf] <- 0
  x <- - (x < 0) * x
}


#' pmax1
#'
#' @param x vector of numerics
#'
#' @return vector x where any value > 1 is replaced by 1
#'
#'
#' @examples pmax1(-2:2) == c(-2,-1,0,1,1)
pmax1 <- function(x) {
  x[x < 1] <- 1
  return(x)
}

# EGMS functions


#' phi
#'
#' @param xi vector of numerics
#'
#' @return vector where any xi>1 is infinity and otherwise is zero
#'
#'
#' @examples
#' phi(-2:2) == c(0,0,0,0,Inf)
phi <- function(xi) {
  xi[xi > 1] <- Inf
  xi[xi <= 1] <- 0

  return(xi)
}

chi <- function(nu, c) {
  neg(nu + c) - neg(c)
}
