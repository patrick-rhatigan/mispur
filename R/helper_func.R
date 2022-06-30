# Positive part/negative part functions
# Take positive part
pos <- function(x) { #[ ]_+
  (x > 0) * x
}

# Take negative part
neg <- function(x) { #[ ]_-
  x[x == Inf] <- 0
  x <- - (x < 0) * x
}

# Treshold at 1
pmax1 <- function(x) {
  x[x < 1] <- 1
  return(x)
}

# EGMS functions
phi <- function(xi) {
  xi[xi > 1] <- Inf
  xi[xi <= 1] <- 0

  return(xi)
}

chi <- function(nu, c) {
  neg(nu + c) - neg(c)
}
