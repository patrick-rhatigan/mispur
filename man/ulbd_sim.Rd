% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_sim.R
\name{ulbd_sim}
\alias{ulbd_sim}
\title{Runs simulations for the generic upper bound/lower bound models}
\usage{
ulbd_sim(Nsim = 1000, n = 250, B = 1000, sdS = 250, k1, k2, mu, ncores)
}
\arguments{
\item{k1}{This gives the vector k1 which corresponds to
the moment "lower bounds"}

\item{k2}{This gives the vector k2 which corresponds to
the moment "upper bounds"}

\item{mu}{This gives the mean vector}

\item{...}{These arguments are passed to \code{\link{momineq}}}
}
\description{
Here, rather than providing r or length, the user provides the mean vector
and then r and length are calculated
}
