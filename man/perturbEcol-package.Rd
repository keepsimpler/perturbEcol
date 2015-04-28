\name{perturbEcol-package}
\alias{perturbEcol-package}
\alias{perturbEcol}
\docType{package}
\title{
perturbEcol: Simulation of ecological dynamic systems under environmental perturbations
}
\description{
\strong{perturbEcol} simulates the dynamics of networked ecological systems that undergo environmental perturbations such as climate change etc. 

The simulation model is implemented as S4 class with the following slots:
  \itemize{
    \item \code{main = function(time, init, params, ...):} a function holding the main equations of the model, such as ODE, SDE.
    \item \code{times:} vector of time steps or vector with three named values \code{from, to, by} specifying the simulation time steps.
    \item \code{init:} initial state (start values) of the simulation. This is a named vector of state variables such as species abundances.
    \item \code{params:} a list of model parameters whose values would be changed by iterated perturbations.
    \item \code{perturb = function(params, ...):} a function describing(?) the perturbation or pressure by environments on ecological systems. The results of perturbations are the changing of model parameters.
    \item \code{perturbNum:} the number of repeated perturbations.
    \item \code{solver = function(main, times, init, params, perturb, perturbNum, ...):} a function specifying the numerical algorithm of simulation. The algorithms model the dynamics of ecological systems undergoing environmental perturbations.
    \item \code{out:} the output of simulation.
  }
}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
\seealso{
  \code{\link[<pkg>:<pkg>-package]{<perturbEcol>}}
}
\examples{
~~ simple examples of the most important functions ~~
}
