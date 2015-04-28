#####################Competition-Antagonism-Mutualism##########################

#' @title Lotka-Volterra (LV) Holling type II model for a competition-antagonism-mutualism hybrid community
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param params, parameters passed to LV model, a list of:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a matrix of intra-species and inter-species competitions}
#'   \item{A}{a matrix of antagonistic interactions among species}
#'   \item{M}{a matrix of mutualistic interactions among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#'   \item{g}{the conversion efficiency of antagonistic interactions}
#' }
#' @return the derivation
#' @import deSolve
model_lv2_cam <- function(time, init, params, ...) {
  r = params[[1]]  # intrinsic growth rates
  C = params[[2]]  # the competition matrix
  A = params[[3]]  # the antagonism matrix
  M = params[[4]]  # the mutualism matrix
  h = params[[5]]  # handling time
  g = params[[6]]  # conversion efficiency
  N = init  # initial state values (abundances) of speices
  # seperate the antagonistic matrix to two sub-matrice :
  # First sub-matrix describes the positive interactions among species
  # Second sub-matrix describes the negative interactions among species
  AP = A
  AP[AP < 0] = 0 # all the negative elements are replaced by 0 for positive sub-matrix
  AN = A
  AN[AN > 0] = 0 # all the positive elements are replaced by 0 for negative sub-matrix
  AN = abs(AN)
  H = matrix(h, nrow = length(h), ncol = length(h), byrow = T)
  # For the positive sub-matrix, 
  dN <- N * ( r - C %*% N  # competitive interactions 
              + (M %*% N) / (1 + h * M %*% N)  # mutualistic interactions
              + g * (AP %*% N) / (1 + h * AP %*% N)  # positive part of antagonistic interactions
              - (AN %*% N) / (1 + (H * AN) %*% N)  # negative part of antagonistic interactions
  )
  list(c(dN))  
}

##########################Competition-Antagonism##############################

#' @title Lotka-Volterra (LV) Equations of Holling type II for a community mixed by Competition and Antagonistic interactions
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param params, parameters passed to LV model, a list of:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a matrix of intra-species and inter-species competitions}
#'   \item{A}{a matrix of antagonistic interactions among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#'   \item{g}{the conversion efficiency of antagonistic interactions}
#' }
#' @return the derivation
#' @import deSolve
model_lv2_ca <- function(time, init, params, ...) {
  r = params[[1]]  # intrinsic growth rates
  C = params[[2]]  # the competitive matrix
  A = params[[3]]  # the antagonistic matrix
  h = params[[4]]  # handling time
  g = params[[5]]  # conversion efficiency
  N = init  # initial state
  # seperate the antagonistic matrix to two sub-matrice :
  # First sub-matrix describes the positive interactions among species
  # Second sub-matrix describes the negative interactions among species
  AP = A
  AP[AP < 0] = 0 # all the negative elements are replaced by 0 for positive sub-matrix
  AN = A
  AN[AN > 0] = 0 # all the positive elements are replaced by 0 for negative sub-matrix
  AN = abs(AN)
  H = matrix(h, nrow = length(h), ncol = length(h), byrow = T)
  # For the positive sub-matrix, 
  dN <- N * ( r - C %*% N + g * (AP %*% N) / (1 + h * AP %*% N) - (AN %*% N) / (1 + (H * AN) %*% N) )
  list(c(dN))
}

##########################Competition-Mutualism##############################

#' @title Lotka-Volterra (LV) Equations of Holling type II for a community mixed by Competition and Mutualism interactions
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param params, parameters passed to LV model, a list of:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a matrix of intra-species and inter-species competitions}
#'   \item{M}{a matrix of mutualism interactions among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#' }
#' @return the derivation
#' @import deSolve
model_lv2_cm <- function(time, init, params, ...) {
  r = params[[1]]  # intrinsic growth rates
  C = params[[2]]  # the competitive matrix
  M = params[[3]]  # the mutualistic matrix
  h = params[[4]]  # handling time
  N = init  # initial state
  dN <- N * ( r - C %*% N + (M %*% N) / (1 + h * M %*% N) )
  list(c(dN))
}

