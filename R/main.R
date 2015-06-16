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


#' @title Simulate ODE dynamics of non-autonomous systems.A example is ecosystems under "press" perturbations. The dynamic is iteration of successive ODE dynamics of automous sytems (\code{\link{sim_ode_auto}}), while at each iterating step, the parameters and/or state values of systems are changed to reflect "press" perturbations.
#' @inheritParams sim_ode_auto
#' @param perturb a function that change the parameters and state values after each iteration step
#' @param iter_steps possiblely maximum iteration steps
#' @param isout if output the transiting trajectory of each ODE iterate step
#' @param ... any arguments which are transfered to perturbation function
#' @return a list of lists :
#' \describe{
#'   \item{out}{output of one ODE simulation, including the trajectory of values of state variables}
#'   \item{nstar}{the values of state variables in equilibrium}
#'   \item{Phi}{the Jacobian matrix in equilibrium}
#'   \item{params}{parameters assigned to the model}
#'   \item{species.survived}{a vector of species that survived}
#' }
sim_ode_press <- function(model, params, init, times, perturb, perturbNum = 500, isout = FALSE, extinct_threshold = 1e-10, ...) {
  ode.outs = list()
  for(i in 1:perturbNum) {
    ode.out = ode(init, times, model, params) 
    nstar = as.numeric(ode.out[nrow(ode.out), 2:ncol(ode.out)]) # species biomass at equilibrium
    nstar[nstar < extinct_threshold] = 0  # species with biomass less than extinct threshold is considered to be extinct
    species.survived = which(nstar > 0)  # survived species
    
    flag = 0
    # if all species are extinct, will end the simulation
    if (length(species.survived) == 0) flag = 1
    # if any species' abundance is NaN, that means the ODE dynamic is unstable, the simulation will also be ended
    if (any(is.nan(nstar))) flag = 2
    
    Phi = jacobian.full(y = nstar, func = model, params = params) # community matrix, Jacobian matrix at equilibrium
    if (isout) {
      ret = list(out = ode.out, nstar = nstar, Phi = Phi, params = params, species.survived = species.survived, flag = flag)
    }
    else {
      ret = list(nstar = nstar, Phi = Phi, params = params, species.survived = species.survived, flag = flag)
    }
    ode.outs[[length(ode.outs) + 1]] = ret
    # if all species are extinct, end the simulation
    if (flag == 1 || flag == 2)
      break;
    
    # perturbation that returns new parameters and initial values
    perturb.ret = perturb(params, nstar, ...)
    params = perturb.ret$params
    init = perturb.ret$nstar
  }
  return(ode.outs)
}


#     if (length(extinct.species) > 0) {
#       ret = remove.species(params, nstar, extinct.species)
#       params = ret$params
#       nstar = ret$nstar
#     }
#     if (length(nstar) == 0) break  # if all species are extinct, then stop and quit
