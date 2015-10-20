#' @title initial values of state variables, i.e., abundances of species
#' @description Assign initial values according to two criteria: 1. using the equilibrium values of LV1 model as initial values. 2. If any of the initial values is less than 0, using the average of equilibrium values as initial values.
#' @param params, the parameters assigned to LV2 model
init_cr2 <- function(params) {
  params$A[params$A < 0] = 0
  init = solve(diag(params$C) - params$G * (params$A + params$M) + t(params$E * (params$A + params$M))) %*% params$r
  if (any(init < 0)) {
    warning('Initial state values is less than 0 !!')
    init = abs(params$r)
  }
  c(init)  
}


#' @title initial values of state variables, i.e., abundances of species
#' @description Assign initial values according to two criteria: 1. using the equilibrium values of LV1 model as initial values. 2. If any of the initial values is less than 0, using the average of equilibrium values as initial values.
#' @param params, the parameters assigned to LV2 model
init_lv2_cam <- function(params) {
  init = solve(params$C - params$A - params$M) %*% params$r
  if (any(init < 0)) {
    warning('Initial state values is less than 0 !!')
    init = params$r
  }
  c(init)  
}


#' @title initial values of state variables, i.e., abundances of species
#' @description Assign initial values according to two criteria: 1. using the equilibrium values of LV1 model as initial values. 2. If any of the initial values is less than 0, using the average of equilibrium values as initial values.
#' @param params, the parameters assigned to LV2 model
init_lv2_ca <- function(params) {
  init = solve(params$C - params$A) %*% params$r
  if (any(init < 0)) {
    warning('Initial state values is less than 0 !!')
    s = length(params$r)
    init = rep(sum(init) / s, s)
  }
  c(init)  
}

#' @title initial values of state variables, i.e., abundances of species
#' @description Assign initial values according to two criteria: 1. using the equilibrium values of LV1 model as initial values. 2. If any of the initial values is less than 0, using the intrinsic growth rates as initial values.
#' @param params, the parameters assigned to LV2 model
init_lv2_cm <- function(params) {
  init = solve(params$C - params$M) %*% params$r
  if (any(init < 0)) {
    warning('Initial state values is less than 0 !!')
    #stop('Initial state values is less than 0 !!', init(LV2))
    init = params$r
  }
  c(init)  
}

