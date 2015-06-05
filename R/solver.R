setGeneric("sim", function(obj, initialize = TRUE, ...) standardGeneric("sim"))

setMethod("sim", "perturbModel",
          function(obj, initialize = TRUE, ...) {
            if (initialize) obj <- initialize(obj)
            sim.out <- do.call(obj@solver, list(obj@main, obj@params, obj@init, obj@times, obj@perturb, obj@perturbNum, ...))
            #obj@out <- obj@solver(obj@main, obj@params, obj@init, obj@times, obj@perturb, obj@perturbNum, ...)
            fragility <- fragility(sim.out)
            obj@out <- list(sim.out = sim.out, fragility = fragility)
            invisible(obj)
          }
)


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
#'   \item{extinct.species}{a vector of species that extincted}
#' }
sim_ode_press <- function(model, params, init, times, perturb, perturbNum = 500, isout = TRUE, extinct_threshold, ...) {
  #times = seq(from = 0, to = steps * stepwise, by = stepwise)
  ode.outs = list()
  for(i in 1:perturbNum) {
    #print(i)
    ode.out = ode(init, times, model, params) 
    nstar = as.numeric(ode.out[nrow(ode.out), 2:ncol(ode.out)]) # species biomass at equilibrium
    nstar[nstar < extinct_threshold] = 0  # species with biomass less than extinct threshold is considered to be extinct
    extinct.species = which(nstar == 0)  # extinct species
    
    flag = 0
    # if all species are extinct, will end the simulation
    if (length(nstar) == length(extinct.species)) flag = 1
    # if any species' abundance is NaN, that means the ODE dynamic is unstable, the simulation will also be ended
    if (any(is.nan(nstar))) flag = 2
    
    Phi = jacobian.full(y = nstar, func = model, params = params) # community matrix, Jacobian matrix at equilibrium
    if (isout) {
      ret = list(out = ode.out, nstar = nstar, Phi = Phi, params = params, extinct.species = extinct.species, flag = flag)
    }
    else {
      ret = list(nstar = nstar, Phi = Phi, params = params, extinct.species = extinct.species, flag = flag)
    }
    ode.outs[[length(ode.outs) + 1]] = ret
    # if all species are extinct, end the simulation
    if (flag == 1 || flag == 2)
      break;
    
    #     if (length(extinct.species) > 0) {
    #       ret = remove.species(params, nstar, extinct.species)
    #       params = ret$params
    #       nstar = ret$nstar
    #     }
    #     if (length(nstar) == 0) break  # if all species are extinct, then stop and quit
    
    perturb.res = perturb(params, nstar, ...)
    params = perturb.res$params
    init = perturb.res$nstar
  }
  return(ode.outs)
}

