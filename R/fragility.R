#' @title compute fragility of ecological communities in gradual pressed context
#' @param sim.out output of simulation under gradual pressed conditions (\code{\link{sim_ode_press}})
#' @return resistance measured by the length of community trajectory
#' @return fragility measured by the variance of community trajectory
fragility <- function(sim.out) {
  # trajectory of survived species at each step
  trajectory = laply(sim.out, function(one) {
    length(one$species.survived)
  })
  # the area of the trajectory which reflects the resistance of system
  resistance.speciesnum = sum(trajectory) 
  # trajectory of NEW extinct species at each step
  trajectory.diff = trajectory[-length(trajectory)] - trajectory[-1]
  fragility.variance = sum(trajectory.diff^2)
  trajectory.positive = trajectory.diff[trajectory.diff > 0]
  fragility.entropy = sum(trajectory.positive * log(trajectory.positive))

  # add species abunance
  trajectory.abund <- laply(sim.out, function(one) {
    sum(one$nstar)
  })
  resistance.abund = sum(trajectory.abund)
  list(trajectory = trajectory, variance = fragility.variance, entropy = fragility.entropy, resistance = length(sim.out), resistance.speciesnum = resistance.speciesnum,
       resistance.abund = resistance.abund)
}

fragility.abund <- function(sim.out) {
  trajectory.abund <- laply(sim.out, function(one) {
    sum(one$nstar)
  })
  resistance = sum(trajectory.abund)
  list(resistance = resistance)
}
#' @title compute resistance of mutualistic communities in gradual pressed context
#' @param sim.out output of simulation under gradual pressed conditions (\code{\link{sim_ode_press}})
#' @return resistance measured by the length of community trajectory
resistance <- function(sim.out) {
  length(sim.out)
}


# random simulation times

