#' @title compute fragility of ecological communities in gradual pressed context
#' @param sim.out output of simulation under gradual pressed conditions (\code{\link{sim_ode_press}})
#' @return resistance measured by the length of community trajectory
#' @return fragility measured by the variance of community trajectory
fragility <- function(sim.out) {
  trajectory = laply(sim.out, function(one) {
    length(one$extinct.species)
  })
  resistance2 = sum(trajectory) # the complement area of the trajectory
  trajectory = trajectory[-1] - trajectory[-length(trajectory)]
  fragility.variance = sum(trajectory^2)
  trajectory.positive = trajectory[trajectory > 0]
  fragility.entropy = sum(trajectory.positive * log(trajectory.positive))

  trajectory.abund <- laply(sim.out, function(one) {
    sum(one$nstar)
  })
  resistance.abund = sum(trajectory.abund)
  list(trajectory = trajectory, variance = fragility.variance, entropy = fragility.entropy, resistance = length(sim.out), resistance.abund = resistance.abund, resistance2 = resistance2)
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

