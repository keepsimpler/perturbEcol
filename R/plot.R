## plot method for cmModel class
setMethod("plot", c("perturbModel", "missing"),
          function(x, y, ...) {
            if (is.null(x@out) || length(x@out) == 0)
              stop("Please simulate the model before plotting", call. = FALSE)
            sim.out = x@out
            fragility = x@fragility
            ode.nstars = laply(sim.out, function(one) {
              one$nstar
            })
            par(mfrow=c(2,1)) 
            matplot(ode.nstars, type = 'l', lwd = 1., xlab = 'Steps of gradual pressures', ylab = 'Abundances of species', ...)
            matplot(fragility$trajectory, type = 'l', lwd = 1., col = c(3), xlab = 'Steps of gradual pressures', ylab = 'Number of survived species', ...)
            
          }
)

