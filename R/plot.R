## plot method for cmModel class
setMethod("plot", c("perturbModel", "missing"),
          function(x, y, ...) {
            if (is.null(x@out))
              stop("Please simulate the model before plotting", call. = FALSE)
            sim.out = x@out$sim.out
            fragility = x@out$fragility
            ode.nstars = laply(sim.out, function(one) {
              one$nstar
            })
            
            matplot(ode.nstars, type = 'l', lwd = 1.5, xlab = 'steps of gradual presses', ylab = 'Abundances of species', ...)  
            
          }
)

