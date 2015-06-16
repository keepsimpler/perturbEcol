setGeneric("sim", function(obj, initialize = TRUE, ...) standardGeneric("sim"))

setMethod("sim", "perturbModel",
          function(obj, initialize = TRUE, ...) {
            if (initialize) obj <- initialize(obj)
            sim.out <- do.call(obj@solver, list(obj@main, obj@params, obj@init, obj@times, obj@perturb, obj@perturbNum, ...))
            #obj@out <- obj@solver(obj@main, obj@params, obj@init, obj@times, obj@perturb, obj@perturbNum, ...)
            fragility <- fragility(sim.out)
            obj@out <- sim.out
            obj@fragility <- fragility
            invisible(obj)
          }
)

