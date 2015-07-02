############################################################
#  The S4 classes of perturbEcol                           #
############################################################

## root class of perturbed model
setClass("perturbModel",
         representation(   # define the slots (variables and functions) of class
           main = "function",
           times     = "numeric",
           init      = "numeric",
           params     = "list",
           perturb   = "function",
           perturbNum= "numeric",
           solver    = "function",
           out       = "ANY",
           fragility = "ANY"
         ),
         prototype(  # initial values of the slots
           solver = sim_ode_press
         )
)

## class of perturbed model for competition and mutualism mixed commnities
setClass("cmModel",
         representation(
           out       = "list",
           fragility = "list"
           ),
         prototype(  # initial values of the slots
           main = model_lv2_cm
         ),
         contains    = "perturbModel"  # inherit from super-class
)

## class of perturbed model for competition, antagonism and mutualism mixed commnities
setClass("camModel",
         representation(
           out       = "list",
           fragility = "list"
         ),
         prototype(
           main = model_lv2_cam
         ),
         contains    = "perturbModel"
)
