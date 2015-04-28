############################################################
#  The S4 classes of perturbEcol                           #
############################################################

## root class of perturbed model
setClass("perturbModel",
         representation(
           main = "function",
           times     = "numeric",
           init      = "numeric",
           params     = "list",
           perturb   = "function",
           perturbNum= "numeric",
           solver    = "function",
           out       = "ANY"
         )
)

## class of perturbed model for competition and mutualism mixed commnities
setClass("cmModel",
         representation(
           out       = "list"
           ),
         prototype(
           main = model_lv2_cm
         ),
         contains    = "perturbModel"
)

## class of perturbed model for competition, antagonism and mutualism mixed commnities
setClass("camModel",
         representation(
           out       = "list"
         ),
         prototype(
           main = model_lv2_cam
         ),
         contains    = "perturbModel"
)
