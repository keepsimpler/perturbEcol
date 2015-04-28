## accessor and replacement functions for params slot of perturbModel
setGeneric("params", function(obj, ...) standardGeneric("params"))
setGeneric("params<-", function(obj, value) standardGeneric("params<-"))

setMethod("params", "perturbModel",
          function(obj, ...) {obj@params}
)

setMethod("params<-", "perturbModel", 
          function(obj, value) {
            obj@params <- value
            obj@out <- NULL
            invisible(obj)
          }
)

#' @title parmaters for hybrid LV2 model according to the network and the coefficients
#' @param hybrid_graph the hybrid interaction topology of communities, which includes three sub-graphs: competition, antagonism and mutualism
#' @param coeff, a list of coefficients:
#' \describe{
#'    \item{alpha.mu, alpha.sd}{coefficients of the intrinsic growth rates of species}
#'    \item{beta0.mu, beta0.sd}{the intra-species competition coefficients which determin a uniform distribution in [beta.mu - beta.sd, beta.mu + beta.sd]}
#'    \item{beta1.mu, beta1.sd}{the inter-species competition coefficients}
#'    \item{antago.mu, antago.sd}{the inter-species antagonism coefficients}
#'    \item{gamma.mu, gamma.sd}{the inter-species mutualism coefficients}
#'    \item{g.mu, g.sd}{conversion efficiency of antagonistic interactions}
#'    \item{h.mu, h.sd}{coefficients of the handling time of species}
#'    \item{delta}{trade-off coefficients of mutualistic interaction strengths}
#' }
#' @return a list of parameters for ode model:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a matrix of intra-species and inter-species competitions}
#'   \item{A}{a matrix of antagonistic interactions among species}
#'   \item{M}{a matrix of mutualistic interactions among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#'   \item{g}{the conversion efficiency of antagonistic interactions}
#' }
params_lv2_cam <- function(hybrid_graph, coeff) {
  competitive_graph = hybrid_graph$competitive_graph
  antago_graph = hybrid_graph$antago_graph
  mutual_graph = hybrid_graph$mutual_graph
  stopifnot(dim(antago_graph)[1] == dim(antago_graph)[2], dim(competitive_graph)[1] == dim(competitive_graph)[2], dim(antago_graph)[1] == dim(competitive_graph)[1])
  s = dim(antago_graph)[1]
  with(as.list(coeff), {
    C <- params_competitive_interactions(competitive_graph, beta0.mu, beta0.sd, beta1.mu, beta1.sd) # generate a competitive interaction matrix
    A <- params_antago_interactions(antago_graph, antago.mu, antago.sd)
    M <- params_mutual_interactions(mutual_graph, gamma.mu, gamma.sd, delta)
    h = runif2(s, h.mu, h.sd)
    g = runif2(s, g.mu, g.sd)
    r <- runif2(s, alpha.mu, alpha.sd)
    list(r = r, C = C, A = A, M = M, h = h, g = g)     
  })
}

#' @title generate the antagonistic interaction matrix according to the antagonistic network and coefficients
#' @param graph the antagonistic interaction topology of communities, which is a bi-directional network
#' @param gamma.mu
#' @param gamma.sd coefficients that determin a uniform distribution of antagonistic interaction strengths
params_antago_interactions <- function(graph, gamma.mu, gamma.sd) {
  stopifnot(dim(graph)[1] == dim(graph)[2]) # if not a adjacency matrix, stop
  stopifnot(gamma.mu >= 0, gamma.sd >= 0, gamma.mu >= gamma.sd)
  diag(graph) <- 0  # ensure the diagonal elements of antagonistic matrix equal 0
  edges = sum(graph > 0)  # !leak! number of interactions, number of positive and negative interactions are same with each other
  foodweb = graph
  foodweb[lower.tri(foodweb) & foodweb != 0] = foodweb[lower.tri(foodweb) & foodweb != 0] * runif(edges, min = gamma.mu - gamma.sd, max = gamma.mu + gamma.sd) # lower tringle of matrix and not equal zero are assigned random variable values
  foodweb[upper.tri(foodweb)] = - t(foodweb)[upper.tri(t(foodweb))] # upper tringle of matrix are negative of lower tringle
  foodweb
}

#' @title parmaters for antagonism LV2 model according to the network and the coefficients
#' @param antago_graph the antagonistic interaction topology of communities, which is the adjacency matrix of a network
#' @param competitive_graph the competitive interaction topology of communities
#' @param coeff, a list of coefficients:
#' \describe{
#'    \item{basal.alpha.mu, basal.alpha.sd}{coefficients of the intrinsic growth rates of basal species, which should be positive}
#'    \item{nobasal.alpha.mu, nobasal.alpha.sd}{coefficients of the intrinsic growth rates of not-basal species, which should be negative}
#'    \item{beta0.mu, beta0.sd}{the intra-species competition coefficients which determin a uniform distribution in [beta.mu - beta.sd, beta.mu + beta.sd]}
#'    \item{beta1.mu, beta1.sd}{the inter-species competition coefficients}
#'    \item{gamma.mu, gamma.sd}{the inter-species antagonism coefficients}
#'    \item{g.mu, g.sd}{conversion efficiency of antagonistic interactions}
#'    \item{h.mu, h.sd}{coefficients of the handling time of species}
#' }
#' @return a list of parameters for ode model:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a matrix of intra-species and inter-species competitions}
#'   \item{A}{a matrix of antagonistic interactions among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#'   \item{g}{the conversion efficiency of antagonistic interactions}
#' }
params_lv2_ca <- function(antago_graph, competitive_graph, coeff) {
  stopifnot(dim(antago_graph)[1] == dim(antago_graph)[2], dim(competitive_graph)[1] == dim(competitive_graph)[2], dim(antago_graph)[1] == dim(competitive_graph)[1])
  s = dim(antago_graph)[1]
  with(as.list(coeff), {
    r <- params_intrinsic_growth_rates(antago_graph, basal.alpha.mu, basal.alpha.sd, nobasal.alpha.mu, nobasal.alpha.sd)
    C <- params_competitive_interactions(competitive_graph, beta0.mu, beta0.sd, beta1.mu, beta1.sd) # generate a competitive interaction matrix
    A <- params_antago_interactions(antago_graph, gamma.mu, gamma.sd)
    h = runif2(s, h.mu, h.sd)
    g = runif2(s, g.mu, g.sd)
    list(r = r, C = C, A = A, h = h, g = g)     
  })
}

#' @title generate the mutualistic interaction matrix according to the mutualistic network and coefficients
#' @param graph the mutualistic interaction topology of communities, which is the adjacency matrix of a network
#' @param gamma.mu
#' @param gamma.sd coefficients that determin a uniform distribution of mutualistic interaction strengths
#' @param delta coefficient that determin the trade-off between the interaction strength and width(node degree) of species
params_mutual_interactions <- function(graph, gamma.mu, gamma.sd, delta) {
  stopifnot(dim(graph)[1] == dim(graph)[2]) # if not a adjacency matrix, stop
  stopifnot(gamma.mu >= 0, gamma.sd >= 0, gamma.mu >= gamma.sd, delta >= 0)
  diag(graph) <- 0  # ensure the diagonal elements of mutualistic matrix equal 0
  edges = sum(graph > 0)  # number of all interactions
  M = graph
  degrees = rowSums(M)  # the width (node degree) of species
  M[M > 0] = runif2(edges, gamma.mu, gamma.sd)  # assign inter-species mutualistic interaction strengths
  old_total_strength = sum(M)
  ## !leak!
  M = M / degrees^delta  # trade-off of mutualistic strengths
  new_total_strength = sum(M) # ensure the total strength constant before and after trade-off
  M = M * old_total_strength / new_total_strength
  M = t(M)  ## !leak! ??
  M
}

#' @title generate the competitive interaction matrix according to the competitive network and coefficients
#' @param graph the competitive interaction topology of communities, which is the adjacency matrix of a network
#' @param beta0.mu
#' @param beta0.sd coefficients that determin a uniform distribution of intra-species interaction strengths
#' @param beta1.mu
#' @param beta1.sd coefficients that determin a uniform distribution of inter-species interaction strengths
params_competitive_interactions <- function(graph, beta0.mu, beta0.sd, beta1.mu, beta1.sd) {
  stopifnot(dim(graph)[1] == dim(graph)[2]) # if not a adjacency matrix, stop
  stopifnot(beta0.mu > 0, beta0.sd >= 0, beta0.mu >= beta0.sd, beta1.mu >= 0, beta1.sd >= 0, beta1.mu >= beta1.sd)
  diag(graph) <- 0
  edges = sum(graph > 0)  # number of all interactions
  s <- dim(graph)[1] # number of total Species
  C <- graph
  C[C > 0] = runif2(edges, beta1.mu, beta1.sd) # assign inter-species competitive interaction strengths
  diag(C) <- runif2(s, beta0.mu, beta0.sd) # assign intra-species competitive interaction strengths
  C
}

#' @title parmaters for mutualism LV2 model according to the network and the coefficients
#' @param mutual_graph the mutualistic interaction topology of communities, which is the adjacency matrix of a network
#' @param competitive_graph the competitive interaction topology of communities
#' @param coeff, a list of coefficients:
#' \describe{
#'    \item{alpha.mu, alpha.sd}{coefficients of the intrinsic growth rates of species}
#'    \item{beta0.mu, beta0.sd}{the intra-species competition coefficients which determin a uniform distribution in [beta.mu - beta.sd, beta.mu + beta.sd]}
#'    \item{beta1.mu, beta1.sd}{the inter-species competition coefficients}
#'    \item{gamma.mu, gamma.sd}{the inter-species mutualism coefficients}
#'    \item{delta}{trade-off coefficients of mutualistic interaction strengths}
#'    \item{h.mu, h.sd}{coefficients of the handling time of species}
#' }
#' @return a list of parameters for ode model:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a competitive interaction matrix }
#'   \item{M}{a mutualistic interaction matrix among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#' }
params_lv2_cm <- function(mutual_graph, competitive_graph, coeff) {
  stopifnot(dim(mutual_graph)[1] == dim(mutual_graph)[2], dim(competitive_graph)[1] == dim(competitive_graph)[2], dim(mutual_graph)[1] == dim(competitive_graph)[1])
  s = dim(mutual_graph)[1]
  with(as.list(coeff), {
    r <- runif2(s, alpha.mu, alpha.sd)
    C <- params_competitive_interactions(competitive_graph, beta0.mu, beta0.sd, beta1.mu, beta1.sd) # generate a competitive interaction matrix
    M <- params_mutual_interactions(mutual_graph, gamma.mu, gamma.sd, delta)
    h = runif2(s, h.mu, h.sd)
    list(r = r, C = C, M = M, h = h)     
  })
}

