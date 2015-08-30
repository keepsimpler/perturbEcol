## accessor and replacement functions for slots of perturbModel
setGeneric("main", function(obj, ...) standardGeneric("main"))
setGeneric("main<-", function(obj, value) standardGeneric("main<-"))

setMethod("main", "perturbModel",
          function(obj, ...) {obj@main}
)

setMethod("main<-", "perturbModel", 
          function(obj, value) {
            obj@main <- value
            obj@out <- list()
            invisible(obj)
          }
)

setGeneric("times", function(obj, ...) standardGeneric("times"))
setGeneric("times<-", function(obj, value) standardGeneric("times<-"))

setMethod("times", "perturbModel",
          function(obj, ...) {obj@times}
)

setMethod("times<-", "perturbModel", 
          function(obj, value) {
            obj@times <- value
            obj@out <- list()
            invisible(obj)
          }
)

setGeneric("init", function(obj, ...) standardGeneric("init"))
setGeneric("init<-", function(obj, value) standardGeneric("init<-"))

setMethod("init", "perturbModel",
          function(obj, ...) {obj@init}
)

setMethod("init<-", "perturbModel", 
          function(obj, value) {
            obj@init <- value
            obj@out <- list()
            invisible(obj)
          }
)

setGeneric("params", function(obj, ...) standardGeneric("params"))
setGeneric("params<-", function(obj, value) standardGeneric("params<-"))

setMethod("params", "perturbModel",
          function(obj, ...) {obj@params}
)

setMethod("params<-", "perturbModel", 
          function(obj, value) {
            obj@params <- value
            obj@out <- list()
            invisible(obj)
          }
)

setGeneric("perturb", function(obj, ...) standardGeneric("perturb"))
setGeneric("perturb<-", function(obj, value) standardGeneric("perturb<-"))

setMethod("perturb", "perturbModel",
          function(obj, ...) {obj@perturb}
)

setMethod("perturb<-", "perturbModel", 
          function(obj, value) {
            obj@perturb <- value
            obj@out <- list()
            invisible(obj)
          }
)

setGeneric("perturbNum", function(obj, ...) standardGeneric("perturbNum"))
setGeneric("perturbNum<-", function(obj, value) standardGeneric("perturbNum<-"))

setMethod("perturbNum", "perturbModel",
          function(obj, ...) {obj@perturbNum}
)

setMethod("perturbNum<-", "perturbModel", 
          function(obj, value) {
            obj@perturbNum <- value
            obj@out <- list()
            invisible(obj)
          }
)

setGeneric("solver", function(obj, ...) standardGeneric("solver"))
setGeneric("solver<-", function(obj, value) standardGeneric("solver<-"))

setMethod("solver", "perturbModel",
          function(obj, ...) {obj@solver}
)

setMethod("solver<-", "perturbModel", 
          function(obj, value) {
            obj@solver <- value
            obj@out <- list()
            invisible(obj)
          }
)

setGeneric("out", function(obj, ...) standardGeneric("out"))
setGeneric("out<-", function(obj, value) standardGeneric("out<-"))

setMethod("out", "perturbModel",
          function(obj, ...) {obj@out}
)

setMethod("out<-", "perturbModel", 
          function(obj, value) {
            obj@out <- list()
            invisible(obj)
          }
)

#' @title parameters for consumer-resource model according to the hybrid network and the coefficients
#' @param hybrid_graph the hybrid interaction topology of communities, which includes three sub-graphs: competition, antagonism and mutualism
#' @param coeff, a list of coefficients:
#' \describe{
#'    \item{alpha.mu, alpha.sd}{coefficients of the intrinsic growth rates of species}
#'    \item{beta0.mu, beta0.sd}{the intra-species competition coefficients which determin a uniform distribution in [beta.mu - beta.sd, beta.mu + beta.sd]}
#'    \item{antago.mu, antago.sd}{the inter-species antagonism coefficients}
#'    \item{h.mu, h.sd}{coefficients of the handling time of species}
#'    \item{g.mu, g.sd}{conversion efficiency of antagonistic interactions from consumer side}
#'    \item{e.mu, e.sd}{conversion efficiency of antagonistic interactions from resource side}
#' }
#' @return a list of parameters for consumer-resource model:
#' \describe{
#'   \item{r}{a vector, the intrinsic growth rates of species}
#'   \item{C}{a vector, the intraspecies self-regulation of species}
#'   \item{M}{a matrix, consumer-resource interactions among species}
#'   \item{H}{a matrix, handling time of consumer species i on resource species j}
#'   \item{G}{a matrix, conversion rate of resource j to the gain of abundance of species i}
#'   \item{E}{a matrix, conversion rate of resource j for species i to the lost of abundance of species j}
#' }
params_cr2 <- function(hybrid_graph, coeff) {
  competitive_graph = hybrid_graph$competitive_graph
  antago_graph = hybrid_graph$antago_graph
  mutual_graph = hybrid_graph$mutual_graph
  stopifnot(dim(antago_graph)[1] == dim(antago_graph)[2], dim(competitive_graph)[1] == dim(competitive_graph)[2], dim(antago_graph)[1] == dim(competitive_graph)[1])
  s = dim(antago_graph)[1]
  with(as.list(coeff), {
    r <- runif2(s, alpha.mu, alpha.sd)
    C <- runif2(s, beta0.mu, beta0.sd) # assign intra-species competitive interaction strengths
    A <- params_antago_interactions(antago_graph, antago.mu, antago.sd)
    M <- params_mutual_interactions(mutual_graph, antago.mu, antago.sd, 0)
    #A[A < 0] = 0  # delete negative links
    #M <- A + M  # combine antagonism and mutualism interactions
    H = matrix(runif2(s * s, h.mu, h.sd), nrow = s, ncol = s)
    G = matrix(runif2(s * s, g.mu, g.sd), nrow = s, ncol = s)
    E = matrix(runif2(s * s, e.mu, e.sd), nrow = s, ncol = s)
    list(r = r, C = C, A = A, M = M, H = H, G = G, E = E)     
  })
}

params_cr2_4 <- function(hybrid_graph, coeff) {
  mutual_graph = hybrid_graph$mutual_graph
  s = dim(mutual_graph)[1]
  with(as.list(coeff), {
    r <- runif2(s, alpha.mu, alpha.sd)
    C <- runif2(s, beta0.mu, beta0.sd) # assign intra-species competitive interaction strengths
    
    # mutual- interaction strengths are symmetric
    tmp = matrix(rep(0, s * s), nrow = s)
    tmp[upper.tri(tmp)] = runif2(s * (s - 1) / 2, antago.mu, antago.sd)
    tmp = tmp + t(tmp)
    M = tmp * mutual_graph
    
    A = matrix(0, nrow = s, ncol = s)

    H = matrix(runif2(s * s, h.mu, h.sd), nrow = s, ncol = s)
    
    shapes = betaShapes(mean1, var1, mean2, var2)
    X = rBivarBetas(s * s, shapes['alpha1'], shapes['beta1'], shapes['alpha2'], shapes['beta2'], rho = rho)
    G = matrix(X[, 1], nrow = s, ncol = s) * mutual_graph
    E = matrix(X[, 2], nrow = s, ncol = s) * mutual_graph
    
    list(r = r, C = C, A = A, M = M, H = H, G = G, E = E)     
  })
}
params_cr2_3 <- function(hybrid_graph, coeff) {
  antago_graph = hybrid_graph$antago_graph
  mutual_graph = hybrid_graph$mutual_graph
  stopifnot(dim(antago_graph)[1] == dim(antago_graph)[2], dim(mutual_graph)[1] == dim(mutual_graph)[2])
  s = dim(antago_graph)[1]
  A = antago_graph
  M = mutual_graph
  with(as.list(coeff), {
    r <- runif2(s, alpha.mu, alpha.sd)
    C <- runif2(s, beta0.mu, beta0.sd) # assign intra-species competitive interaction strengths
    A[A > 0] = runif2(sum(A > 0), antago.mu, antago.sd)  # assign inter-species antagonism interaction strengths
    # mutual- interaction strengths are symmetric
    tmp = matrix(rep(0, s * s), nrow = s)
    tmp[upper.tri(tmp)] = runif2(s * (s - 1) / 2, antago.mu, antago.sd)
    tmp = tmp + t(tmp)
    M = tmp * mutual_graph
    
    H = matrix(runif2(s * s, h.mu, h.sd), nrow = s, ncol = s)
    
    # assign consumer-side conversion rates ([g]) and
    # resource-side conversion rates ([e]) to
    # mutual- interactions
    # for mutual- interactions, 1 >= [g] > [e] >0
    # so we need to generate uniform random points in a tringle:
    # A--(0,0)  B--(1,0)  C--(1,1) for ([g],[e])
    # tringle = matrix(c(0, 0, 1, 0, 1, 1), nrow = 3, byrow = T)
    if (rho.mutual < 0) {
      rho.mutual = 1 + rho.mutual
      tringle1 = matrix(c(0, 0, 1, 1 - rho.mutual, 1, 1), nrow = 3, byrow = T)
      tringle2 = matrix(c(0, 0, rho.mutual, 0, 1, 1 - rho.mutual), nrow = 3, byrow = T)
      P1 = runifTringle(tringle1, floor(s * s / (1 + 1 - rho.mutual)))
      P2 = runifTringle(tringle2, s * s - floor(s * s / (1 + 1 - rho.mutual)))
      P = rbind(P1, P2)
      stopifnot(nrow(P) == s * s)
    } else {
      tringle = matrix(c(rho.mutual, 0, 1, 0, 1, 1 - rho.mutual), nrow = 3, byrow = T)
      P = runifTringle(tringle, s * s)
    }
    G.mutual = matrix(P[, 1], ncol = s)
    E.mutual = matrix(P[, 2], ncol = s)
    #E.mutual = t(E.mutual) # trick with t(E * M) in model_cr2
    G.mutual = G.mutual * mutual_graph
    E.mutual = E.mutual * mutual_graph
    
    # assign [g] and [e] to antago- interactions
    # for antago- interactions, 1 > [g] > 0, [e] = 1
    E.antago = antago_graph
    if (rho.antago < 0) {
      rho.antago = - rho.antago
      G.antago = matrix(runif(s * s, min = 1 - rho.antago, max = 1), ncol = s) * antago_graph
    } else {
      G.antago = matrix(runif(s * s, min = 0, max = rho.antago), ncol = s) * antago_graph
    }
    
    E = E.mutual + E.antago
    G = G.mutual + G.antago
    
    list(r = r, C = C, A = A, M = M, H = H, G = G, E = E)     
  })
}

params_cr2_2 <- function(hybrid_graph, coeff) {
  antago_graph = hybrid_graph$antago_graph
  mutual_graph = hybrid_graph$mutual_graph
  stopifnot(dim(antago_graph)[1] == dim(antago_graph)[2], dim(mutual_graph)[1] == dim(mutual_graph)[2])
  s = dim(antago_graph)[1]
  A = antago_graph
  M = mutual_graph
  with(as.list(coeff), {
    r <- runif2(s, alpha.mu, alpha.sd)
    C <- runif2(s, beta0.mu, beta0.sd) # assign intra-species competitive interaction strengths
    A[A > 0] = runif2(sum(A > 0), antago.mu, antago.sd)  # assign inter-species antagonism interaction strengths
    M[M > 0] = runif2(sum(M > 0), antago.mu, antago.sd)  # assign inter-species mutualism interaction strengths
    H = matrix(runif2(s * s, h.mu, h.sd), nrow = s, ncol = s)
    
    # number of mutual- and antago- edges
    edges.antago = sum(antago_graph > 0)
    edges.mutual = sum(mutual_graph > 0)
    edges.total = edges.antago + edges.mutual
    
    # assign consumer-side conversion rates of interactions. 
    # [g] is divided according to interaction types (antago- or mutual-).
    # for antago- interactions, the consumer-side conversion rates is argued to be very small, we fix it to be 0.1 mean and 0.1 relative sd.
    g.antago.mu = 0.1
    g.sd = 0.1
    # the mean of consumer-side conversion rates for mutual- interactions is assigned under the constraint(restriction) that keeping the mean of consumer-side conversion rates for all interactions constant
    if (edges.mutual == 0) {
      g.mutual.mu = 0
    }
    else {
      g.mutual.mu = (edges.total * g.mu - edges.antago * g.antago.mu) / edges.mutual         
    }
    if (g.mutual.mu <= 0) return(list())
    G.mutual = matrix(runif2(s * s, g.mutual.mu, g.mutual.mu * g.sd), nrow = s, ncol = s) * mutual_graph
    G.antago = matrix(runif2(s * s, g.antago.mu, g.antago.mu * g.sd), nrow = s, ncol = s) * antago_graph
    G = G.mutual + G.antago
    
    # assign resource-side conversion rates of interactions. 
    # [e] is divided according to interaction types (antago- or mutual-).
    # for antago- interactions, the resource-side conversion rates is argued to be 1.0.
    e.entago.mu = 1.
    e.sd = 0.1
    if (edges.mutual == 0) {
      e.mutual.mu = 0
    }
    else {
      e.mutual.mu = (edges.total * e.mu - edges.antago * e.antago.mu) / edges.mutual         
    }
    # if the resource-side conversion rates of mutual- interactions is 0, then fail
    if (e.mutual.mu <= 0) return(list())
    E.mutual = matrix(runif2(s * s, e.mutual.mu, e.mutual.mu * e.sd), nrow = s, ncol = s) * mutual_graph
    E.antago = matrix(runif2(s * s, e.antago.mu, e.antago.mu * 0), nrow = s, ncol = s) * antago_graph
    E = E.mutual + E.antago
    
    list(r = r, C = C, A = A, M = M, H = H, G = G, E = E)     
  })
}


#' @title parameters for consumer-resource model according to the hybrid network and the coefficients
#' @param hybrid_graph the hybrid interaction topology of communities, which includes three sub-graphs: competition, antagonism and mutualism
#' @param coeff, a list of coefficients:
#' \describe{
#'    \item{alpha.mu, alpha.sd}{coefficients of the intrinsic growth rates of species}
#'    \item{beta0.mu, beta0.sd}{the intra-species competition coefficients which determin a uniform distribution in [beta.mu - beta.sd, beta.mu + beta.sd]}
#'    \item{antago.mu, antago.sd}{the inter-species antagonism coefficients}
#'    \item{h.mu, h.sd}{coefficients of the handling time of species}
#'    \item{g.mu, g.sd}{conversion efficiency of antagonistic interactions}
#' }
#' @return a list of parameters for consumer-resource model:
#' \describe{
#'   \item{r}{a vector, the intrinsic growth rates of species}
#'   \item{C}{a vector, the intraspecies self-regulation of species}
#'   \item{M}{a matrix, consumer-resource interactions among species}
#'   \item{H}{a matrix, handling time of consumer species i on resource species j}
#'   \item{G}{a matrix, conversion rate of resource j to the gain of abundance of species i}
#' }
params_cr <- function(hybrid_graph, coeff) {
  competitive_graph = hybrid_graph$competitive_graph
  antago_graph = hybrid_graph$antago_graph
  mutual_graph = hybrid_graph$mutual_graph
  stopifnot(dim(antago_graph)[1] == dim(antago_graph)[2], dim(competitive_graph)[1] == dim(competitive_graph)[2], dim(antago_graph)[1] == dim(competitive_graph)[1])
  s = dim(antago_graph)[1]
  with(as.list(coeff), {
    r <- runif2(s, alpha.mu, alpha.sd)
    C <- runif2(s, beta0.mu, beta0.sd) # assign intra-species competitive interaction strengths
    A <- params_antago_interactions(antago_graph, antago.mu, antago.sd)
    M <- params_mutual_interactions(mutual_graph, antago.mu, antago.sd, 0)
    #A[A < 0] = 0  # delete negative links
    #M <- A + M  # combine antagonism and mutualism interactions
    H = matrix(runif2(s * s, h.mu, h.sd), nrow = s, ncol = s)
    G = matrix(runif2(s * s, g.mu, g.sd), nrow = s, ncol = s)
    list(r = r, C = C, A = A, M = M, H = H, G = G)     
  })
}


#' @title parameters for hybrid LV2 model according to the network and the coefficients
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
  if (new_total_strength != 0)
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

