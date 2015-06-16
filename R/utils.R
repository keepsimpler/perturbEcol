#' @title test if an object is null according to its type
is.empty <- function(obj) {
  if (is.null(obj)) {
    return(TRUE)
  }
  else if ( (is.numeric(obj) || is.list(obj)) && (length(obj) == 0) ) {
    return(TRUE)
  }
  else if (is.function(obj) && is.null(body(obj))) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}


#' @title transfer an incidence matrix to an adjacency matrix
#' @param inc, an incidence matrix
#' @return adj, an adiacency matrix
inc_to_adj <- function(inc){
  p <- dim(inc)[1]  # number of Plants
  a <- dim(inc)[2]  # number of Animals
  s <- p + a  # number of all Species
  adj <- matrix(0, s, s)  # initialize the adjacency matrix as a zero-matrix
  adj[1:p, (p + 1):s] <- inc  # the upper right sub-matrix is the incidence matrix
  adj <- adj + t(adj)  # the lower left sub-matrix is transpose of the incidence matrix
  return(adj)
}


###############################################################################
#' @title Generate a connected graph using package [igraph]
#'
#' @param s, size of network. 
#' if graph type is bipartite, s[1], s[2] represent size of two groups; else s is size of network
#' @param k, average degree for the network.
#' 1 < k < s for unipartite network, 1 < k < s[1]*s[2]/(s[1]+s[2]) for bipartite network.
#' @param gtype, Graph type generated: 'bipartite', 'sf', 'er', 'regular'.
#' @param maxtried, the maximum number of tried times. 
#' If have tried [maxtried] times, the function will return no matter whether the connected graph is generated.
#' @param expower exponent coefficient of Scale-Free network
#' @param ... the params conform to the implementation functions of [igraph]
#' @return the connected graph
#' @details .  
#' @import igraph
gen_connected_graph <- function(s, k, gtype, maxtried = 100, expower = 2.5, ...) {
  #library(igraph)
  if (gtype == 'bipartite' && is.na(s[2])) {  # the bipartite graph need size of two groups of nodes
    warning('sizes of TWO groups of nodes should be designated. 
            we have assumed the size of second group equal to the size of first group.')
    s[2] = s[1]  # if missed second size, we assume it equal to the first size.
  }
  count = 0
  repeat {  # generate a connected graph
    if (gtype == 'bipartite') {
      G = bipartite.random.game(s[1], s[2], type = 'gnm', m = ceiling(k * (s[1] + s[2])))
    } else if (gtype == 'sf') {
      G = static.power.law.game(s, k * s, exponent.out = expower)
    }
    else if (gtype == 'er') {
      G = erdos.renyi.game(s, p.or.m = k * s, type = 'gnm')
    }
    else if (gtype == 'regular') {
      G = k.regular.game(s, k)
    }
    else if (gtype == 'complete') {
      G = graph.full(s)
    }
    if (igraph::is.connected(G)) break  # until a connected graph is generated
    count = count + 1
    if (count == maxtried) {
      warning(paste('Tried', maxtried, 'times, But connected graph still cannot be generated.'))
      break
    }
  }
  G
}
#plot(G, layout = layout.bipartite)

#' @title generate different type of food webs
#' @param s number of species
#' @param k average interactions (feeding and being feeded) of species
#' @param type different type of food webs generated:
#' \describe{
#'   \item{full}{a full food web, where every species interact with each other, and feed or being feeded randomly}
#'   \item{bipartite}{a bipartite food web, where species are seperated to two groups, the first group have all preys, the second group have all predators.}
#'   \item{niche}{a food web generated using niche model}
#' }
gen_foodweb <- function(s, k, type = 'full') {
  if (type == 'bipartite' && is.na(s[2])) {  # the bipartite graph need size of two groups of nodes
    warning('sizes of TWO groups of nodes should be designated. 
            we have assumed the size of second group equal to the size of first group.')
    s[2] = s[1]  # if missed second size, we assume it equal to the first size.
  }
  
  if (type == 'full') {
    G = graph.full(s)
    foodweb = as.matrix(get.adjacency(G))
    for (i in 2:s) {
      for (j in 1:(i-1)) {
        if (foodweb[i, j] == 1) { # if an undirected edge exists between i and j
          if (runif(1) < 0.5)
            foodweb[i, j] = -1
          else
            foodweb[j, i] = -1
        }
      }
    }
  }
  else if (type == 'bipartite') {
    G = bipartite.random.game(s[1], s[2], type = 'gnm', m = ceiling(k * (s[1] + s[2]) / 2))
    foodweb = as.matrix(get.adjacency(G))
    # the predators feed on preys and have negative effect on preys
    foodweb[upper.tri(foodweb)] = - foodweb[upper.tri(foodweb)]
  }
  foodweb
}

#' @title generate a hybrid network that include competition, antagonism, mutualism interactions
#' @param s number of species
#' @param k average degree of species
#' @param type network type, 'er':random graph, 'sf':scale-free, 'bipartite':bipartite graph, 'niche':niche model
#' @param pc probability of competition interactions
#' @param pa probability of antagonism interactions
#' @param pm probability of mutualism interactions
#' @param ... additional arguments transformed to graph generate such as [expower]
gen_hybrid_network <- function(s, k, type = 'er', pc = 0., pa = 0., pm = 1., ...) {
  stopifnot(pc >= 0., pa >= 0., pm >= 0., pc + pa + pm == 1)
  G = gen_connected_graph(s, k, type, ...)  # generate a connected graph
  graph = as.matrix(get.adjacency(G))  # transform to matrix form
  # split the graph to three sub-graph - competition, antagonism and mutualism graphs according to the probability of occurance of three different types of interactions
  competitive_graph = matrix(0, nrow = s, ncol = s)
  antago_graph = matrix(0, nrow = s, ncol = s)
  mutual_graph = matrix(0, nrow = s, ncol = s)
  ps = runif(sum(graph))
  cursor = 0
  for (i in 2:s) {
    for (j in 1:(i-1)) {
      if (graph[i, j] == 1) { # if an undirected edge exists between i and j
        cursor = cursor + 1
        p = ps[cursor]
        if (p < pc) {
          competitive_graph[i, j] = 1
          competitive_graph[j, i] = 1
        }
        else if (p < pc + pa / 2) {
          antago_graph[i, j] = 1
          antago_graph[j, i] = - 1
        }
        else if (p < pc + pa) {
          antago_graph[i, j] = - 1
          antago_graph[j, i] = 1
        }
        else if ( p < pc + pa + pm) {
          mutual_graph[i, j] = 1
          mutual_graph[j, i] = 1
        }
      }
    }
  }
  list(competitive_graph = competitive_graph, antago_graph = antago_graph, mutual_graph = mutual_graph)
}
#' @title a niche model food web generator according to Williamns and Martinez nature 2000
#' copy from https://gist.github.com/emhart/1503428
#' @param s, # of species
#' @param c, connectivity
niche_model <- function(s, c) {
  new.mat <- matrix(0, nrow = s, ncol = s) # initialize a zero matrix
  ci <- vector()
  niche <- runif(s, 0, 1)
  r <- rbeta(s, 1, ( ( 1 / (2 * c) ) - 1 )) * niche
  for(i in 1:s) {
    ci[i] <- runif(1, r[i] / 2, niche[i])
  }
  #now set the smallest species niche value to have an n of 0
  r[which(niche==min(niche))] <- .00000001
  for(i in 1:s){
    for(j in 1:s){
      if(niche[j] > (ci[i] - .5 * r[i]) &&
           niche[j] < (ci[i] + .5 * r[i])) {
        new.mat[j,i] <- 1
      }
    }
  }
  new.mat <- new.mat[, order( apply(new.mat, 2, sum) )]
  return(new.mat)
}


#' @title another form of uniform distribution between [mean - sd, mean + sd]
runif2 <- function(n, mean, sd) {
  runif(n) * 2 * sd + (mean - sd)
}

display_output <- function(sim.out, coeff, fragility) {
  ode.nstars = laply(sim.out, function(one) {
    one$nstar
  })
  # fragility = fragility(sim.out)
  
  fig_title <- paste('h.mu = ', coeff[['h.mu']], ', gamma.mu = ', coeff[['gamma.mu']], ', entropy = ', sprintf("%.2f", fragility['entropy']), ', variance = ', fragility['variance'], ', resist = ', fragility['resistance'], sep = '')
  matplot(ode.nstars, main = fig_title, type = 'l', lwd = 1.5, xlab = 'steps of gradual presses', ylab = 'Abundances of species')  
}