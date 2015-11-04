#' @title get proportions of competition, antagonism, mutualism interactions accroding to consumer-side and resource-side conversion rates
#' @param params, model params including all conversion rates
proportions <- function(params) {
  G = params$G
  E = params$E
  GE = G - E
  GE[GE > 0] = 1
  GE[GE < 0] = -2
  GE = GE + t(GE)
  pc = sum(GE == -4) / sum(GE != 0)
  pa = sum(GE == -1) / sum(GE != 0)
  pm = sum(GE == 2) / sum(GE != 0)
  c(pc = pc, pa = pa, pm = pm)
}
#' @title generate polygon distribution in rectangle [(0,0), (1,1)]
#' @param rho, control sample subspace in rectangle [(0,0), (1,1)]
#' @param n, number of random variables
polygon_distribution <- function(rho, n) {
  if (rho < 0) {
    rho = - rho
    tringle1 = matrix(c(0, 0, 1, 0, 1, 1), nrow = 3, byrow = T)
    tringle2 = matrix(c(0, 0, 1, 1, 1 - rho, 1), nrow = 3, byrow = T)
    tringle3 = matrix(c(0, 0, 1 - rho, 1, 0, rho), nrow = 3, byrow = T)
    P1 = runifTringle(tringle1, floor(n / (1 + rho + rho * (1 - rho))))
    P2 = runifTringle(tringle2, floor(rho * n / (1 + rho + rho * (1 - rho))))
    P3 = runifTringle(tringle3, n - floor(n / (1 + rho + rho * (1 - rho))) -
                        floor(rho * n / (1 + rho + rho * (1 - rho))) )
    P = rbind(P1, P2, P3)
    stopifnot(nrow(P) == n)
  } else {
    tringle = matrix(c(rho, 0, 1, 0, 1, 1 - rho), nrow = 3, byrow = T)
    P = runifTringle(tringle, n)
  }
  P
}

#' @title generate tringle distribution in [(0,0), (1,0), (1,1)]
#' @param rho, control the sample space
#' @param n, number of random variables
tringle_distribution <- function(rho, n) {
  if (rho < 0) {
    rho = 1 + rho
    tringle1 = matrix(c(0, 0, 1, 1 - rho, 1, 1), nrow = 3, byrow = T)
    tringle2 = matrix(c(0, 0, rho, 0, 1, 1 - rho), nrow = 3, byrow = T)
    P1 = runifTringle(tringle1, floor(n / (1 + 1 - rho)))
    P2 = runifTringle(tringle2, n - floor(n / (1 + 1 - rho)))
    P = rbind(P1, P2)
    stopifnot(nrow(P) == n)
  } else {
    tringle = matrix(c(rho, 0, 1, 0, 1, 1 - rho), nrow = 3, byrow = T)
    P = runifTringle(tringle, n)
  }
  P
  #c(apply(P, 2, mean), apply(P, 2, var))
}


# a new kind of matrix multiplication?
# it differ from matrix multiplication based on inner products
# and Kronecker multiplication based on outer products.
# it expand two matrix A(m*n) and B(n*p) to new C(m*n*p)
# A = matrix(c(0,1,1,0,1,0,0,1,1,0,0,0,0,1,0,0), nrow = 4, byrow = T)
multiply.matrix <- function(A, B) {
  stopifnot(dim(A)[2] == dim(B)[1])
  m = dim(A)[1]
  n = dim(A)[2]
  p = dim(B)[2]
  C = array(0, dim = c(m, n, p))
  for (i in 1:m) {
    for (j in 1:p) {
      for (k in 1:n) {
        C[i, k, j] = A[i, k] * B[k, j]
      }
    }
  }
  C
}

#' @title get abundance sum of length 2 interactions
#' @param C, arrary(i,j,k), i-->j-->k interactions
#' @param N, species abundances
abund.length2 <- function(C, N) {
  apply(C, c(1,2), function(x) sum(x * N))
}

abund.length <- function(A, N) {
  m = dim(A)[1]
  n = dim(A)[2]
  N2 = matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    for (j in 1:n) {
      if (A[i, j] == 1)
        N2[i, j] = sum(A[j, ] * N)
    }
  }
  N2
}
# repeat a vector [n] times to create a matrix
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

# http://www.evanmiller.org/bayesian-ab-testing.html#cite1   Probability of one Beta variable is charger than another
#' @title generate two correlated Beta random variables
#' @param n, number of samples
#' @param mean1, var1, mean2, var2
#' @param rho, correlation coefficient
#' 
rBivarBetas <- function(n, alpha1, beta1, alpha2, beta2, rho) {
  # simulate bivariate normal data with the specified correlation
  #require(MASS) # use [mvnorm] function
  Z = mvrnorm(n, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2), empirical = T) # empirical: samples or population
  # transform the normal variates into uniform variates
  U = apply(Z, 2, pnorm)
  X1 = qbeta(U[,1], alpha1, beta1)
  X2 = qbeta(U[,2], alpha2, beta2)
  X = cbind(X1, X2)
  X
}

#' @title calculate shape parameters alpha and beta from means and variances
#' @param mean1, var1, mean2, var2
#' @return a list of (alpha1, beta1, alpha2, beta2)
betaShapes <- function(mean1, var1, mean2, var2) {
  stopifnot(mean1 > 0 & mean1 < 1 & mean2 > 0 & mean2 < 1 &
              var1 > 0 & var1 <= mean1 * (1 - mean1) &
              var2 > 0 & var2 <= mean2 * (1 - mean2))
  alpha1 = (mean1 * (1 - mean1) / var1 - 1) * mean1
  alpha2 = (mean2 * (1 - mean2) / var2 - 1) * mean2
  beta1 = alpha1 * (1 / mean1 - 1)
  beta2 = alpha2 * (1 / mean2 - 1)
  c(alpha1 = alpha1, beta1 = beta1, alpha2 = alpha2, beta2 = beta2)
}


#' @title generate uniformly random points within a tringle
#' @description \deqn{P = (1 - \sqrt{r_1}) A + (\sqrt{r_1} (1 - r_2))  B + (r_2 \sqrt{r_1}) C}
#'              where \deqn{r_1, r_2 \sim U[0, 1]}
#' @references http://www.cs.princeton.edu/~funk/tog02.pdf (Section 4.2)
#' @param tringle, a tringle represented by three points, a 3*2 matrix
#' @param n, number of uniformly random points 
#' @examples
#'  tringle = matrix(c(0, 0, 1, 0, 1, 1), nrow = 3, byrow = T)
#'  unifTringle(tringle, 10)
runifTringle <- function(tringle, n) {
  A = tringle[1, ]
  B = tringle[2, ]
  C = tringle[3, ]
  r1 = runif(n) # uniform in [0, 1]
  r2 = runif(n) # uniform in [0, 1]
  P = (1 - r1 ^ 0.5) %*% t(A) + (r1 ^ 0.5 * (1 - r2)) %*% t(B) +
    (r2 * r1 ^ 0.5) %*% t(C)
  colnames(P) <- c('x', 'y')
  P
}


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
#' @param gtype, Graph type generated: 'bipartite', 'sf', 'er', 'dag', 'regular'.
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
      G = sample_gnm(s, k * s)
    }
    else if (gtype == 'regular') {
      G = k.regular.game(s, k)
    }
    else if (gtype == 'complete') {
      G = graph.full(s)
    }
    else if (gtype == 'dag') {
      require('spacejam')  # generate random directed Acyclic graphs
      G = rdag(s, s * k * 2)
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

#' @title Rewiring one link from mutualism to antagonism
#' @param A, the antagonism sub-graph
#' @param M, the mutualism sub-graph
#' @details pick one mutualism interaction randomly, then split it to two antagonism interactions
rewire_mutual_to_antago <- function(A, M) {
  A[A < 0] = 0 # remove negative elements of antagonism sub-graph
  s = dim(A)[1]  # number of species
  if (sum(M != 0) > 0) { # if have mutualism interactions to split
    # pick one mutualism interaction randomly, then split it to two antagonism interactions
    indx = which(M > 0, arr.ind = T)
    indx.one = indx[sample(1:nrow(indx), 1), ]
    # split one mutualism interaction to two antagonism interactions
    # one is left at the same place, the one is random
    A[indx.one[1], indx.one[2]] = M[indx.one[1], indx.one[2]]
    repeat {
      newlink = sample(1:s, 2)
      if (A[newlink[1], newlink[2]] == 0 || A[newlink[2], newlink[1]] == 0 || M[newlink[1], newlink[2]] == 0 || M[newlink[2], newlink[1]] == 0) {  # if there is no existed interaction at the place of new link, then move half of mutual interaction to here
        A[newlink[1], newlink[2]] = M[indx.one[2], indx.one[1]]
        break
      }
    }
    M[indx.one[1], indx.one[2]] = 0
    M[indx.one[2], indx.one[1]] = 0
  }
  else {
    warning("No more mutualism interaction available!")
  }
  list(A = A, M = M)
}

#' @title Rewiring one link from antagonism to mutualism
#' @param A, the antagonism sub-graph
#' @param M, the mutualism sub-graph
#' @param need.connected, does the new graph need to be (weakly) connected
#' @param ntry, number of times to try
#' @return A M, the antagonism and mutualism sub-graphs after rewiring
rewire_antago_to_mutual <- function(A, M, need.connected = T, ntry = 100) {
  A[A < 0] = 0 # remove negative elements of antagonism sub-graph
  s = dim(A)[1]  # number of species
  graph <- A + M  # the combined graph
  graph.bin = graph
  graph.bin[graph.bin > 0] = 1  # transfer to binary graph
  if (!is.connected(graph.adjacency(graph.bin, mode = 'directed'), mode = 'weak') && need.connected) {
    stop("The hybrid graph including both mutualism and antagonism interactions MUST be weakly connected!")
  }
  
  if (sum(A != 0) >= 2) {  # if at least two antagonism interactions available
    flag = F # is the rewiring succeed
    for (i in 1:ntry) {
      # pick two random antagonism interactions, then combine them into one mutualism interaction
      indx = which(A > 0, arr.ind = T)  # index of links, e.g. (1,2) (3,4)
      indx.two = indx[sample(1:nrow(indx), 2), ]  # index of two random antagonism interactions
      graph.tmp = graph.bin
      graph.tmp[indx.two[2, 1], indx.two[2, 2]] = 0  # imagine the second chosen antagonism interaction is removed from the hybrid graph
      if (is.connected(graph.adjacency(graph.tmp, mode = 'directed'), mode = 'weak') || !need.connected) {
        M[indx.two[1, 1], indx.two[1, 2]] = A[indx.two[1, 1], indx.two[1, 2]] # move the first chosen antago interaction to mutual
        M[indx.two[1, 2], indx.two[1, 1]] = A[indx.two[2, 1], indx.two[2, 2]]  # move the second chosen antago interaction to mutual
        A[indx.two[1, 1], indx.two[1, 2]] = 0
        A[indx.two[2, 1], indx.two[2, 2]] = 0
        flag = T
        break
      }
    }
    if (flag == F) {
      warning(paste(ntry, 'times has been tried, but the rewiring still donot succeed.'))
    }
    ret = list(A = A, M = M, flag = flag, tried = i)      
  }
  else {
    warning("No enough antagonism interactions available!")
    ret = list(A = A, M = M, flag = F, tried = 0)
  }
  return(ret)
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
        if (p < pc) {  # competitive sub-graph
          competitive_graph[i, j] = 1
          competitive_graph[j, i] = 1
        }
        else if (p < pc + pa / 2) { # antagonism sub-graph
          antago_graph[i, j] = 1
          antago_graph[j, i] = - 1
        }
        else if (p < pc + pa) { # antagonism sub-graph
          antago_graph[i, j] = - 1
          antago_graph[j, i] = 1
        }
        else if ( p < pc + pa + pm) { # mutualism sub-graph
          mutual_graph[i, j] = 1
          mutual_graph[j, i] = 1
        }
      }
    }
  }
  list(competitive_graph = competitive_graph, antago_graph = antago_graph, mutual_graph = mutual_graph)
}

#' @title generate a hybrid network that includes both mutualistic and antagonistic interactions
#' @param s number of species
#' @param k average degree of species
#' @param type network type, 'er':random graph, 'sf':scale-free, 'bipartite':bipartite graph, 'niche':niche model
#' @param pa probability of antagonism interactions
#' @param pm probability of mutualism interactions
#' @param ... additional arguments transformed to graph generate such as [expower]
gen_hybrid_network_2 <- function(s, k, type = 'er', pa = 0., pm = 1., ...) {
  stopifnot(pa >= 0., pm >= 0., pa + pm == 1)
  # ONE mutual interaction is splited to TWO antago interactions
  edges.mutual = floor(s * k * pm)
  #edges.antago = sum(graph > 0) / 2 - edges.mutual
  k = k * (pm + 2 * pa)  # 
  G = gen_connected_graph(s, k, type, ...)  # generate a connected graph
  graph = as.matrix(get.adjacency(G))  # transform to matrix form
  # split the graph to two sub-graphs, antagonism and mutualism graphs, according to the probability of occurance of two different types of interactions
  antago_graph = matrix(0, nrow = s, ncol = s)
  mutual_graph = matrix(0, nrow = s, ncol = s)
  # edges indexed by node pairs, for example (1, 2)
  indx = which(lower.tri(graph) & graph > 0, arr.ind = T )
  # randomly chose mutual edges from edges
  tmp = sample(1:nrow(indx), as.numeric(edges.mutual))
  indx.mutual = indx[tmp, ]
  # the other edges are antago- edges
  if (length(tmp) == 0 | edges.mutual == 0) {
    indx.antago = indx
  }
  else {
    indx.antago = indx[-tmp, ]    
  }
  # the antago- edges should further be partitioned by direction
  tmp = sample(1:nrow(indx.antago), nrow(indx.antago) / 2)
  indx.antago.lower = indx.antago[tmp, ]
  indx.antago.upper = indx.antago[-tmp, ]
  # extract mutual- partition
  mutual_graph = graph
  mutual_graph[upper.tri(mutual_graph)] = 0
  mutual_graph[indx.antago] = 0
  mutual_graph = mutual_graph + t(mutual_graph)
  # generate antago- partition
  antago_graph = graph
  antago_graph[upper.tri(antago_graph)] = 0
  antago_graph[indx.mutual] = 0
  tmp = antago_graph
  tmp[indx.antago.upper] = 0
  antago_graph[indx.antago.lower] = 0
  antago_graph = antago_graph + t(tmp)
  
  list(antago_graph = antago_graph, mutual_graph = mutual_graph)
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