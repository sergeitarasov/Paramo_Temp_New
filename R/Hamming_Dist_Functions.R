
#   ____________________________________________________________________________
#   Functions to calculate Hamming distance from root to tips               ####

# get edges id from tip
get_path_edges <- function(tree.merge, node){
  E <- tree.merge$edge
  E.des <- c()
  
  while (length(node)>0){
    e.id <- which(E[,2]==node)
    E.des <-c(E.des,e.id)
    node <- E[e.id, 1]
  }
  
  return(E.des)
}

# Path <- get_states_path(tree.merge, node=3)
#tree.merge <- tree.list[[1]]
#node=1 # node is a tip
get_states_path <- function(tree.merge, node){
  Maps <- tree.merge$maps
  # get absolute ages
  #H <- phytools:::nodeHeights(tree.merge)
  #age.loc <- lapply(Maps, function(x) cumsum(x))
  #Maps.abs <-mapply(function(x,y) x+y, age.loc, H[,1] )
  
  # get state path
  E.id <- get_path_edges(tree.merge, node)
  edges.foc <- Maps[E.id]
  state.path <- lapply(edges.foc, function(x) rev(x)) %>% unlist %>% rev
  
  # make edge ids
  E.as.id <- mapply(function(MM, EE) rep(EE, length(MM) ), MM=edges.foc, EE=E.id) %>% unlist %>% rev
  #--
  
  Cum <- c(0, state.path)
  Cum <- cumsum(Cum)
  Int <- cbind(Cum[-length(Cum)], Cum[-1]) %>% as_tibble()
  colnames(Int)<-c('t.start', 't.end')
  Int <- mutate(Int,  States=names(state.path) )
  
  Int <- mutate(Int,  Edge.id=E.as.id )
  
  #### Remove duplicated successive states
  #Int <-tibble(t.start=c(1:11), t.end=c(2:12), States=c('0001', '0001', '0001', '0002', '0001', '0003', '0003', '0004', '0004', '0004', '0004') )
  SS <- 1
  while (SS<nrow(Int)) {
    if (Int$States[SS]==Int$States[SS+1]){
      # merge two the same states
      Int$t.end[SS] <- Int$t.end[SS+1] # reassign time
      Int <- Int[-(SS+1),] # remove duplicate
      SS <- SS
    } else
      SS <- SS+1
  }
  ###-------------------
  
  Int <-mutate(Int,  delta.t=t.end-t.start )
  
  return(Int) 
}


#Path <- path_hamming(Path)
# Note some states may be the same if they are separated by split
path_hamming <- function(Path){
  st <- Path$States
  ham.d <- stringdist::stringdist(st[1], st, method = c("hamming") )
  # normalized hamming
  str.len <- nchar(st[1])
  ham.d.norm <-ham.d/str.len
  
  return(mutate(Path, Ham.dist=ham.d, Ham.dist.n=ham.d.norm))
}


#path_hamming_over_tips(tree.merge)
#tree.merge <- tree.list[[1]]
path_hamming_over_tips <- function(tree.merge){
  n.tips <- length(tree.merge$tip.label)
  Dist.tips <- tibble()
  i <- 1
  for (i in 1:n.tips){
    
    P <- get_states_path(tree.merge, node=i)
    Dist <- path_hamming(P)
    
    # add Poissson cumulative change
    Dist <-mutate(Dist, Pois.count=c(1:nrow(Dist)) ) 
    
    Dist.tips <-bind_rows(Dist.tips, mutate(Dist, tip.id=i) )
  }
  return(Dist.tips)
}

#path_hamming_over trees
path_hamming_over_trees <- function(tree.list){
  
  Dist.trees <- tibble()
  i <- 1
  for (i in 1:length(tree.list)){
    tr <- tree.list[[i]]
    tr.i <- path_hamming_over_tips(tr)
    tr.i <-mutate(tr.i, tree.id=i )
    tr.i <-mutate(tr.i, tree.tip.id=paste0(tree.id,'-', tip.id) )
    Dist.trees <-bind_rows(Dist.trees, tr.i )
  }
  return(Dist.trees)
}


# get exact probability of Hamming distanc through time
# Hypercube is a combination of 1d cubes (i.e. two-state Markov process)


# Rate in two state Markov chain
q01 <- function(rate, t) 0.5-0.5*exp(-2*rate*t)


#' Title Probability Distribution function for Hamming distance for binary Hypecube Pr(dist| L, q, t)
#'
#' @param dist distance (Natural N)
#' @param L Dimension of Hypecube
#' @param q rate in 1d Hypercube; q=rate.Hepercube/L
#' @param t time
#'
#' @return Probability of dist
#' @export
#'
#' @examples
#' dist <- seq(0,394, 1)
#' plot(dist, dHamming_dist(dist, L=394, q=3/394, t=7000), type='l' )
#' ## 394*q01(rate=3/394, t=20000) # expectation = L*q01

dHamming_dist <- function(dist, L, q, t){
  # Pr(dist| L, q, t)
  prob <- dbinom(dist, L, q01(rate=q, t) )
  return(prob)
}


# Quantile function
#qHamming_dist(quant=0.025, L=394, q=3/394, t=200)
#qHamming_dist(quant=0.975, L=394, q=3/394, t=200)
qHamming_dist <- function(quant, L, q, t){
  # Pr(dist| L, q, t)
  qq <- qbinom(quant, L, q01(rate=q, t) )
  return(qq)
}


# Log curve
#curve(Log_curve(Range=1, Slope=2, Infl=50, Lowb=0.1, x=t),  0, 200, xname = 't')
Log_curve <- function(Range, Slope, Infl, Lowb, x){
  ((Range-Lowb)/(1+exp(-Slope*(x-Infl))))+Lowb
}
