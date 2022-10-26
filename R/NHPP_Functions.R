# read trees batch
read_STMP_500trees_batch <- function(prefix){
  
  n.trees <- c(1,101,201,301,401)
  
  tree.list <- c()
  for (i in 1:length(n.trees) ){
    cat('Reading from ', prefix, ' from tree ',  n.trees[i], 'to ', n.trees[i]+99, '\n')
    
    tree.list1 <-lapply(c(n.trees[i]:(n.trees[i]+99)), function(x) rds_from_zip(zip.name=paste0(prefix, '.zip'), rds.prefix=paste0(prefix, '_'), rds.id=x)[[1]])
    closeAllConnections()
    tree.list <- c(tree.list, tree.list1)
  }
  
  cat('Merging trees', '\n')
  tree.list <-lapply(tree.list, function(x) merge_tree_cat(x) )
  return(tree.list )
}


# Make data (chnages over one tip) for NHPP stan model
#Tb.trees <- path_hamming_over_trees(tree.list)
#tip_id <- 1
make_data_NHPP_over_tip <- function(Tb.trees, tip_id){
  f.tip <- filter(Tb.trees, tip.id==tip_id)
  change.time <- f.tip$t.start
  change.time <- change.time[-which(change.time==0)]
  return(change.time)
}

make_data_NHPP_over_edge_MarkovKDE <- function(Tb.trees, Focal.Edge){
  #Focal.Edge <- 1
  f.tip <- filter(Tb.trees, Focal.Edge.id==Focal.Edge)
  change.time <- f.tip$t.start
  change.time <- change.time[-which(change.time==0)]
  return(change.time)
}


# Get NHPP data fro stan across all tips
#make_data_NHPP(Tb.trees)
make_data_NHPP <- function(Tb.trees){
  Ntips <- Tb.trees$tip.id %>% unique() %>% length()
  
  dt.out <- list()
  for (i in 1:Ntips){
    dt.out[[i]] <- make_data_NHPP_over_tip(Tb.trees, tip_id=i)
  }
  
  return(dt.out)
}


#Path.data <-make_data_NHPP_KDE_Markov_kernel(Tb.trees)
#Path.data[[1]]
make_data_NHPP_KDE_Markov_kernel <- function(Tb.trees){
  Nedges <- Tb.trees$Focal.Edge.id %>% unique() %>% length()
  
  dt.out <- list()
  i <- 1
  for (i in 1:Nedges){
    #dt.out[[i]] <- make_data_NHPP_over_tip(Tb.trees, tip_id=i)
    dt.out[[i]] <- make_data_NHPP_over_edge_MarkovKDE(Tb.trees, Focal.Edge=i)
  }
  
  return(dt.out)
}


#   ____________________________________________________________________________
#   Diagnostocs for NHPP inference                                          ####

#' Title Dignostocs of Stan output for NHPP across all tips
#'
#' @param stan.post stan output for all tips (length=ntips)
#' @param round.factor factor to round Rhat and then compaer it to 1.00. For example, if Rhat is 1.0001, we round it using factor 2 that gives 1.00 that, in turn, 
#' is equal to 1
#'
#' @return Boolean if all Rhats are equal to 1 given the magnitude of raound facvtor
#' @export
#'
#' @examples
Rhat_diagnostics_boolean <- function(stan.post, round.factor=2){
  
  dt.out <- list()
  for (i in 1:length(stan.post)){
    ss <- summary(stan.post[[i]])$summary
    ss <-tibble(pars=row.names(ss), Rhat=ss[,10])
    dt.out[[i]] <- all(round(ss$Rhat,round.factor)==1)
  }
  return(dt.out)
}

#' Title Return all Rhats for all paramters across all tips
#'
#' @param stan.post  stan output for all tips (length=ntips)
#'
#' @return tibble
#' @export
#'
#' @examples
#' 
Rhat_diagnostics_table <- function(stan.post){
  
  ss <- summary(stan.post[[1]])$summary
  dt.out <- tibble(pars=row.names(ss))
  for (i in 1:length(stan.post)){
    ss <- summary(stan.post[[i]])$summary
    dt.out <-bind_cols(dt.out, Rhat=ss[,10])
  }
  return(dt.out)
}


##  ............................................................................
##  Beta Mixture distr                                                      ####

#Tmax <- 280.8679
#param <- 'mean'
# Beta mixture
dMix_Beta <- function(t, pi, Beta_alpha, Beta_beta, Tmax, param='mean'){
  sum(pi[[param]]*(1/Tmax)*dbeta(t/Tmax, Beta_alpha[[param]], Beta_beta[[param]]))
}

# Beta mixture for vrctors
dMix_Beta_vector <- function(t, pi, Beta_alpha, Beta_beta, Tmax, param){
  sapply(t, function(i) dMix_Beta(i, pi, Beta_alpha, Beta_beta, Tmax, param))
}


##  ............................................................................
##  Discretize Simmap                                                       ####
# for one tree
discr_Simmap<-function(tree, res){
  
  steps <- 0:res/res * max(phytools:::nodeHeights(tree))
  H <- phytools:::nodeHeights(tree)
  maps.n <- vector(mode = "list", length = nrow(tree$edge))
  
  # i=170
  for (i in 1:nrow(tree$edge)) {
    YY <- cbind(c(H[i, 1], steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))]), 
                c(steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))], H[i, 2])) -  H[i, 1]
    
    
    TR<-cumsum(tree$maps[[i]])
    # TR[length(TR)]<-YY[nrow(YY), 2] # this to make the length equal as it sometiems does not hold
    # YY[,1]-YY[,2]
    ######
    #sprintf("%.54f", c(TR[length(TR)], YY[nrow(YY), 2]) )
    # TR[1]==TR[2]
    # all.equal(TR[1], TR[2])
    # all.equal(TR)
    # duplicated(TR)
    #length(int.out)
    #all.equal(TR[length(TR)], YY[nrow(YY), 2])
    #TR[length(TR)]==YY[nrow(YY), 2]
    ######
    
    #TR[length(TR)]==YY[nrow(YY), 2]
    int.out=findInterval(YY[,2], c(0,TR), left.open=T, rightmost.closed = FALSE, all.inside = TRUE)
    #int.out=findInterval(YY[,2], c(0,TR), left.open=T, rightmost.closed = FALSE)
    #findInterval(seq(0.1, 4, .1), c(0, 0.5, 0.7, 1.5, 1.6, 4 ), left.open=T, rightmost.closed = F)
    maps.n[[i]]<-setNames(YY[,2]-YY[,1], names(tree$maps[[i]])[int.out])
  }
  tree$maps<-maps.n
  return(tree)        
}
##
#

# for multiphylo
discr_Simmap_all<-function(tree, res){
  
  if (class(tree)[1]=="simmap") {
    tree<-discr_Simmap(tree, res)
  }
  
  if (class(tree)[1]=="multiSimmap") {
    
    for (j in 1:length(tree)){
      tree[[j]]<-discr_Simmap(tree[[j]], res)
    }
  }
  return(tree)
}
##

# Get Mixture parameters from Stan obj

get_dBetaMix_pars <- function(post.tip){
  Beta_alpha <-summary(post.tip, pars = c("alpha_cl"), probs = c(0.025, 0.975))$summary
  nm <- row.names(Beta_alpha)
  Beta_alpha <-as_tibble(Beta_alpha) %>% bind_cols(pars=nm, .)
  
  Beta_beta <-summary(post.tip, pars = c("beta_cl"), probs = c(0.025, 0.975))$summary
  nm <- row.names(Beta_beta)
  Beta_beta <-as_tibble(Beta_beta) %>% bind_cols(pars=nm, .)
  
  pi <-summary(post.tip, pars = c("pi"), probs = c(0.025, 0.975))$summary
  nm <- row.names(pi)
  pi <-as_tibble(pi) %>% bind_cols(pars=nm, .)
  
  out <- list(Beta_alpha=Beta_alpha, Beta_beta=Beta_beta, pi=pi)
  return(out)
}


#   ____________________________________________________________________________
#   make contMap object                                                     ####

#' Title make contMap object for plotting NHPP
#'
#' @param tree.discr descritized tree to use for mapping
#' @param posterior.NHPP list of stan objects
#' @param Tmax maximum time
#' @param param which paramter to use as an estimate from stan object
#' @param lambda.posterior a vector of lambda estimates for each tip (i.e. density multipliers for NHPP). NULL means 1.
#'
#' @return
#' @export
#'
#' @examples
make_contMap_NHPP <- function(tree.discr, posterior.NHPP, Tmax, param ='mean', lambda.posterior=NULL){
  
  if (is.null(lambda.posterior)) lambda.posterior <- rep(1, Ntip(tree.discr)) 
  
  Maps <- tree.discr$maps
  Maps.mean <- Maps
  Maps.sd<- Maps
  # absolute ages
  H <- phytools:::nodeHeights(tree.discr)
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <-mapply(function(x,y) x+y, age.loc, H[,1] )
  
  # loop over edges
  i <- 1
  for (i in 1:length(age.glob)){
    Edge.focal <- age.glob[[i]]
    # get tips asscoated with this edge
    Node <- tree.discr$edge[i,2]
    Tips <- geiger::tips(tree.discr, Node)
    
    
    # loop over tips descending from this edge and extract probs for each time bin
    Probs <- matrix(ncol=length(Edge.focal), nrow=length(Tips))
    j <- 1
    for (j in 1:length(Tips)){
      tip.id <- which(tree.discr$tip.label==Tips[j])
      post.tip <- posterior.NHPP[[tip.id]]
      Mix.pars <- get_dBetaMix_pars(post.tip)
      
      Probs[j,] <- dMix_Beta_vector(t=Edge.focal, Mix.pars$pi, 
                  Mix.pars$Beta_alpha, Mix.pars$Beta_beta, Tmax=Tmax, param=param)*lambda.posterior[[tip.id]]
    }
    
    # average over tips and get mean ans sd
    Maps.mean[[i]] <- apply(Probs, 2, mean)
    Maps.sd[[i]] <-apply(Probs, 2, sd)
  }
  
  
  #------- assign colors to Maps.mean
  Maps.col <- Maps
  MM <- Maps.mean %>% unlist
  lims <- c(min(MM), max(MM))
  
  cols <- rainbow(1001, start = 0, end = 0.7)
  names(cols) <- 0:1000
  trans <- 0:1000/1000 * (lims[2] - lims[1]) + lims[1]
  names(trans) <- 0:1000
  
  i <- 1
  for (i in 1:length(Maps.col)){
    b <- Maps.mean[[i]]
    d <- sapply(b, phytools:::getState, trans = trans)
    names(Maps.col[[i]]) <- d
  }
  #---------------
  tree.new <- tree.discr
  tree.new$maps <- Maps.col
  tree.new$maps.sd <- Maps.sd
  
  xx <- list(tree = tree.new, cols = cols, lims = lims)
  class(xx) <- "contMap"
  
  return(xx)
}


#   ____________________________________________________________________________
#   make contMap KDE object                                               ####

#' Title make contMap object for plotting NHPP
#'
#' @param tree.discr descritized tree to use for mapping
#' @param posterior.NHPP list of stan objects
#' @param Tmax maximum time
#' @param param which paramter to use as an estimate from stan object
#' @param lambda.posterior a vector of lambda estimates for each tip (i.e. density multipliers for NHPP). NULL means 1.
#'
#' @return
#' @export
#'
#' @examples
#' #NHPP.Map <- make_contMap_NHPP(tree.discr=tree.H.discr, post.test, Tmax=280.8679, param = 'mean', lambda.posterior=lambda.post$Mean)
#' #posterior.NHPP <- post.test
#' 

#Edge.KDE$Maps.mean.norm
#lambda.post.stat='lambda.mean'
make_contMap_KDE <- function(tree.discr, Edge.KDE.stat ){
  
  Maps <- tree.discr$maps
  
  #-- lambda
  #lambda <- lambda.post[[lambda.post.stat]] # !!! PErhaps to include lambda right in Edge.KDE!!!!!!
  #Maps.mean <-lapply(Edge.KDE$Maps.mean.norm, function(x) x*lambda)
  #Maps.mean <-Edge.KDE[[lambda.post.stat]]
  Maps.mean <-Edge.KDE.stat
  
  # absolute ages
  #H <- phytools:::nodeHeights(tree.discr)
  #age.loc <- lapply(Maps, function(x) cumsum(x))
  #age.glob <-mapply(function(x,y) x+y, age.loc, H[,1] )
  
  
  #------- assign colors to Maps.mean
  Maps.col <- Maps
  MM <- Maps.mean %>% unlist
  lims <- c(min(MM), max(MM))
  
  cols <- rainbow(1001, start = 0, end = 0.7) %>% rev
  names(cols) <- 0:1000
  trans <- 0:1000/1000 * (lims[2] - lims[1]) + lims[1]
  names(trans) <- 0:1000
  
  i <- 1
  for (i in 1:length(Maps.col)){
    b <- Maps.mean[[i]]
    d <- sapply(b, phytools:::getState, trans = trans)
    names(Maps.col[[i]]) <- d
  }
  #---------------
  tree.new <- tree.discr
  tree.new$maps <- Maps.col
  #tree.new$maps.sd <- Maps.sd
  
  xx <- list(tree = tree.new, cols = cols, lims = lims)
  class(xx) <- "contMap"
  
  return(xx)
}


make_postPois_KDE <- function(Edge.KDE.stat, lambda.post, lambda.post.stat='Mean'){
  
  #Maps <- tree.discr$maps
  #-- lambda
  lambda <- lambda.post[[lambda.post.stat]] # !!! PErhaps to include lambda right in Edge.KDE!!!!!!
  #XX <-lapply(Edge.KDE$Maps.mean.norm, function(x) x*lambda)
  XX <-lapply(Edge.KDE.stat, function(x) x*lambda)
  return(XX)
}

#   ____________________________________________________________________________
#   Analytical posterioir for lambda NHPP                                   ####

# lambda~Gamma(a+sum(k_i), b+n), when a,b -> 0, then posterior approaches MLn

# setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/EP")
# tree.list <-lapply(c(1:2), function(x) rds_from_zip(zip.name='whole_organism.zip', rds.prefix='whole_organism_', rds.id=x)[[1]])
# tree.list <-lapply(tree.list, function(x) merge_tree_cat(x) )
# closeAllConnections()
# #tree.merge <- tree.list[[1]]
# 
# Tb.trees <- path_hamming_over_trees(tree.list)
# uniq.trees <- Tb.trees$tree.id %>% unique() %>% length()
# X <- filter(Tb.trees, tip.id==1)
# 
# k.i <- sapply(1:uniq.trees, function(y) filter(X, tree.id==y) %>% nrow)
# curve(dgamma(x, 0.01+sum(k.i), 0.01+length(k.i)), 0, 400)
# mean.lambda <- ( 0.01+sum(k.i))/(0.01+length(k.i))

#' Title Analytical posterior of lambda parameter for HNPP. lambda~Gamma(a+sum(k_i), b+n)
#'
#' @param Tb.trees path_hamming_over_trees(tree.list) 
#'
#' @return a list that includes n changes per each path-tip (row) and tree (col), mean, sd and quentiles of lambda for each path-tip
#' @export
#'
#' @examples
posterior_lambda <- function(Tb.trees){
  uniq.trees <- Tb.trees$tree.id %>% unique() %>% length()
  uniq.tips <- Tb.trees$tip.id %>% unique() %>% length()
  
  out <- matrix(ncol = uniq.trees, nrow=uniq.tips)
  # loop over tips
  i <- 1
  for (i in 1:uniq.tips){
    X <- filter(Tb.trees, tip.id==i)
    k.i <- sapply(1:uniq.trees, function(y) filter(X, tree.id==y) %>% nrow)
    out[i,] <- k.i-1
  }
  
  # get rquired statistics
    Mean <- apply(out, 1, function(x) ( 0.01+sum(x))/(0.01+length(x)) )
    Q_2.5 <- apply(out, 1, function(x) qgamma(p=0.025, 0.01+sum(x), 0.01+length(x)) )
    Q_97.5 <- apply(out, 1, function(x) qgamma(p=0.975, 0.01+sum(x), 0.01+length(x)) )
    SD <- apply(out, 1, function(x) sqrt(( 0.01+sum(x))/(0.01+length(x))^2)  )
    
    XX <- list(n.chages.path=out, Mean=Mean, SD=SD, Q_2.5=Q_2.5, Q_97.5=Q_97.5)
    return(XX)

}

#------- KDE




posterior_lambda_KDE <- function(tree.list){
  # uniq.trees <- Tb.trees$tree.id %>% unique() %>% length()
  # uniq.tips <- Tb.trees$tip.id %>% unique() %>% length()
  # 
  # # n changes
  # n.changes <- c()
  # for (i in 1:uniq.trees){
  #   n.changes[i] <-filter(Tb.trees, tree.id==i) %>% nrow-1
  # }
  
  #class(tree.list) <- append(class(tree.list),"multiSimmap")
  n.changes <-lapply(tree.list, function(x) countSimmap(x)$N) %>% unlist
  
  
  # get rquired statistics
  Mean <-  (0.01+sum(n.changes))/(0.01+length(n.changes)) 
  SD <- sqrt(( 0.01+sum(n.changes))/(0.01+length(n.changes))^2)
  Q_2.5 <- qgamma(p=0.025, 0.01+sum(n.changes), 0.01+length(n.changes))
  Q_97.5 <- qgamma(p=0.975, 0.01+sum(n.changes), 0.01+length(n.changes)) 
  
  XX <- list(Mean=Mean, SD=SD, Q_2.5=Q_2.5, Q_97.5=Q_97.5)
  return(XX)
  
}


#tree.list <- Trees.focal.BRs[[1]]
#BR.name <- 'br'
posterior_lambda_KDE_Distr <- function(tree.list, n.sim=10, BR.name){
  # uniq.trees <- Tb.trees$tree.id %>% unique() %>% length()
  # uniq.tips <- Tb.trees$tip.id %>% unique() %>% length()
  # 
  # # n changes
  # n.changes <- c()
  # for (i in 1:uniq.trees){
  #   n.changes[i] <-filter(Tb.trees, tree.id==i) %>% nrow-1
  # }
  
  #class(tree.list) <- append(class(tree.list),"multiSimmap")
  n.changes <-lapply(tree.list, function(x) countSimmap(x)$N) %>% unlist
  
  #----
  sim <- rgamma(n.sim, 0.01+sum(n.changes), 0.01+length(n.changes))
  out <- tibble(BR=BR.name, sim)

  return(out)
  
}

##  ............................................................................
##  NHPP Finite Mixture Stan Model                                         ####

NHPP_finite_Mix <-"
data{
int<lower=0> C; //num of cludter
int<lower=0> N; //data num
real<lower=0> y[N]; //count data
real<lower=0> Tmax; //time max
}

parameters {
simplex [C] pi;
//vector<lower=0, upper=Tmax> [C] mu_cl; //cluster mu
//vector<lower=0> [C] tau_cl; //cluster tau
positive_ordered [C] mu_cl; //cluster mu
vector<lower=0> [C] tau_cl; //cluster tau
}

transformed parameters{

vector<lower=0>[C] beta_cl;
vector<lower=0>[C] alpha_cl;
alpha_cl=(mu_cl .* tau_cl ) ./ Tmax;
beta_cl=(1 - (mu_cl ./ Tmax) ) .* tau_cl;

}

model {

real ps[C];
pi ~ dirichlet(rep_vector(1.0, C));
mu_cl ~ uniform(0, Tmax);
tau_cl ~ inv_gamma(2,8);

for(i in 1:N){
for(c in 1:C){
ps[c]=log(pi[c]) + log(1/Tmax) + beta_lpdf((y[i]/Tmax) | alpha_cl[c], beta_cl[c]);
}
target += log_sum_exp(ps);
}

}
" 
#--------





#   ____________________________________________________________________________
#   NHPP KDE                                                                ####

KDE_unnormalized_scalar <- function(x, h, dat){
  sum(dnorm(x, dat, h))
}

KDE_unnorm <- Vectorize(KDE_unnormalized_scalar, vectorize.args='x')

#my_den_phy(c(3), h=1, dat = c(e1,e2,e3))
#v_my_den_phy(c(1,1,2,3), h=1, dat = c(e1,e2,e3))

# truncated normal
library("truncnorm")
KDE_unnormalized_scalar_trunc <- function(x, h, dat, start, end){
  dtruncnorm(x, a=start, b=end, mean = dat, sd = h) %>% sum()
  #sum(dnorm(x, dat, h))
}

KDE_unnorm_trunc <- Vectorize(KDE_unnormalized_scalar_trunc, vectorize.args='x')


# truncated normal experimental
#curve(dtruncnorm(x, a=-10, b=10, mean = 0, sd = 5), -20, 20)

library("truncnorm")
KDE_unnormalized_scalar_trunc_exp <- function(x, h, dat, start, end){
  dtruncnorm(x, a=start, b=end, mean = dat, sd = h) %>% sum()
  #sum(dnorm(x, dat, h))
}

KDE_unnorm_trunc <- Vectorize(KDE_unnormalized_scalar_trunc_exp, vectorize.args='x')



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Calculate KDE over edges                                                ####

#' Title Estimates unnormilized KDE over each edge of a discretized phylogeny. KDE for each edge is averages across all possible root-tip paths
#'
#' @param tree.discr descritezed tree
#' @param root.taxon path root-tip to use as a rooting
#' @param Path.data path.data object
#'
#' @return mean KDE for edge (Maps.mean), and mean absolute error (Maps.error)
#' @export
#'
#' @examples
#' 
#' estimate_edge_KDE_unnorm(tree.discr, root.taxon=87, Path.data)
#' 
estimate_edge_KDE_unnorm <- function(tree.discr, root.taxon=87, Path.data, band.width='bw.ucv'){
  
  Maps <- tree.discr$maps
  #Maps.mean <- Maps
  #Maps.error<- Maps
  # absolute ages
  H <- phytools:::nodeHeights(tree.discr)
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <-mapply(function(x,y) x+y, age.loc, H[,1] )
  #--
  #h <- 10 # bandwidth
  
  # Root taxon for KDE
  # in our case it's otu 87 that is edge 172
  #root.taxon <- 87
  root.ed <- get_path_edges(tree.discr, root.taxon)
  root.changes <- -1*Path.data[[root.taxon]]
  #hist(root.changes, breaks = 50)
  #root.time <- -1*age.glob[[root.ed]]
  
  
  #tree <- tree.list[[1]]
  # tree <- read.tree(
  #   '/Users/taravser/Documents/My_papers/onto_phylo/data/working/dataset/Data_final/Hymenoptera_br_resolved.tre')
  
  # loop over edges
  
  Maps.mean <- list()
  Maps.error <- list()
  h.list <- list()
  
  Edge2search <- c(1:length(Maps))
  Edge2search <-Edge2search[-which(Edge2search==root.ed)]
  Band.width <- vector(length = length(Maps), mode='list')
  
  #Edge <- 1 
  #for (Edge in 1:length(Maps)){
  for (Edge in Edge2search){
    
    cat('Working on edge ', Edge, ' \n')
    # get node corresponding to edge
    Node <- tree.discr$edge[Edge,2]
    # start and end times of an edge
    #edge.times <- phytools:::nodeHeights(tree)[Edge,]
    
    Tips.descen <- phangorn::Descendants(tree.discr, node=Node, type = c("tips"))[[1]]
    
    #--- Get density unnormalized over edge interval
    Dens.tips <- vector(length = length(Tips.descen), mode='list')
    #Win <- seq(edge.times[1], edge.times[2], len=100)
    Win <- age.glob[[Edge]]
    i <- 1
    for(i in 1:length(Tips.descen)){
      dt <- Path.data[[Tips.descen[i]]]
      # add rot to dt
      dt <- c(root.changes, dt)
      #hist(dt, breaks = 50)
      
      if (band.width=='bw.nrd')
        h <- bw.nrd0(dt)
      if (band.width=='bw.ucv')
        h <-bw.ucv(dt, lower = 0.01, upper = 20)
      if (band.width=='bw.bcv')
        h <-bw.bcv(dt, lower = 0.01, upper = 20)
      if (band.width=='bw.SJ')
        h <-bw.SJ(dt, lower = 0.01, upper = 20)
      
      Band.width[[Edge]] <- c(Band.width[[Edge]], h)
      
      Dens.tips[[i]] <- KDE_unnorm(Win, h=h, dat=dt) #/length(dt)
      #h.list[[i]] <-h
      # unlist(h.list) %>% hist
    }
    
    # sum denisties over interval and nomalize them
    # normalizing constants per each tip-path
    #Norm.const <- lapply(Path.data[Tips.descen], length) %>% unlist
    DD <- Reduce('+', Dens.tips)/length(Dens.tips) #/sum(Norm.const)
    
    #Error <- lapply(Dens.tips, function(x) (DD-x)^2 )
    Error <- lapply(Dens.tips, function(x) abs(DD-x) )
    Error <-Reduce('+', Error)/length(Dens.tips) # mean sq error
    
    Maps.mean[[Edge]] <- DD
    Maps.error[[Edge]] <- Error
    #Maps.error[[Edge]]
    
  }
  
  # getting KDE for root edge
  # Edge2search <- c(1:length(Maps))
  # Edge2search <-Edge2search[-which(Edge2search==root.ed)]
  
  Tips <- (1:Ntip(tree.discr))
  Tips.descen <-Tips[-which(Tips==root.taxon)]
  #Tips.descen <- Descendants(tree.discr, node=Node, type = c("tips"))[[1]]
  
  Dens.tips <- vector(length = length(Tips.descen), mode='list')
  #Win <- seq(edge.times[1], edge.times[2], len=100)
  Win <- age.glob[[root.ed]]
  
  cat('Working on edge ', root.ed, ' \n')
  i <- 1
  for(i in 1:length(Tips.descen)){
    dt <- Path.data[[Tips.descen[i]]]
    # add rot to dt
    dt <- -1*c(root.changes, dt)
    #hist(dt, breaks = 50)
    
    if (band.width=='bw.nrd')
      h <- bw.nrd0(dt)
    if (band.width=='bw.ucv')
      h <-bw.ucv(dt, lower = 0.01, upper = 20)
    if (band.width=='bw.bcv')
      h <-bw.bcv(dt, lower = 0.01, upper = 20)
    if (band.width=='bw.SJ')
      h <-bw.SJ(dt, lower = 0.01, upper = 20)
    
    Band.width[[root.ed]] <- c(Band.width[[root.ed]], h)
    
    Dens.tips[[i]] <- KDE_unnorm(Win, h=h, dat=dt) #/length(dt)
  }
  
  DD <- Reduce('+', Dens.tips)/length(Dens.tips) #/sum(Norm.const)
  
  #Error <- lapply(Dens.tips, function(x) (DD-x)^2 )
  Error <- lapply(Dens.tips, function(x) abs(DD-x) )
  Error <-Reduce('+', Error)/length(Dens.tips) # mean sq error
  
  Maps.mean[[root.ed]] <- DD
  Maps.error[[root.ed]] <- Error
  
  return(list(Maps.mean=Maps.mean, Maps.error=Maps.error, Band.width=Band.width))

}



##  ............................................................................
##  Markov Kernel                                                           ####


library("truncnorm")
KDE_unnormalized_scalar_Markov_kernel <- function(x, h, dat){
  dtruncnorm(x, a=dat, b=Inf, mean = dat, sd = h) %>% sum()
  #sum(dnorm(x, dat, h))
}

KDE_unnorm_trunc_Markov <- Vectorize(KDE_unnormalized_scalar_Markov_kernel, vectorize.args='x')


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Estimate bandwidth                                                      ####
# Uses path root-tips for estimation
#data.path <- Pseudo.path.data
#length(data.path)

Estimate_band_W <- function(tree.discr, data.path, band.width=c('bw.nrd0', 'bw.nrd0', 'bw.ucv', 'bw.bcv', 'bw.SJ') ){
  
  otus <- c(1:length(tree.discr$tip.label))
  #tree.discr$edge[,2] %in% otus
  Edges <- match(otus, tree.discr$edge[,2])
  
  Band.width <- c()
  i <- 14
  for (i in Edges){
    cat('Working on edge ', i, ' \n')
    dt <- data.path[[i]]
    
    if (band.width=='bw.nrd')
      h <- bw.nrd(dt)
    if (band.width=='bw.nrd0')
      h <- bw.nrd0(dt)
    if (band.width=='bw.ucv')
      h <-bw.ucv(dt, lower = 0.01, upper = 20)
    if (band.width=='bw.bcv')
      h <-bw.bcv(dt, lower = 0.01, upper = 20)
    if (band.width=='bw.SJ')
      h <-bw.SJ(dt, lower = 0.01, upper = 20)
    
    Band.width <- c(Band.width, h)
  }

  return(Band.width)
}


#   ____________________________________________________________________________
#   Loess Smoothing                                                         ####

library("fANCOVA")

#Edges <-c(1:length(Edge.KDE[[2]]$Maps.mean) )
Loess_smooting_KDE <- function(tree.discr, Edge.KDE ){
  
  # plot edge profiles
  Maps <- tree.discr$maps
  # absolute ages
  H <- phytools:::nodeHeights(tree.discr)
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <-mapply(function(x,y) x+y, age.loc, H[,1] )
  
  Edges <-c(1:length(Edge.KDE$Maps.mean) )
  

  Maps.mean.loess <- vector(length = length(Maps), mode = 'list' )
  i <- 106
  for (i in Edges){
    #for (i in 105:105){
    cat('Working on edge ', i, ' \n')
    dt <- tibble(X=age.glob[[i]], Y=Edge.KDE$Maps.mean[[i]] )
    #dt <- tibble(X=age.glob[[i]], Y=Edge.KDE[[2]]$Maps.mean[[i]] )
    
    #loessMod <- loess(Y ~ X, data=dt, span=span)
    #loessMod <- loess.as(dt$X, dt$Y, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F)
    loessMod <- loess.as(dt$X, dt$Y, degree = 1, 
      criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F, 
      control = loess.control(surface = "direct") )
    
    smoothed <- predict(loessMod)
    #str(smoothed)
    Maps.mean.loess[[i]] <- smoothed
  }
  
  return(Maps.mean.loess)
}


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Add Pseudodata                                                          ####

# ADD PSEUDODATA

add_pseudodata <- function(Edge.goups, Pseudo.data, Path.data){
  
  Pseudo.path.data <- vector(length = length(Path.data), mode = 'list')

  for (i in 1:length(Edge.goups)){
    j <- 1
    for(j in 1:length(Edge.goups[[i]])){
      Pseudo.path.data[[Edge.goups[[i]][j]]] <- c(Pseudo.data[[i]], Path.data[[Edge.goups[[i]][j]]] )
    }
  }

return(Pseudo.path.data)
}
  



#' Title Estimates unnormilized KDE over each edge of a discretized phylogeny. KDE for each edge is averages across all possible root-tip paths
#'
#' @param tree.discr descritezed tree
#' @param root.taxon path root-tip to use as a rooting
#' @param Path.data path.data object
#'
#' @return mean KDE for edge (Maps.mean), and mean absolute error (Maps.error)
#' @export
#'
#' @examples
#' 
#' estimate_edge_KDE_unnorm(tree.discr, root.taxon=87, Path.data)
#' 
#Edge.KDE <-estimate_edge_KDE_Markov_kernel_unnorm(tree.discr, Path.data=Pseudo.path.data , h=15)
#Path.data <- Pseudo.path.data
estimate_edge_KDE_Markov_kernel_unnorm <- function(tree.discr, Path.data, h=10){
  #root.taxon=87
  Maps <- tree.discr$maps
  #Maps.mean <- Maps
  #Maps.error<- Maps
  # absolute ages
  #tree.discr$edge
  H <- phytools:::nodeHeights(tree.discr)
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <-mapply(function(x,y) x+y, age.loc, H[,1] )
  

  
  Density <- vector(length = length(Maps), mode='list')
  Edge <- 1
  for (Edge in 1:nrow(H)){
    cat('Working on edge ', Edge, ' \n')
    X <- age.glob[[Edge]]
    Y <- Path.data[[Edge]]
    #-----------
    #hist(Y, breaks=50)
    
    #Y <-Y[Y>0] !!! UNMASK
    
    #--
    #if (length(Y)==0){Y <- 0}
    #Y <- 0
    #density(Y)
    #unique(Y)
    #------------
    if (length(Y)>0){
    Density[[Edge]] <- KDE_unnorm_trunc_Markov(x=X, h=h, dat=Y)
    }
    
    if (length(Y)==0){
      Density[[Edge]] <-rep(0, length(X) )  
      }
  }
  
  return(list(Maps.mean=Density))
}
  

  

  
  
  
### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Trash                                                                   ####

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Normilized KDE over edges                                              ####

estimate_edge_KDE <- function(tree.discr, Path.data, h){
  
  #Edge.KDE <- estimate_edge_KDE_unnorm(tree.discr, root.taxon=87, Path.data, band.width)
  Edge.KDE <-estimate_edge_KDE_Markov_kernel_unnorm(tree.discr, Path.data=Path.data , h=h)
  #-- normilize KDE
  #res <- 1000
  #Tmax <- max(phytools:::nodeHeights(tree.discr))
  #ss <- 0:res/res * max(phytools:::nodeHeights(tree.discr))
  #ss[1]-ss[2]
  #delta <- Tmax/1000
  
  Maps.dx <- tree.discr$maps
  Maps.mean <- Edge.KDE$Maps.mean # dy
  
  # estimations
  dy.dx<-mapply(function(x,y) y*x, y=Maps.mean, x=Maps.dx)
  sum.dy.dx <- lapply(dy.dx, sum)
  total.sum <- unlist(sum.dy.dx) %>% sum()
  
  Maps.mean.norm<- lapply(Maps.mean, function(x) x/total.sum)
  
  Edge.KDE$Maps.mean.norm <- Maps.mean.norm
  
  return(Edge.KDE)
}
  

normilize_KDE <- function(tree.discr, Maps.mean.loess){
  
  Maps.dx <- tree.discr$maps
  Maps.mean <-  Maps.mean.loess # dy
  
  # estimations
  dy.dx<-mapply(function(x,y) y*x, y=Maps.mean, x=Maps.dx)
  sum.dy.dx <- lapply(dy.dx, sum)
  total.sum <- unlist(sum.dy.dx) %>% sum()
  
  Maps.mean.norm<- lapply(Maps.mean, function(x) x/total.sum)
  
  #Edge.KDE$Maps.mean.norm <- Maps.mean.norm
  
  return(Maps.mean.norm)
}


#Edge.KDE.stat <- Edge.KDE$loess.lambda.mean
derivative_KDE <- function(tree.discr, Edge.KDE.stat){
  
  # dx <- tree.discr$maps
  # y <-  Edge.KDE.stat # y
  # 
  # # Claculate dy
  # #Maps.mean <- list(c(1,2,3,5), c(8,5,6))
  # dy <- lapply(Maps.mean, function(x) c(x[-1]-x[-length(x)], x[length(x)]-x[(length(x)-1)]) )
  # 
  # # estimations
  # dy.dx<-mapply(function(x,y) y/x, y=dy, x=Maps.dx)
  
  #------------
  dx.list <- tree.discr$maps
  y.list <-  Edge.KDE.stat # y
  E.map<-tree.discr$edge
  
  Map.deriv <- vector(length = nrow(E.map), mode='list')
  
  # getting non root edge
  roots <- c(1:nrow(E.map))[!(E.map[,1] %in% E.map[,2])]
  not.roots <- c(1:nrow(E.map))[(E.map[,1] %in% E.map[,2])]
  
  # Working on non roots
  i <- 7
  for (i in not.roots){
    y <- y.list[[i]]
    dx <- dx.list[[i]]
    
    # add first element from combination ancestor + current edge
    ances.edge <- which(E.map[,2] == E.map[i,1])
    #y.anc <- y.list[[ances.edge]][length(y.list[[ances.edge]])]
    #dx.anc <- dx.list[[ances.edge]][length(dx.list[[ances.edge]])]
    y.anc <- y.list[[ances.edge]]
    dy.last.anc <- y.anc[[length(y.anc)]]-y.anc[[length(y.anc)-1]]
    dx.last.anc <- dx.list[[ances.edge]][length(dx.list[[ances.edge]])]
    Der.last.anc <- dy.last.anc/dx.last.anc
    
    #--- average between acestor and fisrt derivatives
    dy <- y[-1]-y[-length(y)]
    dx <- dx[-1]
    Der <- dy/dx
    Der <- c((Der[1]+Der.last.anc)/2, Der)
    #dy.first <- y[1]-y.anc
    #dy.first <- (dy.last.anc+dy[[1]])/2
    #dy <- c(dy.first, dy)
    
    #dx.first <- dx[1]+dx.anc
    #dx[[1]] <- dx.first
    #----
    
    Map.deriv[[i]] <- Der
    #Map.deriv[[i]] <- dy/dx
    
  }
  
  # Working on roots, just duplicating n+1 derivative
  i <- 1
  for (i in roots){
    y <- y.list[[i]]
    dx <- dx.list[[i]]
    
    dy <- y[-1]-y[-length(y)]
    dy <- c(dy[[1]], dy)
    
    Map.deriv[[i]] <- dy/dx
    
  }
  
  return(Map.deriv)
}
  
### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### calculate KDE integral over edges                                              ####

integrate_edge_KDE <- function(tree.discr, Edge.KDE.list){
  
  #Edge.KDE <- estimate_edge_KDE_unnorm(tree.discr, root.taxon=87, Path.data)
  
  #-- normilize KDE
  #res <- 1000
  #Tmax <- max(phytools:::nodeHeights(tree.discr))
  #ss <- 0:res/res * max(phytools:::nodeHeights(tree.discr))
  #ss[1]-ss[2]
  #delta <- Tmax/1000
  
  Maps.dx <- tree.discr$maps
  #Maps.mean <- Edge.KDE$Maps.mean # dy
  Maps.mean <- Edge.KDE.list #dy
  
  # estimations
  dy.dx<-mapply(function(x,y) y*x, y=Maps.mean, x=Maps.dx)
  sum.dy.dx <- lapply(dy.dx, sum)
  total.sum <- unlist(sum.dy.dx) %>% sum()
  
  names(total.sum) <- 'Sum of all edges'
  
  return(total.sum)
}


integrate_path <- function(delta, Y){
  sum(Y*delta)
}


# Make edge profiles for plotting 

edge_profiles4plotting <- function(tree.discr, Edge.KDE.stat){
  # plot edge profiles
  Maps <- tree.discr$maps
  # absolute ages
  H <- phytools:::nodeHeights(tree.discr)
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <-mapply(function(x,y) x+y, age.loc, H[,1] )
  
  # Something is wrong here !!!!!
  v.x <- age.glob %>% enframe(name = "edge.id", value = "X") %>% unnest()
  #v.x$edge.id %>% unique()
  #v.y <-Edge.KDE$lambda.mean  %>% enframe(name = "edge.id", value = "Y") %>% unnest()
  v.y <-Edge.KDE.stat  %>% enframe(name = "edge.id", value = "Y") %>% unnest()
  
  edge.profs <- bind_cols(v.x,v.y)
  #edge.profs <-select(edge.profs, -edge.id1) %>% mutate(edge.type='main')
  edge.profs <- mutate(edge.profs, edge.type='main')
  
  #all(edge.profs$edge.id==edge.profs$edge.id1)
  # add joins for edges
  edge.profs.joint <-join_edges(tree.discr, edge.profs)
  edge.profs.joint <-mutate(edge.profs.joint,edge.type='joint')
  
  edge.profs <- bind_rows(edge.profs, edge.profs.joint)
  
  #-------------
  # v.x <- age.glob %>% enframe(name = "edge.id", value = "X") %>% unnest()
  # #v.x$edge.id %>% unique()
  # v.y <-Edge.KDE$lambda.mean %>% enframe(name = "edge.id", value = "Y") %>% unnest()
  # 
  # edge.profs <- bind_cols(v.x,v.y)
  # 
  # # add joins for edges
  # edge.profs.joint <-join_edges(tree.discr, edge.profs)
  # 
  # edge.profs <- bind_rows(edge.profs, edge.profs.joint)
  
  return(edge.profs)
}

#- join neighboring edges in edge profiles
join_edges<-function(tree.discr, edge.profs){
  
  E.map<-tree.discr$edge
  
  # getting non tip edges
  tips <- c(1:Ntip(tree.discr))
  E.map.nontip <- c(1:nrow(E.map))[!(E.map[,2] %in% tips)]
  
  #Non.tip.edge <- edge.profs$edge.id %>% unique()
  
  Joints<-tibble()
  i <- 14
  for (i in E.map.nontip ){
    
    # focal edge, last X,Y coords
    parent.edge <- filter(edge.profs, edge.id==i)
    Par <- parent.edge[nrow(parent.edge),]
    
    # descnedant eges
    descen.edges <- which(E.map[,1] == E.map[i,2])
    #descen.edges <- E.map[(E.map[,1] %in% E.map[i,2]),]
    Des1 <- filter(edge.profs, edge.id==descen.edges[1])[1,]
    Des2 <- filter(edge.profs, edge.id==descen.edges[2])[1,]
    
    # joints<- bind_rows(Par, Des1, Par, Des2)
    # joints<-mutate(joints, groups=paste0(i, c('a', 'a', 'b', 'b')))
    joints<- bind_rows(Par, Des1, Par, Des2)
    #joints<-mutate(joints, edge.id=paste0(i, c('a', 'a', 'b', 'b')))
    joints<-mutate(joints, edge.id=(i + c(0.1,0.1, .2,.2)) )
    
    Joints<- bind_rows(Joints, joints)
  }
  
  return(Joints)
}

# Helping functions

Deriv_one_edge <- function(X, Y){
  delta.x <- X[-1]-X[-length(X)]
  delta.y <-Y[-1]-Y[-length(X)]
  D <- delta.y/delta.x
  Tib <- tibble(X=X[-length(X)], Y=D)
  return(Tib)
}


Deriv_edges <- function(edge.profs, power=1){
  uniq.edges <- edge.profs$edge.id1 %>% unique() %>% length()
  
  Tib <- tibble()
  i <- 6
  for (i in 1:uniq.edges){
    X <- filter(edge.profs, edge.id==i)$X
    Y <- filter(edge.profs, edge.id==i)$Y
    #plot(X, Y, type='l')
    Dd <- Deriv_one_edge(X,Y^power)
    Dd <- mutate(Dd, edge.id=rep(i, nrow(Dd)) )
    Tib <-bind_rows(Tib, Dd)
    
  }
  
  Tib <-bind_rows(Tib, join_edges(tree.discr, Tib) )
  return(Tib)
}

Deriv_path <- function(delta, Y){
  Y1 <-Y[-1]
  Y0 <- Y[-length(Y)]
  delta.y <- Y1-Y0
  D <- delta.y/delta
  return(D)
}


#------ Edge Paths

edge_paths4plotting <- function(tree.discr, Tb.trees, delta=0.1, root.taxon=87, band.width='bw.ucv'){
  # ??? root separately
  
  # Plot all paths
  lambda.paths <- posterior_lambda(Tb.trees)
  
  # normilize paths
  Tmax <- max(phytools:::nodeHeights(tree.discr))
  #delta <- 0.1
  X <- seq(0, Tmax, delta)
  
  #Path.tibble <- tibble(Edge.n=NA, Y.norm=NA)
  #root.taxon <- 87
  root.ed <- get_path_edges(tree.discr, root.taxon)
  root.changes <- -1*Path.data[[root.taxon]]
  
  Path.tibble <- tibble()
  
  #Edge.n <- 87
  for (Edge.n in 1:(length(Path.data))){ #-1 means except root
    
    cat('Working on path ', Edge.n, ' \n')
    
    dt <- c(root.changes, Path.data[[Edge.n]] )
    #hist(dt)
    if (band.width=='bw.nrd')
      h <- bw.nrd0(dt)
    if (band.width=='bw.ucv')
      h <-bw.ucv(dt, lower = 0.01, upper = 20)
    if (band.width=='bw.bcv')
      h <-bw.bcv(dt, lower = 0.01, upper = 20)
    if (band.width=='bw.SJ')
      h <-bw.SJ(dt, lower = 0.01, upper = 20)
    #h <- bw.nrd0(dt)
    #h <-bw.bcv(dt, lower = 0.01, upper = 20)
    Y <-KDE_unnorm(X, h=h, dt)
    Const <- integrate_path(delta=delta, Y=Y)
    #integrate_path(delta=delta, Y=Y/Const)
    Y.norm <- Y/Const
    
    # get lambda using analytical lambda
    Y.lambda <- Y.norm*lambda.paths$Mean[Edge.n]
    
    Path.tibble <- bind_rows(Path.tibble, as_tibble(cbind(Edge.n, Y.lambda, X)) )
    
  }
  
  
  return(Path.tibble)
  
}


# Derivative
edge_paths_derivative <- function(Path.tibble, delta=0.1, power=1){
  uniq.tips <- Path.tibble$Edge.n %>% unique() %>% length()
  
  Path.tibble.deriv <- tibble()
  # loop over tips
  i <- 1
  for (i in 1:uniq.tips){
    Y.fil <- filter(Path.tibble, Edge.n==i)
    Y.D <- Deriv_path(delta = delta, Y=Y.fil$Y.lambda^power)
    X <- Y.fil$X[-length(Y.fil$X)]
    Path.tibble.deriv <- bind_rows(Path.tibble.deriv, cbind(Edge.n=i, Y.lambda.Deriv= Y.D, X=X) %>% as_tibble() )
  }
  
  return(Path.tibble.deriv)
}


# Filter points for plots
filter_points4plot <- function(edge.profs, edge.profs.derivPOW, round.factor=6){
  
  filtered.vals <- round(edge.profs.derivPOW$Y, round.factor)
  filtered.vals <-  which(filtered.vals!=0)
  Filtered.derivative <- edge.profs.derivPOW[filtered.vals,]
  
  # filte in edge profiles
  uniq.tips <- edge.profs$edge.id1 %>% unique() %>% length()
  Tib <- tibble()
  # loop over tips
  i <- 1
  for (i in 1:uniq.tips){
    Y.fil <-filter(edge.profs, edge.id==i) 
    Y.fil <- Y.fil[-nrow(Y.fil),]
    Tib <-bind_rows(Tib, Y.fil)
  }
  Filtered.function <-Tib[filtered.vals,]
  
  return(list(Filtered.derivative=Filtered.derivative, Filtered.function=Filtered.function))
}


filter_points4plot_withJoints <- function(edge.profs, edge.profs.derivPOW, round.factor=6){
  
  filtered.vals <- round(edge.profs.derivPOW$Y, round.factor)
  filtered.vals <-  which(filtered.vals!=0)
  Filtered.derivative <- edge.profs.derivPOW[filtered.vals,]
  
  # filte in edge profiles
  uniq.tips <- edge.profs$edge.id1 %>% unique() %>% length()
  Tib <- tibble()
  # loop over tips
  i <- 1
  for (i in 1:uniq.tips){
    Y.fil <-filter(edge.profs, edge.id==i) 
    Y.fil <- Y.fil[-nrow(Y.fil),]
    Tib <-bind_rows(Tib, Y.fil)
  }
  
  Filtered.function <-Tib[filtered.vals,]
  
  return(list(Filtered.derivative=Filtered.derivative, Filtered.function=Filtered.function))
}
#max(filtered.vals)




#   ____________________________________________________________________________
#   Path hamming for NHPP                                                   ####
#get_path_edges(tree.merge, node=90)
# Path <- get_states_path(tree.merge, node=3)
#tree.merge <- tree.list[[1]]
#node=90 # node is a tip
#get_states_path_v2 is based on get_states_path
get_states_path_v2 <- function(tree.merge, node){
  Maps <- tree.merge$maps
  # get absolute ages
  #H <- phytools:::nodeHeights(tree.merge)
  #age.loc <- lapply(Maps, function(x) cumsum(x))
  #Maps.abs <-mapply(function(x,y) x+y, age.loc, H[,1] )
  
  # get state path
  tree.merge$edge
  #get_path_edges(tree.merge, node=95)
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
  # Int <-tibble(t.start=c(1:11), t.end=c(2:12), 
  #              States=c('0001', '0001', '0001', '0002', '0001', '0003', '0003', '0004', '0004', '0004', '0004'),
  #              Edge.id=c(1,1,1,1, 2,2, 3, 3,3,3,4))
  # 
  SS <- 1
  while (SS<nrow(Int)) {
    if (Int$States[SS]==Int$States[SS+1]){
      # merge two the same states
      Int$t.end[SS] <- Int$t.end[SS+1] # reassign time
      Int$Edge.id[SS] <- Int$Edge.id[SS+1] # reassign edge.id; edge.id shows the id of edge when change happens
      Int <- Int[-(SS+1),] # remove duplicate
      SS <- SS
    } else
      SS <- SS+1
  }
  ###-------------------
  
  Int <-mutate(Int,  delta.t=t.end-t.start )
  
  return(Int) 
}

# get_states_path(tree.merge, node=88)
# get_states_path(tree.merge, node=95)
#path_hamming_over_tips(tree.merge)

path_hamming_over_all_edges <- function(tree.merge){
  #n.tips <- length(tree.merge$tip.label)
  #Nodes <- tree.merge$edge[,2] %>% unique()
  #n.tips <- nrow(tree.merge$edge)
  Dist.tips <- tibble()
  i <- 2
  #for (i in 1:n.tips){
  for (i in 1:nrow(tree.merge$edge)){
    
    Node <- tree.merge$edge[i,2]
    P <- get_states_path_v2(tree.merge, node=Node)
    Dist <- path_hamming(P)
    
    # add Poissson cumulative change
    Dist <-mutate(Dist, Pois.count=c(1:nrow(Dist)) ) 
    
    Dist.tips <-bind_rows(Dist.tips, mutate(Dist, Focal.Edge.id=i) )
  }
  return(Dist.tips)
}

#Tb.trees <- path_hamming_over_trees_KDE(tree.list)
path_hamming_over_trees_KDE <- function(tree.list){
  
  Dist.trees <- tibble()
  i <- 1
  for (i in 1:length(tree.list)){
    tr <- tree.list[[i]]
    #tr.i <- path_hamming_over_tips(tr)
    tr.i <-path_hamming_over_all_edges(tr)
    tr.i <-mutate(tr.i, tree.id=i )
    tr.i <-mutate(tr.i, tree.tip.id=paste0(tree.id,'-', Focal.Edge.id) )
    Dist.trees <-bind_rows(Dist.trees, tr.i )
  }
  return(Dist.trees)
}


# Get NHPP data fro stan across all tips
#make_data_NHPP(Tb.trees)

make_data_NHPP <- function(Tb.trees){
  Ntips <- Tb.trees$tip.id %>% unique() %>% length()
  
  dt.out <- list()
  for (i in 1:Ntips){
    dt.out[[i]] <- make_data_NHPP_over_tip(Tb.trees, tip_id=i)
  }
  
  return(dt.out)
}




#   ____________________________________________________________________________
#   Some KDE statistics                                                     ####


# pairwise_pearson_corr_per_edge(n.edges=172, stat='loess.lambda.mean.deriv', Edge.KDE)
# pairwise_pearson_corr_per_edge(n.edges=172, stat='Maps.mean.loess.norm', Edge.KDE)

# Pairwise in respect to BRs correlation for each edge using pearson's coff.

pairwise_pearson_corr_per_edge <- function(n.edges=172, stat='loess.lambda.mean', Edge.KDE){
  
  brs <- names(Edge.KDE)
  len <- length(Edge.KDE)
  
  Comb <- combn(c(1:len), 2, simplify = T)
  Comb.names <- rbind(brs[Comb[1,]], brs[Comb[2,]])
  
  M.out <-matrix(0,nrow = n.edges, ncol = ncol(Comb)) 
  M.out <-as_tibble(M.out)
  
  i <- 2
  for (i in 1:ncol(M.out)){
    
    focal.x <- Comb[1,i]
    focal.y <- Comb[2,i]
    x <- Edge.KDE[[focal.x]][[stat]]
    y <- Edge.KDE[[focal.y]][[stat]]
    
    #---
    #cor.test(x[[1]],y[[1]], method = c("pearson") )
    #---
    
    cor.estimate <- mapply(function(x,y) cor.test(x,y, method = c("pearson") )$estimate, 
                           x=x, y=y)
    
    cor.pval <- mapply(function(x,y) cor.test(x,y, method = c("pearson") )$p.value, 
                       x=x, y=y)
    
    # assign 0 to lage pvals
    cor.estimate[cor.pval>0.5]<-0
    M.out[[i]] <- cor.estimate
    
  }
  
  return(list(corr=M.out, pairs=Comb.names) )
  
}


# Pairwise mean sq error for each edge between densities

#M.err <- pairwise_MSE_over_edges(n.edges=172, stat='Maps.mean.loess.norm', Edge.KDE)
pairwise_MSE_over_edges <- function(n.edges=172, stat='Maps.mean.loess.norm', Edge.KDE){
  
  brs <- names(Edge.KDE)
  len <- length(Edge.KDE)
  
  Comb <- combn(c(1:len), 2, simplify = T)
  Comb.names <- rbind(brs[Comb[1,]], brs[Comb[2,]])
  
  M.diff<-matrix(0,nrow = n.edges, ncol = ncol(Comb))
  M.diff <-as_tibble(M.diff)
  
  for (i in 1:ncol(M.diff)){
    
    focal.x <- Comb[1,i]
    focal.y <- Comb[2,i]
    x <- Edge.KDE[[focal.x]][[stat]]
    y <- Edge.KDE[[focal.y]][[stat]]
    
    # difference Mean sq err
    err <- mapply(function(x,y) mean(c(abs(x-y)^2) ), x=x, y=y)
    M.diff[[i]] <-err
    
    
    
  }
  
  return(list(corr=M.diff, pairs=Comb.names) )
  
}

# Make distance matrix for Pairwise mean sq error for each edge between densities

make_dist_M <- function(M.err, stat.f, Edge.KDE, Comb.names){
  len <- length(Edge.KDE)
  brs <- names(Edge.KDE)
  dist <- matrix( nrow=len, ncol=len)
  colnames(dist) <- rownames(dist) <- brs
  
  for (i in 1:length(stat.f)){
    dist[which(Comb.names[2,i]==colnames(dist)), which(Comb.names[1,i]==row.names(dist) ) ] <- stat.f[i]
  }
  
  dist <- as.dist(dist)
  return(dist)
  
}


#   ____________________________________________________________________________
#   Edge Entropy                                                           ####

edge_entropy <- function(n.edges=172, stat='Maps.mean.loess.norm', Edge.KDE){
  
  brs <- names(Edge.KDE)
  len <- length(Edge.KDE)
  
  #Comb <- combn(c(1:len), 2, simplify = T)
  #Comb.names <- rbind(brs[Comb[1,]], brs[Comb[2,]])
  
  Maps.out <- vector(n.edges, mode='list')
  j <- 1
  i <- 1
  for (j in 1:n.edges){
    
    # get the same edge across all brs
    E.all.brs <- vector(len, mode='list')
    for (i in 1:len){
      E.all.brs[[i]] <- Edge.KDE[[i]][[stat]][[j]]
    }
    
    #E.all.brs <-lapply(E.all.brs, function(x) abs(x))
    
    # calculate Entropy
    Norm.const <- Reduce('+', E.all.brs)
    Prob <- lapply(E.all.brs, function(x) x/Norm.const)
    Prob.mod <- lapply(Prob, function(x) abs(x)*log(abs(x)) )
    Maps.out[[j]] <- Reduce('+', Prob.mod)*-1
    #lapply(E.all.brs, length)
    
  }
  
  return(Maps.out)
  
}


#--------------------
# Edge contribution

edge_contribution <- function(n.edges=172, stat='Maps.mean.loess.norm', Edge.KDE){
  
  brs <- names(Edge.KDE)
  len <- length(Edge.KDE)
  
  #Comb <- combn(c(1:len), 2, simplify = T)
  #Comb.names <- rbind(brs[Comb[1,]], brs[Comb[2,]])
  
  M.out<-matrix(0, nrow = n.edges, ncol = len)
  M.out <-as_tibble(M.out)
  colnames(M.out) <- brs
  
  #Maps.out <- vector(n.edges, mode='list')
  j <- 1
  i <- 1
  for (j in 1:n.edges){
    
    # get the same edge across all brs
    E.all.brs <- vector(len, mode='list')
    for (i in 1:len){
      E.all.brs[[i]] <- Edge.KDE[[i]][[stat]][[j]]
    }
    
    Den.sum <- lapply(E.all.brs, function(x) mean(x) ) %>% unlist # average rate
    M.out[j,] <- Den.sum/sum(Den.sum)
    
  }
  
  return(M.out)
  
}

# Create N Mpas, where N is n(brs)

edge_contribution_maps <- function(n.edges=172, stat='Maps.mean.loess.norm', Edge.KDE, scale=FALSE){
  
  brs <- names(Edge.KDE)
  len <- length(Edge.KDE)
  scale.const <- 1/len
  
  #Comb <- combn(c(1:len), 2, simplify = T)
  #Comb.names <- rbind(brs[Comb[1,]], brs[Comb[2,]])
  
  M.out<-vector(len, mode='list')
  names(M.out) <- brs
  M.out<-lapply(M.out, function(x) x <- vector(n.edges, mode='list') )
  
  
  #Maps.out <- vector(n.edges, mode='list')
  j <- 1
  i <- 1
  for (j in 1:n.edges){
    
    # get the same edge across all brs
    E.all.brs <- vector(len, mode='list')
    for (i in 1:len){
      E.all.brs[[i]] <- Edge.KDE[[i]][[stat]][[j]]
    }
    
    Norm.const <- Reduce('+', E.all.brs)
    Prob <- lapply(E.all.brs, function(x) x/Norm.const)
    
    if (scale==T){
      Prob <-lapply(Prob, function(x) x-scale.const)
    }
    
    for (i in 1:len){
      M.out[[i]][[j]] <- Prob[[i]]
    }
    
    
  }
  
  return(M.out)
  
}
