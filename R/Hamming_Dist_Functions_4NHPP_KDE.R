#   ____________________________________________________________________________
#                  ####
#phangorn::Descendants(tree.discr, node=Node, type = c("tips"))[[1]]


#tree.merge$edge
# get descnedant edges ids from node
#tree.merge <- tree.list[[1]]
# Ancestors(tree.merge, node=90, type = c("all"))
# tree.merge$edge
# get_descehnpath_edges <- function(tree.merge, node){
#   E <- tree.merge$edge
#   E.des <- c()
#   
#   while (length(node)>0){
#     e.id <- which(E[,2]==node)
#     E.des <-c(E.des,e.id)
#     node <- E[e.id, 1]
#   }
#   
#   return(E.des)
# }

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

#tree.merge <- tree.list[[1]]
path_hamming_over_all_edges <- function(tree.merge){
  #n.tips <- length(tree.merge$tip.label)
  #Nodes <- tree.merge$edge[,2] %>% unique()
  #n.tips <- nrow(tree.merge$edge)
  Dist.tips <- tibble()
  i <- 3
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
#tree.list <- Trees.focal.BRs$head
#tree.list <-Trees.focal.BRs$female_genitalia
path_hamming_over_trees_KDE <- function(tree.list){
  
  Dist.trees <- tibble()
  i <- 1
  for (i in 1:length(tree.list)){
    #---
    cat('Working on tree', i, '\n')
    #---
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