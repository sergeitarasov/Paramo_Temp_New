


#' Title Multidimensianl scaling of states from a stachastim map
#'
#' @param tree.merge stochastim map with merged discrte bins
#' @param add.noise a vector of length 2 or NULL.should nois be added? Useful if there are many the same states which occupy the same point in 2d coordinates. the noise 
#' is calculated as var(V)*add.noise
#'
#' @return list of two tibbles -- points, and lines which correpond to tree branches
#' @export
#'
#' @examples
#' MultiScale_1tree(tree.merge, add.noise=c(.5,.5))
#' MultiScale_1tree(tree.merge, add.noise=NULL)
#' 
MultiScale.simmap <- function(tree.merge){
  
  cat('MultiScale.simmap')
  # ---- recode states
  Maps <- tree.merge$maps
  lengths <- lapply(Maps, function(x) length(x)) %>% unlist()
  
  start <- lengths 
  start <- cumsum(start)+1
  start <-c(1, start[-length(start)])
  end <- cumsum(lengths)
  
  states_id <- mapply(function(x,y) c(x:y), start, end)
  # get Map
  Map.new <-mapply(function(x,y) set_names(x,y), Maps, states_id)
  
  # ---- Mape 2 edges
  Map.new.names <- lapply(Map.new, function(x) names(x) %>% as.numeric  )
  X <- lapply(Map.new.names, function(y) cbind( y[-length(y)], y[-1]) )
  Ed.int <- c()
  for (i in 1:length(X)){Ed.int <-rbind(Ed.int, X[[i]])}
  
  Ed.ext.end <- start #lapply(Map.new.names, function(x) x[1]) %>% unlist()
  Ed.ext.start <-end[match(tree.merge$edge[,1], tree.merge$edge[,2])]
  
  # make root
  root <- Ed.ext.end[is.na(Ed.ext.start)]
  # bind colls without NAs
  Ed.ext <- cbind(Ed.ext.start[!is.na(Ed.ext.start)], Ed.ext.end[!is.na(Ed.ext.start)])
  # add root
  Ed.ext <-rbind(root, Ed.ext)
  
  
  Ed.all <- rbind(Ed.int, Ed.ext) %>% as_tibble()
  
  # get those state which are t-1 for rxtant taxa which are t
  tax.id.org <- c(1:length(tree.merge$tip.label))
  tax.id.new <- end[match(tax.id.org, tree.merge$edge[,2] )]
  #
  
  
  # org states
  #Maps <- tree.merge$maps
  states <- lapply(Maps, function(x) names(x)) %>% unlist
  ham.d <- stringdistmatrix(states, states, method = c("hamming") )
  #colnames(ham.d) <- rownames(ham.d) <- states
  mds1 <-  cmdscale(ham.d, k = 2)
  M <- as_tibble(mds1)
  
  # # add noise to separate points with the same coords
  # if (!is.null(add.noise)){
  #   #add.noise=c(.1,.1)
  #   noise.x <- var(M$V1)*add.noise[1]
  #   noise.y <- var(M$V2)*add.noise[2]
  #   noise <- cbind(runif(nrow(M), -noise.x, noise.x), runif(nrow(M), -noise.y, noise.y))
  #   M <- (M+noise) %>% as_tibble()
  # }
  
  
  # get absolute ages of change and join them with M
  H <- phytools:::nodeHeights(tree.merge)
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <-mapply(function(x,y) x+y, age.loc, H[,1] )
  Age <- lapply(age.glob, function(x) x) %>% unlist
  M <- bind_cols(M, time=Age)
  
  # combine extant/estinct species with M
  spp <- rep('no', nrow(M) )
  spp[tax.id.new] <- 'yes'
  M <-bind_cols(M, sp_extant=spp)
  # mark root (2 taxa since tree is not rooted)
  isroot <- rep('no', nrow(M) )
  isroot[root] <- 'yes'
  M <-bind_cols(M, is_root=isroot)
  #Extant.taxa <- bind_cols(M[tax.id.new,], id=tax.id.new)
  
  # get start and end coordinates for lines
  L.start <- M[Ed.all$V1,]
  L.start <- L.start %>% rename(start.V1=V1, start.V2=V2, start.time=time, start.sp_extant=sp_extant)
  L.end <- M[Ed.all$V2,]
  L.end <-L.end %>% rename(end.V1=V1, end.V2=V2, end.time=time, end.sp_extant=sp_extant)
  Lines.coor <-bind_cols(L.start, L.end) 
  
  MD <- list(Points=M, Lines=Lines.coor, Edge.map=Ed.all)
  class(MD) <- append(class(MD), 'md_1tree')
  return(MD)
}

# # general wrapper
# add_noise <- function(MD, add.noise=c(.1,.1)){
#   if ('tree.id' %in% names(MD$Point))
#     return(add_noise_MD_list(MD, add.noise))
#   
#   if (!('tree.id' %in% names(MD$Point)) )
#     return(add_noise_MD(MD, add.noise))
# }

add_noise <- function(MD, ...){
  UseMethod('add_noise_MD', MD)
}

# general wrapper
MD_morpho <- function(MD){
  UseMethod('MultiScale', MD)
}

#' Title Adds nois to Multi Dim scaling for one tree
#'
#' @param MD 
#' @param add.noise 
#'
#' @return
#' @export
#'
#' @examples

#MD <- MS
add_noise_MD.md_1tree <- function(MD, add.noise){
  cat('Noise Simmap')
  
  # Update Points
  M <- MD$Points  
  # noise.x <- var(M$V1)*add.noise[1]
  # noise.y <- var(M$V2)*add.noise[2]
  #noise <- cbind(runif(nrow(M), -add.noise[1], add.noise[1]), runif(nrow(M), -add.noise[2], add.noise[2]))
  M$V1 <- M$V1+runif(nrow(M), -add.noise[1], add.noise[1])
  M$V2 <- M$V2+runif(nrow(M), -add.noise[2], add.noise[2])
  MD$Points <- M
  
  # Update lines
  L.start <- M[MD$Edge.map$V1,]
  L.start <- L.start %>% rename(start.V1=V1, start.V2=V2, start.time=time, start.sp_extant=sp_extant)
  L.end <- M[MD$Edge.map$V2,]
  L.end <-L.end %>% rename(end.V1=V1, end.V2=V2, end.time=time, end.sp_extant=sp_extant)
  Lines.coor <-bind_cols(L.start, L.end) 
  MD$Lines <- Lines.coor
  
  return(MD)
}

#   ____________________________________________________________________________
#   Working on                                                              ####


##  ............................................................................
##  Helper functions                                                        ####

#tree.id <- 1
# get graph of states connections for 1 tree
get_Ed_all <- function(tree.merge, tree.id){
  # ---- recode states
  Maps <- tree.merge$maps
  lengths <- lapply(Maps, function(x) length(x)) %>% unlist()
  
  start <- lengths 
  start <- cumsum(start)+1
  start <-c(1, start[-length(start)])
  end <- cumsum(lengths)
  
  states_id <- mapply(function(x,y) c(x:y), start, end)
  # get Map
  Map.new <-mapply(function(x,y) set_names(x,y), Maps, states_id)
  
  # ---- Mape 2 edges
  Map.new.names <- lapply(Map.new, function(x) names(x) %>% as.numeric  )
  X <- lapply(Map.new.names, function(y) cbind( y[-length(y)], y[-1]) )
  Ed.int <- c()
  for (i in 1:length(X)){Ed.int <-rbind(Ed.int, X[[i]])}
  
  Ed.ext.end <- start #lapply(Map.new.names, function(x) x[1]) %>% unlist()
  Ed.ext.start <-end[match(tree.merge$edge[,1], tree.merge$edge[,2])]
  
  # make root
  root <- Ed.ext.end[is.na(Ed.ext.start)]
  # bind colls without NAs
  Ed.ext <- cbind(Ed.ext.start[!is.na(Ed.ext.start)], Ed.ext.end[!is.na(Ed.ext.start)])
 
   # add root
  #Ed.ext <-rbind(root, Ed.ext)
  
  Ed.all <-rbind(Ed.int, Ed.ext) %>% as_tibble() %>% mutate(., is.root='no') %>% 
    add_row(V1=root[1], V2=root[2], is.root='root') %>% cbind(.,tree.id)
  
  #Ed.all <- rbind(Ed.int, Ed.ext) %>% cbind(.,tree.id) %>% as_tibble()
  return(Ed.all)
}


# get absolute ages of change and join them with M
states_ages <- function(tree.merge){
  H <- phytools:::nodeHeights(tree.merge)
  Maps <- tree.merge$maps
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <-mapply(function(x,y) x+y, age.loc, H[,1] )
  Age <- lapply(age.glob, function(x) x) %>% unlist
  return(Age)
}

# get absolute ages of change for subsampled tree
states_ages_subsmpl <- function(tree.merge){
  #H <- phytools:::nodeHeights(tree.merge)
  Maps <- tree.merge$maps
  #age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <-Maps
  Age <- lapply(age.glob, function(x) x) %>% unlist
  return(Age)
}

# age.glob <- Maps
# Age <- lapply(age.glob, function(x) x) %>% unlist
# M <- bind_cols(M, time=Age)



##  ............................................................................
##  Multi Scale across trees                                                ####

#tree.merge$maps

#rm(tree.list)
#tree.list <- list(tree.merge)

# MS <- MD_morpho(tree.list)
# tree.merge$maps
# plot(tree.merge, get_rough_state_cols(tree.merge))

MultiScale.list <- function(tree.list){

  cat('MultiScale.list')
  ###
  # GET M  shared among all trees
  
  
  Maps <- lapply(tree.list, function(x) x$maps)
  # get states from a single map
  S <- function(MP){lapply(MP, function(x) names(x)) %>% unlist}
  #S(tree.merge$maps)
  # get maps over list
  states <-lapply(Maps, function(y) S(y)) #%>% unlist
  # get tree ids fro states
  n.states.tr <- lapply(states, function(y) length(y)) %>% unlist
  trs4states <- rep(c(1:length(n.states.tr)), n.states.tr)
  # joins ids and states
  states <-tibble(states=unlist(states), tree.id=trs4states)
  
  # hamming
  ham.d <- stringdist::stringdistmatrix(states$states, states$states, method = c("hamming") )
  mds1 <-  cmdscale(ham.d, k = 2) %>% as_tibble()
  M <- bind_cols( mds1, states)
  
  ####
  
  # get edges matrix for all trees
  Ed.list.all <- c()
  for (i in 1:length(tree.list)){Ed.list.all <-rbind(Ed.list.all, get_Ed_all(tree.list[[i]], tree.id =i)  )}
  Ed.list.all <- Ed.list.all %>% as_tibble()

  
  Age <- lapply(tree.list, function(x) states_ages(x))%>% unlist
  M <- bind_cols(M, time=Age)
  
  # get those state which are t-1 for rxtant taxa which are t
    
    M.glob <- tibble() 
    Lines.coor <-tibble()
    i <- 1
    for (i in 1:length(tree.list)){
      

      
      tree.merge <- tree.list[[i]]
     
       # ---- recode states
      MapsI <- tree.merge$maps
      lengths <- lapply(MapsI, function(x) length(x)) %>% unlist()
      
      #start <- lengths 
      #start <- cumsum(start)+1
      #start <-c(1, start[-length(start)])
      end <- cumsum(lengths)
      ###
      
      Mi <- dplyr::filter(M, tree.id==i)
      # get tips
      tax.id.org <- c(1:length(tree.merge$tip.label))
      tax.id.new <- end[match(tax.id.org, tree.merge$edge[,2] )]
      
      # combine extant/estinct species with M
      spp <- rep('no', nrow(Mi) )
      spp[tax.id.new] <- 'yes'
      Mi <-bind_cols(Mi, sp_extant=spp)
      
      
      # get start and end coordinates for lines
      Ed.all.i <- dplyr::filter(Ed.list.all, tree.id==i) 
      L.start <- Mi[Ed.all.i$V1,]
      L.start <- L.start %>% rename(start.V1=V1, start.V2=V2, start.time=time, start.sp_extant=sp_extant)
      L.end <- Mi[Ed.all.i$V2,]
      L.end <-L.end %>% rename(end.V1=V1, end.V2=V2, end.time=time, end.sp_extant=sp_extant)
      
      Lines.coor <-bind_cols(L.start, L.end) %>% bind_rows(Lines.coor, .)
      
      # root
      Root <-Ed.all.i %>% filter(is.root=='root') %>% select(V1,V2) %>% unlist
      isroot <- rep('no', nrow(Mi) )
      isroot[Root] <- 'yes'
      Mi <-bind_cols(Mi, is_root=isroot)
      
      M.glob <-bind_rows(M.glob, Mi)

      
    }
    
  MD <- list(Points=M.glob, Lines=Lines.coor, Edge.map=Ed.list.all)
  class(MD) <- append(class(MD), 'md_tree_list')
  return(MD)
}

# add noise for list tree objects
#MD <- MSm
add_noise_MD.md_tree_list <- function(MD, add.noise){
  cat('Nois List')
  
  # Update Points
  M <- MD$Points  
  M$V1 <- M$V1+runif(nrow(M), -add.noise[1], add.noise[1])
  M$V2 <- M$V2+runif(nrow(M), -add.noise[2], add.noise[2])
  MD$Points <- M
  
  # Update lines
  uniq.trs <- M$tree.id %>% unique() %>% length()
  Lines.coor <-tibble()
  
  #i <- 1
  for (i in 1:uniq.trs){
    # get start and end coordinates for lines
    Mi <- dplyr::filter(M, tree.id==i)
    Ed.all.i <- dplyr::filter(MD$Edge.map, tree.id==i) 
    L.start <- Mi[Ed.all.i$V1,]
    L.start <- L.start %>% rename(start.V1=V1, start.V2=V2, start.time=time, start.sp_extant=sp_extant)
    L.end <- Mi[Ed.all.i$V2,]
    L.end <-L.end %>% rename(end.V1=V1, end.V2=V2, end.time=time, end.sp_extant=sp_extant)
    
    Lines.coor <-bind_cols(L.start, L.end) %>% bind_rows(Lines.coor, .)
  }
  
  
  # L.start <- M[MD$Edge.map$V1,]
  # L.start <- L.start %>% rename(start.V1=V1, start.V2=V2, start.time=time, start.sp_extant=sp_extant)
  # L.end <- M[MD$Edge.map$V2,]
  # L.end <-L.end %>% rename(end.V1=V1, end.V2=V2, end.time=time, end.sp_extant=sp_extant)
  # Lines.coor <-bind_cols(L.start, L.end) 
  MD$Lines <- Lines.coor
  
  return(MD)
}



#   ____________________________________________________________________________
#   Subsample states                                                        ####

# subsample states and return the tree where maps are scale to absolute ages
# subsampling removes a proportion of states in the middle of a branch



# removes states in moddle
remove.states <- function(len, prop, is.root){
  if (is.root){ return(c(TRUE, runif(n=len-2, 0,1)>=prop, TRUE) )}
  if (!is.root){ return(c(runif(n=len-1, 0,1)>=prop, TRUE)) }
}

# check if edge is a root
is_root <- function(RR){
  if (RR==0) 
    return(TRUE)
  else 
    return(FALSE)
}



#' Title subsample states and return the tree where maps are scale to absolute ages. 
#' subsampling removes a proportion of states in the middle of a branch
#'
#' @param tree.merge tree
#' @param prop random proportion of internal states to remove
#'
#' @return
#' @export
#'
#' @examples
#' tr <- subsample_states(tree.merge, prop=.5)
#' plot(tr, get_rough_state_cols(tr))
#' 
subsample_states.simmap <- function(tree.merge, prop){
  
  cat('Simmap')
  
  Maps <- tree.merge$maps
  H <- phytools:::nodeHeights(tree.merge)
  age.loc <- lapply(Maps, function(x) cumsum(x))
  age.glob <-mapply(function(x,y) x+y, age.loc, H[,1] )
  
  E.new <-vector(mode = 'list', length(age.glob))
  
  for (i in 1:nrow(H)){
    
    E <- age.glob[[i]]
    isR <- is_root(H[i,1])
    
    if (length(E)==1){E.new[[i]] <- E}
    
    if (length(E)>1){
      E.new[[i]] <- E[remove.states(len=length(E), prop, is.root=isR)]
      
    }
    
  }
  tree.new <- tree.merge
  tree.new$maps <- E.new
  class(tree.new) <- append('simmap_subsmpl', class(tree.new))
  return(tree.new)

}

# subsample oer tree list
subsample_states.list <- function(tree.list, prop){
  cat('List')
  tree.new <- lapply(tree.list, function(x) subsample_states.simmap(x, prop))
  class(tree.new) <- append('list_subsmpl', class(tree.new))
  return(tree.new)
}

sample_statesINtree <- function(x, ...){
  UseMethod('subsample_states', x)
}

#   ____________________________________________________________________________
#   MutiDim scaling on subsampled tress                                     ####
#' Title Multidimensianl scaling of states from a stachastim map
#'
#' @param tree.merge subsampled stochastim map with merged discrte bins
#' @param add.noise a vector of length 2 or NULL.should nois be added? Useful if there are many the same states which occupy the same point in 2d coordinates. the noise 
#' is calculated as var(V)*add.noise
#'
#' @return list of two tibbles -- points, and lines which correpond to tree branches
#' @export
#'
#' @examples
#' MultiScale_1tree(tree.merge, add.noise=c(.5,.5))
#' MultiScale_1tree(tree.merge, add.noise=NULL)
#' 
MultiScale.simmap_subsmpl <- function(tree.merge){
  
  cat('MultiScale.simmap_subsmpl')
  
  # ---- recode states
  Maps <- tree.merge$maps
  lengths <- lapply(Maps, function(x) length(x)) %>% unlist()
  
  start <- lengths 
  start <- cumsum(start)+1
  start <-c(1, start[-length(start)])
  end <- cumsum(lengths)
  
  states_id <- mapply(function(x,y) c(x:y), start, end)
  # get Map
  Map.new <-mapply(function(x,y) set_names(x,y), Maps, states_id)
  
  # ---- Mape 2 edges
  Map.new.names <- lapply(Map.new, function(x) names(x) %>% as.numeric  )
  X <- lapply(Map.new.names, function(y) cbind( y[-length(y)], y[-1]) )
  Ed.int <- c()
  for (i in 1:length(X)){Ed.int <-rbind(Ed.int, X[[i]])}
  
  Ed.ext.end <- start #lapply(Map.new.names, function(x) x[1]) %>% unlist()
  Ed.ext.start <-end[match(tree.merge$edge[,1], tree.merge$edge[,2])]
  
  # make root
  root <- Ed.ext.end[is.na(Ed.ext.start)]
  # bind colls without NAs
  Ed.ext <- cbind(Ed.ext.start[!is.na(Ed.ext.start)], Ed.ext.end[!is.na(Ed.ext.start)])
  # add root
  Ed.ext <-rbind(root, Ed.ext)
  
  
  Ed.all <- rbind(Ed.int, Ed.ext) %>% as_tibble()
  
  # get those state which are t-1 for rxtant taxa which are t
  tax.id.org <- c(1:length(tree.merge$tip.label))
  tax.id.new <- end[match(tax.id.org, tree.merge$edge[,2] )]
  #
  
  
  # org states
  #Maps <- tree.merge$maps
  states <- lapply(Maps, function(x) names(x)) %>% unlist
  ham.d <- stringdistmatrix(states, states, method = c("hamming") )
  #colnames(ham.d) <- rownames(ham.d) <- states
  mds1 <-  cmdscale(ham.d, k = 2)
  M <- as_tibble(mds1)
  
  
  
  # get absolute ages of change and join them with M
  #H <- phytools:::nodeHeights(tree.merge)
  #age.loc <- lapply(Maps, function(x) cumsum(x))
  #age.glob <-mapply(function(x,y) x+y, age.loc, H[,1] )
  age.glob <- Maps
  Age <- lapply(age.glob, function(x) x) %>% unlist
  M <- bind_cols(M, time=Age)
  
  # combine extant/estinct species with M
  spp <- rep('no', nrow(M) )
  spp[tax.id.new] <- 'yes'
  M <-bind_cols(M, sp_extant=spp)
  # mark root (2 taxa since tree is not rooted)
  isroot <- rep('no', nrow(M) )
  isroot[root] <- 'yes'
  M <-bind_cols(M, is_root=isroot)
  #Extant.taxa <- bind_cols(M[tax.id.new,], id=tax.id.new)
  
  # get start and end coordinates for lines
  L.start <- M[Ed.all$V1,]
  L.start <- L.start %>% rename(start.V1=V1, start.V2=V2, start.time=time, start.sp_extant=sp_extant)
  L.end <- M[Ed.all$V2,]
  L.end <-L.end %>% rename(end.V1=V1, end.V2=V2, end.time=time, end.sp_extant=sp_extant)
  Lines.coor <-bind_cols(L.start, L.end) 
  
  return(list(Points=M, Lines=Lines.coor) )
}

######################




#tree.list <- list.sub
#<- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- 
MultiScale.list_subsmpl <- function(tree.list){
  
  ###
  # GET M  shared among all trees
  cat('MultiScale.list_subsmpl')
  
  Maps <- lapply(tree.list, function(x) x$maps)
  # get states from a single map
  S <- function(MP){lapply(MP, function(x) names(x)) %>% unlist}
  #S(tree.merge$maps)
  # get maps over list
  states <-lapply(Maps, function(y) S(y)) #%>% unlist
  # get tree ids fro states
  n.states.tr <- lapply(states, function(y) length(y)) %>% unlist
  trs4states <- rep(c(1:length(n.states.tr)), n.states.tr)
  # joins ids and states
  states <-tibble(states=unlist(states), tree.id=trs4states)
  
  # hamming
  ham.d <- stringdistmatrix(states$states, states$states, method = c("hamming") )
  mds1 <-  cmdscale(ham.d, k = 2) %>% as_tibble()
  M <- bind_cols( mds1, states)
  
  
  # get edges matrix for all trees
  Ed.list.all <- c()
  for (i in 1:length(tree.list)){Ed.list.all <-rbind(Ed.list.all, get_Ed_all(tree.list[[i]], tree.id =i)  )}
  
  # age.glob <- Maps
  # Age <- lapply(age.glob, function(x) x) %>% unlist
  # M <- bind_cols(M, time=Age)
  
  
  Age <- lapply(tree.list, function(x) states_ages_subsmpl(x))%>% unlist
  M <- bind_cols(M, time=Age)
  
  # get those state which are t-1 for rxtant taxa which are t
  
  M.glob <- tibble() 
  Lines.coor <-tibble()
  
  i <- 1
  for (i in 1:length(tree.list)){
    
    tree.merge <- tree.list[[i]]
    Mi <- dplyr::filter(M, tree.id==i)
    # get tips
    tax.id.org <- c(1:length(tree.merge$tip.label))
    tax.id.new <- end[match(tax.id.org, tree.merge$edge[,2] )]
    
    # combine extant/estinct species with M
    spp <- rep('no', nrow(Mi) )
    spp[tax.id.new] <- 'yes'
    Mi <-bind_cols(Mi, sp_extant=spp)
    
    # mark root (2 taxa since tree is not rooted)
    isroot <- rep('no', nrow(Mi) )
    isroot[root] <- 'yes'
    Mi <-bind_cols(Mi, is_root=isroot)
    
    M.glob <-bind_rows(M.glob, Mi)
    
    # make root
    #root <- Ed.ext.end[is.na(Ed.ext.start)]
    
    # get start and end coordinates for lines
    Ed.all.i <- dplyr::filter(Ed.list.all, tree.id==i) 
    L.start <- Mi[Ed.all.i$V1,]
    L.start <- L.start %>% rename(start.V1=V1, start.V2=V2, start.time=time, start.sp_extant=sp_extant)
    L.end <- Mi[Ed.all.i$V2,]
    L.end <-L.end %>% rename(end.V1=V1, end.V2=V2, end.time=time, end.sp_extant=sp_extant)
    
    Lines.coor <-bind_cols(L.start, L.end) %>% bind_rows(Lines.coor, .)
    
  }
  
  MD <- list(Points=M.glob, Lines=Lines.coor, Edge.map=Ed.list.all)
  class(MD) <- append(class(MD), 'md_tree_list')
  return(MD)
}


#   ____________________________________________________________________________
#   Instruction for MultiDim Sacaling functions                             ####
# 
# MultiScale_1tree(tree.merge) # works on 1 tree class(MD)  append(class(MD), 'tree1'); class simmap
# MultiScale_trees(tree.list) # works on the list of trees; class list
# MultiScale_subspl_1tree(tree.merge) # works on the subsampled tree (in subsampled tr map edges are in absolute ages);  class tree_subsmpl, simmap
# MultiScale_subspl_trees(tree.list) # class tree_subsmpl, list
# 
# add_noise_MD_list(MD, add.noise) 
# add_noise_MD(MD, add.noise)
# 
# class(tree.list)
# class(tree.merge)
# ttt <- c(10:20)
# class(ttt) <- append(class(ttt), c('c2', 'c1') )
# 
# printM <- function(x) UseMethod('printM', x)
# printM.c1 <- function(x) print(x[1])
# printM.c2 <- function(x) print(x[2])
# printM.c1.c2 <- function(x) print(x[3])
# 
# printM(ttt)
# class(tree.list)
# class(tree.merge)


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Function ggplot ThemeBlack                                              ####

theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
      panel.margin = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}
