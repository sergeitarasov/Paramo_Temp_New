stack_stm<-function(stm.list){
  M<-lapply(stm.list, function(x) x$maps)
  M<-lapply(M, function(x) lapply(x, function(y) names(y)))
  M<-Reduce(stack2, M)
  
  M.out<-mapply(function(x,y) 
  {setNames(x, y) },
  x=stm.list[[1]]$maps, y=M )
  
  out<-stm.list[[1]]
  out$maps<-M.out
  return(out)
  
}


#### stack two discrete stm's lists; x,y are the list of state names (i.e. maps)
stack2<-function(x,y){
  mapply(function(x,y) 
  {paste(x,y, sep="") },
  x=x, y=y )
}

# Final stack of maps
# cc chars id to stack
# ntrees number of trees to stack
# dirW directory for zip file

paramo<-function(cc, ntrees=10, dirW=c("") )
{
  tr<-vector("list", ntrees)
  for (i in 1:ntrees){
    
    fl<-paste0(cc, "_", i, ".rds")  
    
    stack.L<-vector("list", length(fl))
    
    for (j in 1:length(fl)){
      
      print(paste0("Reading ", paste0(cc[j], ".zip"), " and ", fl[j]))
      con<-unz(paste0(dirW, cc[j], ".zip"), filename=paste0(dirW, fl[j]) )
      con2 <- gzcon(con)
      stack.L[[j]] <- readRDS(con2)
      close(con)
    }
    
    tr[[i]]<- stack_stm(stack.L)
  }
  return(tr)
}


#   ____________________________________________________________________________
#   Added June 2019                                                         ####

# cc <- paramo.BR1.des$`female genitalia` %>% sub(':', '-', .)
# dir2read <- c("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_discr_new/")
# dir2write <- ("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR1/")
# dirW <- c("") # is the direction inside stack file
# br.name <- c("HAO-test")
# trees2read.in.zip <- 10

# Final stack of maps
# cc chars id to stack
# ntrees number of trees to stack
# dirW directory for zip file


paramo_read_treeID<-function(cc, TreeID=1, dir2read, dirW=c("") )
{
  tr<-vector("list", 1)
  
  fl<-paste0(cc, "_", TreeID, ".rds")  
  
  stack.L<-vector("list", length(fl))
  
  for (j in 1:length(fl)){
    
    print(paste0("Reading ", paste0(cc[j], ".zip"), " and ", fl[j]))
    con<-unz(paste0(dir2read, cc[j], ".zip"), filename=paste0(dirW, fl[j]) )
    con2 <- gzcon(con)
    stack.L[[j]] <- readRDS(con2)
    close(con)
  }
  
  tr[[1]]<- stack_stm(stack.L)
  
  return(tr)
  
}

# 
# cc <- c("C1",   "C4",  "C8")
# dir2read= c("Indv_stack/Discr_maps/")
# dir2write= c("Br_stacks/")
# dirW= c("STEP_5/Discr_maps/") # is the direction inside stack file
# br.name <- c('HAO-test')
# trees2read.in.zip <- 10

#' Title Read discretized maps of individual characters, amalgamate them, and save as one body region character
#'
#' @param cc vector of individual characters (=zip archives) to read 
#' @param dir2read direction where to read them
#' @param dir2write direction where to write body region character
#' @param dirW the direction inside zip zrchive where rds files are located file
#' @param br.name name for the created body region character
#' @param trees2read.in.zip number of trees to read inside zip archive
#'
#' @return
#' @export
#'
#' @examples
#' cc <- c("C1",   "C4",  "C8")
#' dir2read= c("Indv_stack/Discr_maps/")
#' dir2write= c("Br_stacks/")
#' dirW= c("STEP_5/Discr_maps/") # is the direction inside stack file
#' br.name <- c('HAO-test')
#' trees2read.in.zip <- 10


make.chars.BR <- function(
  cc,
  dir2read,
  dir2write,
  dirW = c(""),
  br.name,
  trees2read.in.zip = 1) {
  
  
  for (j in 1:trees2read.in.zip) {
    # read tree
    br.tree <- paramo_read_treeID(cc, TreeID = j, dir2read = dir2read, dirW = dirW)
    
    # write.tree
    saveRDS(br.tree, file = paste0(dir2write, br.name, "_", j, ".rds"))
  }
  
  # putting rds files into archive
  
      # zip works messy if files paths are specified, so we abandon them
            # get current directory
            dir.org <- getwd()
            setwd(dir2write) # set to writing dir to work right there
            # files <- paste0(dir2write, br.name, "_", c(1:trees2read.in.zip), ".rds")
            # zip(zipfile=paste0(dir2write, br.name, ".zip"), files = files, extras = '-i')
            
  files <- paste0(br.name, "_", c(1:trees2read.in.zip), ".rds")
  zip(zipfile=paste0(br.name, ".zip"), files = files)
  file.remove(files)
  
  #get back to the original dir
  setwd(dir.org)

}



#' Title Read rds from zim archive
#'
#' @param zip.name name of zip archive
#' @param rds.prefix  prefix for rds files inside zip
#' @param rds.id file id
#'
#' @return
#' @export
#'
#' @examples
#' zip.name <- 'cranium.zip'
#' rds.prefix <- 'cranium_'
#' 
rds_from_zip <- function(zip.name, rds.prefix, rds.id=1){
  con<-unz(zip.name, filename=paste0(rds.prefix, rds.id, '.rds') )
  con2 <- gzcon(con)
  file.rds <- readRDS(con2)
  return(file.rds)
}


#' Title Get state colors for protting tree when there many states
#'
#' @param tree phyloSimmap
#'
#' @return
#' @export
#'
#' @examples
#' library("RColorBrewer")
#' plotSimmap(tree, get_rough_state_cols(tree),  lwd=3, pts=F,ftype="off", ylim=c(0,100))
get_rough_state_cols <- function(tree) {
  states <- lapply(tree$maps, names) %>%
    unlist() %>%
    unique()
  
  hm.palette <- colorRampPalette(brewer.pal(9, "Set1"), space = "Lab")
  color <- hm.palette(length(states))
  
  return(setNames(color, states))
}




#' Title Merge the same state bins (over one branch) in the desritized stochastic map
#'
#' @param br branch
#'
#' @return
#' @export
#'
#' @examples
#' br=z$maps[[153]]
#' br
#' merge_branch_cat(br)
#' z$maps<-lapply(z$maps, merge_branch_cat)
#' plot(z, setNames(getPalette(length(states)), states))
#' 
merge_branch_cat<-function(br){
  i=2
  while (i<=length(br)){
    if( (names(br[i])) == names( br[i-1] )) {
      br[i-1]<-br[i-1]+br[i]
      br<-br[-i]
    } else{
      i=i+1
    }
  }
  return(br)
}


#' Title Title Merge the same state bins (over a tree) in the desritized stochastic map
#'
#' @param tree tree
#'
#' @return
#' @export
#'
#' @examples
#' 
merge_tree_cat <- function(tree){
  tree$maps<-lapply(tree$maps, merge_branch_cat)
  return(tree)
}


#' Title Splits Simmap tree into n bins (episodic bins)
#'
#' @param tree 
#' @param res nbins
#'
#' @return
#' @export
#'
#' @examples
make_bins_Simmap<-function(tree, res){
  
  
  steps <- 0:res/res * max(phytools:::nodeHeights(tree))
  interv.templ <- cbind(steps[-length(steps)], steps[-1], c(1:(length(steps)-1)) )
  
  H <- phytools:::nodeHeights(tree)
  maps.n <- vector(mode = "list", length = nrow(tree$edge))
  
  
  for (i in 1:nrow(tree$edge)) {
    YY <- cbind(c(H[i, 1], steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))] ), 
                c(steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))], H[i, 2])) -  H[i, 1]
    
    YY.names <- cbind(c(H[i, 1], steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))] ), 
                      c(steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))], H[i, 2]))        #-  H[i, 1]
    
    int.ids <- sapply(YY.names[,2], function(x) which(interv.templ[,2]-x>=0)[1] )
    
    
    TR<-cumsum(tree$maps[[i]])
    
    maps.n[[i]]<-setNames(YY[,2]-YY[,1], int.ids )
  }
  tree$maps<-maps.n
  return(tree)        
}





#   ____________________________________________________________________________
#   Get n changes per bin-edge for episodic model                           ####



#' Title Get sttatistics per bin-edge for episodic model
#'
#' @param tree tree with mapped states (not discritezed, i.e., same categories merged)
#' @param tree.binned binned tree used to scoop changes for bins
#'
#' @return tibble 
#' @export
#'
#' @examples
#' # reda tree from zip
#' tree <- rds_from_zip(zip.name='cranium.zip', rds.prefix='cranium_', rds.id=2)[[1]]
#' # merge discritized bins
#' tree.merge <- merge_tree_cat(tree)
#' # get 10 bins
#' tree.binned <- make_bins_Simmap(tree.merge, res=10)
#' plotSimmap(tree.binned, get_rough_state_cols(tree.binned),  lwd=3, pts=F,ftype="off", ylim=c(0,100))
#' get_bins_stat_Simmap(tree=tree.merge, tree.binned=tree.binned)
#' 
get_bins_stat_Simmap <- function(tree, tree.binned) {
  maps.binned.resc <- tree.binned$maps
  maps.binned.changes <- tree.binned$maps
  maps.focal.resc <- tree$maps


  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### Get number of changes for each bin and write it to maps.binned.changes  ####

  # i <- 145
  for (i in 1:nrow(tree$edge)) {
    maps.focal.resc[[i]] <- c(maps.focal.resc[[i]]) %>% cumsum()
    maps.focal.resc[[i]] <- maps.focal.resc[[i]][-length(maps.focal.resc[[i]])]

    maps.binned.resc[[i]] <- c(maps.binned.resc[[i]]) %>% cumsum()

    int.out <- findInterval(c(maps.focal.resc[[i]]), c(0, maps.binned.resc[[i]]),
      left.open = T, rightmost.closed = FALSE, all.inside = TRUE
    )

    maps.binned.changes[[i]] <- setNames(maps.binned.changes[[i]], rep(0, length(maps.binned.changes[[i]])))


    bins <- table(int.out) %>%
      names() %>%
      as.numeric()
    B <- names(maps.binned.changes[[i]])
    B[bins] <- table(int.out) %>% as.character()
    names(maps.binned.changes[[i]]) <- B
  }

  # tree.out$maps <- maps.binned.changes
  # plotSimmap(tree.out,  lwd=3, pts=F,ftype="off", ylim=c(0,100))

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### Summarize changes as tibble                                           ####

  H <- phytools:::nodeHeights(tree.binned)

  maps.int.id <- tree.binned$maps

  dt.out <- tibble()
  #
  #i <- 171
  for (i in 1:nrow(tree$edge)) {
    res.changes <- c(H[i, 1], maps.binned.changes[[i]]) %>% cumsum()
    n.transition <- maps.binned.changes[[i]] %>%
      names() %>%
      as.numeric()
    inter.id <- maps.int.id[[i]] %>%
      names() %>%
      as.numeric()
    times <- cbind(res.changes[-length(res.changes)], res.changes[-1])

    dt <- tibble(interval_id = inter.id, n_changes = n.transition, t0 = res.changes[-length(res.changes)], t = res.changes[-1])
    dt.out <- bind_rows(dt.out, dt)
  }

  dt.out <- dt.out %>% dplyr::mutate(delta_t = t - t0)
  return(dt.out)
}



#   ____________________________________________________________________________
#   Poission likelihood for episodic model                                  ####


#' Title Poission Likelihood function for Episodic model
#'
#' @param lam lambda parameter
#' @param data output from get_bins_stat_Simmap
#'
#' @return
#' @export
#'
#' @examples
#' nlm(Poiss, p=c(2), data=data)
Poiss <- function(lam, data){
  if (lam<0) return(10^10)
  
  ML <- -1* (dpois(data$n_changes, lambda=lam*data$delta_t)  %>% log %>% sum)
  return(ML)
}



##  ............................................................................
##  Inference for Episodic model                                            ####

#' Title Rum ML inference for episodic model
#'
#' @param dt.out output from get_bins_stat_Simmap
#' @param starting.value starting value for nlm
#'
#' @return estimates of lambda parameter for each bin
#' @export
#'
#' @examples
#' Lambda.est <- runn_episodic_Ln(dt.out, starting.value=1)
#' plot(Lambda.est, type='l')


runn_episodic_Ln <- function(dt.out, starting.value=0.5){
  nbins <- dt.out$interval_id %>% unique() %>% length()
  
  results <- rep(NA, nbins)
  for (i in 1:nbins){
    int <- i
    data <- dt.out %>% dplyr::filter(interval_id==int)
    maxLn <- nlm(Poiss, p=c(starting.value), data=data)
    results[i] <- maxLn$estimate
  }
  return(results)
}

