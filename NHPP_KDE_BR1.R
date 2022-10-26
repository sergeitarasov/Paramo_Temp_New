
#   ____________________________________________________________________________
#   NHPP KDE for BR1                                    ####

#setwd("~/Documents/My_papers/onto_phylo/data/working/R/June2019")
source('R/Multidim_Scaling_Functions.R')
source('R/Hamming_Dist_Functions.R')

#setwd("~/Documents/My_papers/onto_phylo/data/working/R/PARAMO")
source("R/R_PARAMO/Functions_Stack_maps.R")

#setwd("~/Documents/My_papers/onto_phylo/data/working/R/June2019")
source("R/NHPP_Functions.R")
source('R/gggeo_scale.R')

library("dplyr")
library("magrittr")
library("geiger")
library(tidyverse)
library("ggplot2")
library("stringdist")
library("RColorBrewer")
library("ape")
library("phytools")
library("phangorn")
##########################

rootdir <- getwd()

#   ____________________________________________________________________________
#   Make 3 levels of resolution; Name -> HOA-> Fig layers                      ####

#setwd("~/Documents/My_papers/onto_phylo/data/working/dataset/Data_final")
# relates to former F obj
al<-read.csv(file="data/alluvial_plot_data.csv", header=T,  stringsAsFactors = FALSE, na.strings = c("", NA) )
al<-as_tibble(al)
al
# let's in BR names that is necessary for saving files
al$Level_1_names <- sub(' ', '_', al$Level_1_names)
al$Level_2 <- sub(' ', '_', al$Level_2)
# # former C object
# map.layer <- read.csv("Terms4Graphics.csv", header = T, stringsAsFactors = FALSE, na.strings = c("", NA))
# map.layer<-as_tibble(map.layer)

##  ............................................................................
##  Three Levels                                                            ####

# # Body region level 1 terms
# paramo.BR1 <- set_names( al$Level_1_ids, al$Level_1_names)
# Body region level 2 terms
#paramo.BR1 <-al %>% select(Level_1, Level_1_ids) %>% distinct(Level_1, Level_1_ids, .keep_all = FALSE)
paramo.BR1 <-al %>% select(Level_1_names, Level_1_ids) 
focal.BR <- paramo.BR1
# # objects to store data
# episodic.tibbles <- list(BR1=NA, BR1=NA, EP=NA) # store tibble data for each BR
# episodic.EST <- list(BR1=NA, BR1=NA, EP=NA) # # store episodic estimates for each BR


#   ____________________________________________________________________________
#   Read in trees                                                           ####

####
# This tree are already prepered meaning that each BR is stored in zip file after running the 
# published PARAMO pipeline
####

#setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR1")
setwd("data/Run1_all/stm_R500_BR/BR1")

# LIST TO KEEP ALL TREES
Trees.focal.BRs <- vector(length = nrow(focal.BR), mode='list')
names(Trees.focal.BRs) <- focal.BR$Level_1_names

i <- 1
for (i in 1:length(Trees.focal.BRs)){
  Trees.focal.BRs[[i]] <- read_STMP_500trees_batch(prefix=names(Trees.focal.BRs[i]))
}

# chack if all is fine
lapply(Trees.focal.BRs, length)

setwd(rootdir)
#saveRDS(Trees.focal.BRs, 'data/merged_trees/Merged_trees_for_all_BR1.rds')


##  ............................................................................
##  KDE                                                                    ####
#setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR1")
setwd("data/Run1_all/stm_R500_BR/BR1")
#Trees.focal.BRs <- readRDS('Merged_trees_for_all_BR1.rds')
lapply(Trees.focal.BRs, length)
plotSimmap(Trees.focal.BRs[[1]][[1]], get_rough_state_cols(Trees.focal.BRs[[1]][[1]]),  
           lwd=3, pts=F,ftype="off", ylim=c(0,100))

setwd(rootdir)
### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Read backbone tree                                                      ####
Tmax <- 280.8679
tree.H <- read.tree('data/timetree/Hymenoptera_br_resolved.tre')
tree.discr <- discr_Simmap(tree.H, res=2000)
tree.discr$edge


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Make path_hamming                                                      ####

# LIST TO KEEP ALL Tb.trees

###--- !!!! The line below is ver slow, takes several hrs. Perhaps should be optimized?!
# At least re-written with for loop.
#- Try this to see the speed for one BR
# Trees.focal.BRs[[1]]
# x <- path_hamming_over_trees_KDE(Trees.focal.BRs[[1]])
###-- !!!!

#Tb.trees <- vector(length = nrow(focal.BR), mode='list')
Tb.trees <-lapply(Trees.focal.BRs, function(x) path_hamming_over_trees_KDE(x))
names(Tb.trees) <- focal.BR$Level_1_names
#saveRDS(Tb.trees, 'Tb.trees_for_all_BR1.rds')

# Read previous one
Tb.trees <- readRDS('data/merged_trees/Tb.trees_for_all_BR1.rds')


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Make path_data                                                     ####
Path.data <-lapply(Tb.trees, function(x) make_data_NHPP_KDE_Markov_kernel(x))
#names(Path.data) <- focal.BR$Level_2
lapply(Path.data, length)
#saveRDS(Path.data, 'Path.data_for_all_BR1.rds')
#Path.data <- readRDS('Path.data_for_all_BR1.rds')


#   ____________________________________________________________________________
#   ADD Pseudodata                                                          ####

root.node <- 88
root.edges <- which(tree.discr$edge[,1]==root.node)
#Edge.goups <- list(c(1:171), 172)
Edge.goups <- list(c(1:172))

Pseudo.path.data <-vector(length(Path.data), mode='list')
names(Pseudo.path.data) <- names(Path.data)
i <- 1
for (i in 1:length(Path.data)){
  #Pseudo.data <- list(c(-1*Path.data[[i]][[172]]) )
  Pseudo.data <- c(-1*Path.data[[i]][[172]])
  Pseudo.data <- Pseudo.data[Pseudo.data>-100]
  Pseudo.data <-list(Pseudo.data)
    
  Pseudo.path.data[[i]] <-add_pseudodata(Edge.goups, Pseudo.data, Path.data[[i]])
  #hist(Pseudo.data[[1]])
  #hist(Pseudo.path.data[[1]][[172]])
}

#hist(Pseudo.path.data[[2]][[172]])
hist(Pseudo.path.data$cranium[[150]], breaks=30)

#   ____________________________________________________________________________
#   Estimate bind width                                                     ####
# Bin width

h.nrd <- vector(length(Pseudo.path.data), mode='list')
names(h.nrd) <- names(Pseudo.path.data)

for (i in 1:length(Pseudo.path.data)){
  h.nrd[[i]] <- Estimate_band_W(tree.discr, Pseudo.path.data[[i]], band.width='bw.nrd' )
  #h.nrd <- Estimate_band_W(tree.discr, Pseudo.path.data, band.width='bw.nrd' )
  #h.nrd0<- Estimate_band_W(tree.discr, Pseudo.path.data, band.width='bw.nrd0' )
  #h.ucv <- Estimate_band_W(tree.discr, Pseudo.path.data, band.width='bw.ucv' )
}

h <- lapply(h.nrd, mean)
h

#   ____________________________________________________________________________
#   Get KDE                                                                 ####

# Got down to read precooked Edge.KDE

# normalized KDE
Edge.KDE <- vector(length(Pseudo.path.data), mode='list')
names(Edge.KDE) <- names(Pseudo.path.data)

i <- 3
for (i in 1:length(Pseudo.path.data)){
#for (i in 2:length(Pseudo.path.data)){
  #for (i in 4:length(Pseudo.path.data)){
  Edge.KDE[[i]] <- estimate_edge_KDE(tree.discr, Path.data=Pseudo.path.data[[i]] , h=h[[i]])
}

#KDE_test <- estimate_edge_KDE(tree.discr, Path.data=Pseudo.path.data[[2]] , h=h[[2]])

# check if integral over all edges equals 1
for (i in 1:length(Edge.KDE)){
  pr <- integrate_edge_KDE(tree.discr, Edge.KDE[[i]]$Maps.mean.norm)
  cat(pr, '\n')
}

str( Edge.KDE)
Edge.KDE$head


#   ____________________________________________________________________________
#   Load Data                                                               ####

### Read backbone tree                                                      ####
Tmax <- 280.8679
tree.H <- read.tree('data/timetree/Hymenoptera_br_resolved.tre')
tree.discr <- discr_Simmap(tree.H, res=2000)
tree.discr$edge

#setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR1")
#saveRDS(Edge.KDE, 'Edge.KDE_BR1.rds')
Edge.KDE <- readRDS('data/merged_trees/Edge.KDE_BR1.rds')
Trees.focal.BRs <- readRDS('data/merged_trees/Merged_trees_for_all_BR1.rds')

#   ____________________________________________________________________________
#   Loess                                                                   ####

#Edge.KDE.loess <- vector(length(Pseudo.path.data), mode='list')
#names(Edge.KDE.loess) <- names(Pseudo.path.data)

for (i in 1:length(Edge.KDE)){
  Edge.KDE[[i]]$Maps.mean.loess <- Loess_smooting_KDE(tree.discr, Edge.KDE[[i]])
  Edge.KDE[[i]]$Maps.mean.loess.norm <- normilize_KDE(tree.discr, Edge.KDE[[i]]$Maps.mean.loess)
}

# KDE_test.loes <- list()
# KDE_test.loes$Maps.mean.loess <-  Loess_smooting_KDE(tree.discr, KDE_test)
# KDE_test.loes$Maps.mean.loess.norm <- normilize_KDE(tree.discr, KDE_test.loes$Maps.mean.loess)

# check integral
lapply(Edge.KDE, function(x) integrate_edge_KDE(tree.discr, x$Maps.mean.loess.norm))


#   ____________________________________________________________________________
#   Analytical Posterior                                                    ####

lambda.post <- vector(length(Edge.KDE), mode='list')
names(lambda.post) <- names(Edge.KDE)

for (i in 1:length(Edge.KDE)){
  lambda.post[[i]] <- posterior_lambda_KDE(Trees.focal.BRs[[i]])
}

lapply(lambda.post, function(x) x$Mean)


# Get Distriubtions Analytical Posterior
lambda.distr <- tibble()
#names(lambda.post) <- names(Edge.KDE)

for (i in 1:length(Edge.KDE)){
  lambda.distr <- bind_rows(lambda.distr, 
                            posterior_lambda_KDE_Distr(Trees.focal.BRs[[i]], 
                                    n.sim=1000, BR.name=names(Trees.focal.BRs[i]) ) )
}

lambda.distr

# get color legend
source('R/Plot_Images_functions.R')
Stat <- lapply(lambda.post, function(x) x$Mean) %>% unlist
hm.palette <- colorRampPalette(brewer.pal(9, 'Spectral') %>% rev, space='Lab')
cols.maps <-make_colors(Stat, palette=hm.palette(100))
cols.boxpl <- cols.maps[order(Stat)%>% rev]

# Plot Boxes
library(forcats)
ggplot(lambda.distr, aes(x=fct_reorder(BR, sim, 
                                       .fun = mean, .desc =TRUE), y=sim, color=BR)) + 
geom_boxplot(size=2)+
  scale_color_manual(values=cols.boxpl)+
  #scale_fill_brewer(palette="Dark2")+
  xlab('Body region')+ylab('rate')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        #panel.background = element_rect(fill = "transparent",colour = NA), 
        #plot.background = element_rect(fill = "transparent",colour = NA),
        legend.position="none") 
  


#   ____________________________________________________________________________
#   Get Posteriors                                                          ####

#Edge.KDE$head

# Posterioirs and derivative
for (i in 1:length(Edge.KDE)){
  Edge.KDE[[i]]$lambda.mean<- make_postPois_KDE(Edge.KDE[[i]]$Maps.mean.norm, lambda.post[[i]], lambda.post.stat='Mean')
  Edge.KDE[[i]]$loess.lambda.mean<- make_postPois_KDE(Edge.KDE[[i]]$Maps.mean.loess.norm, lambda.post[[i]], lambda.post.stat='Mean')
  
  Edge.KDE[[i]]$loess.lambda.mean.deriv <-derivative_KDE(tree.discr, Edge.KDE[[i]]$loess.lambda.mean)
}

#KDE_test.loes$loess.lambda.mean <- make_postPois_KDE(KDE_test.loes$Maps.mean.loess.norm, lambda.post[[2]], lambda.post.stat='Mean')

#saveRDS(Edge.KDE, 'Edge.KDE_BR1.rds')



# Make Maps
#Edge.KDE$cranium

NHPP.Map <- vector(length(Edge.KDE), mode='list')
names(NHPP.Map) <- names(Edge.KDE)

for (i in 1:length(Edge.KDE)){
  NHPP.Map[[i]]$loess.lambda.mean <- make_contMap_KDE(tree.discr, Edge.KDE[[i]]$loess.lambda.mean)
  NHPP.Map[[i]]$loess.lambda.mean.deriv <- make_contMap_KDE(tree.discr, Edge.KDE[[i]]$loess.lambda.mean.deriv)
  #NHPP.Map[[i]]$loess.lambda.mean.deriv <- make_contMap_KDE(tree.discr, Edge.KDE[[i]]$loess.lambda.mean.deriv)
}

#saveRDS(NHPP.Map, 'NHPP.Map_BR1.rds')

# Plot Simmap

plot.contMap(NHPP.Map$cranium$loess.lambda.mean, lwd=2, outline=F, legend=F, ftype="off", outline=F)
plot.contMap(NHPP.Map$head$loess.lambda.mean.deriv, lwd=2, outline=F, legend=F, ftype="off", outline=F)

plot.contMap(NHPP.Map$wing$loess.lambda.mean, lwd=2, outline=F, legend=F, ftype="off", outline=F)

#----


#   ____________________________________________________________________________
#   Edge profiles                                                           ####

# This gives error due to edge_profiles4plotting() func. Probably earlier version of tibble is needed?

edge.profs <- vector(length(Edge.KDE), mode='list')
names(edge.profs) <- names(Edge.KDE)

i=1
for (i in 1:length(Edge.KDE)){
  edge.profs[[i]]$loess.lambda.mean <- edge_profiles4plotting(tree.discr, Edge.KDE[[i]]$loess.lambda.mean)
  edge.profs[[i]]$loess.lambda.mean.deriv <- edge_profiles4plotting(tree.discr, Edge.KDE[[i]]$loess.lambda.mean.deriv)
  
  edge.profs[[i]]$Maps.mean.loess.norm <- edge_profiles4plotting(tree.discr, Edge.KDE[[i]]$Maps.mean.loess.norm)
}

#Ed.pr.test<- edge_profiles4plotting(tree.discr, KDE_test.loes$loess.lambda.mean)



#   ____________________________________________________________________________
#   Edge Corrs                                                              ####

M.pers <- pairwise_pearson_corr_per_edge(n.edges=172, stat='loess.lambda.mean.deriv', Edge.KDE)
M.pers <- pairwise_pearson_corr_per_edge(n.edges=172, stat='Maps.mean.loess.norm', Edge.KDE)

# PCA
library(factoextra)
res.pca <- prcomp(M.pers$corr, scale = TRUE)
fviz_eig(res.pca)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,    # Avoid text overlapping
             addEllipses = F
             #geom='point'
)

stat.f <- apply(M.pers$corr, 2, function(x) abs(x)%>%mean)
stat.f <- 1/stat.f
DD <- make_dist_M (M.pers$corr, stat.f, Edge.KDE,  Comb.names=M.pers$pairs)
cdmd.dist <- cmdscale(DD)
plot(cdmd.dist)
text(cdmd.dist[,1], cdmd.dist[,2], rownames(cdmd.dist), cex = 0.6)

M.err <- pairwise_MSE_over_edges(n.edges=172, stat='Maps.mean.loess.norm', Edge.KDE)
stat.f <- apply(M.err$corr, 2, function(x) x%>%mean)
DD <- make_dist_M (M.err$corr, stat.f, Edge.KDE, Comb.names=M.pers$pairs)
cdmd.dist <- cmdscale(DD)
plot(cdmd.dist)
text(cdmd.dist[,1], cdmd.dist[,2], rownames(cdmd.dist), cex = 0.6)


#   ____________________________________________________________________________
#   Edge Entropy                                                           ####

Entr <- c()
Entr$Maps.mean.loess.norm <- edge_entropy(n.edges=172, stat='Maps.mean.loess.norm', Edge.KDE)
Entr$loess.lambda.mean <- edge_entropy(n.edges=172, stat='loess.lambda.mean', Edge.KDE)

unlist(Entr$Maps.mean.loess.norm) %>% hist(., breaks=30)
unlist(Entr$loess.lambda.mean) %>% hist(., breaks=30)

TT.map1 <- make_contMap_KDE(tree.discr, Entr$Maps.mean.loess.norm)
TT.map2 <- make_contMap_KDE(tree.discr, Entr$loess.lambda.mean)

plot.contMap(TT.map1, lwd=2, outline=F, legend=T, ftype="off", outline=F)
plot.contMap(TT.map2, lwd=2, outline=F, legend=T, ftype="off", outline=F)

#---

Edge.contr <- edge_contribution(n.edges=172, stat='Maps.mean.loess.norm', Edge.KDE)
hist(Edge.contr[[5]], breaks=30)
plot(Edge.contr[[6]], Edge.contr[[1]], type='l' )
lines(Edge.contr[[6]], Edge.contr[[2]], type='l', col='red')
cor.test(Edge.contr[[3]], Edge.contr[[5]], method = c("pearson") )

Edge.contr <-mutate(Edge.contr, row=c(1:172) )
mm <- gather(Edge.contr, "key", value = "value", -row)
ggplot(mm, aes(x=row, y=value, group=key, color=key))+
  geom_line()

# Maps
Edge.contr.maps <- edge_contribution_maps(n.edges=172, stat='Maps.mean.loess.norm', Edge.KDE, scale=T)

for (i in 1:length(Edge.KDE)){
  NHPP.Map.contr[[i]]$Maps.mean.loess.norm <- make_contMap_KDE(tree.discr, Edge.contr.maps[[i]])
}

plot.contMap(M.w, lwd=2, outline=F, legend=T, ftype="off", outline=F)

# Profiles
edge.profs.contr<- vector(length(Edge.KDE), mode='list')
names(edge.profs.contr) <- names(Edge.KDE)

for (i in 1:length(Edge.KDE)){
  edge.profs.contr[[i]]$Maps.mean.loess.norm <- edge_profiles4plotting(tree.discr, Edge.contr.maps[[i]])
}

Edge.contr.maps$metasoma %>% unlist() %>% min()
Edge.contr.maps$metasoma[[1]]%>% min()
lapply(Edge.contr.maps$metasoma, min)
#   ____________________________________________________________________________
#   Combine edge profiles                                                   ####


edge.profs.comb.lambda <- tibble() # lambda mean
edge.profs.comb.lambda.deriv <- tibble() # derivative
edge.profs.Maps.mean.loess.norm <- tibble() # density
i <- 1
for (i in 1:length(Edge.KDE)){
  edge.profs.comb.lambda <-bind_rows(edge.profs.comb.lambda,
                                     mutate(edge.profs[[i]]$loess.lambda.mean, BR=names(edge.profs[i]), edge.id.BR=paste0(edge.id,'-', names(edge.profs[i])) )
  )
  
  edge.profs.comb.lambda.deriv <-bind_rows(edge.profs.comb.lambda.deriv,
                                           mutate(edge.profs[[i]]$loess.lambda.mean.deriv, BR=names(edge.profs[i]), edge.id.BR=paste0(edge.id,'-', names(edge.profs[i])) )
  )
  
  edge.profs.Maps.mean.loess.norm <-bind_rows(edge.profs.Maps.mean.loess.norm,
                                              mutate(edge.profs[[i]]$Maps.mean.loess.norm, BR=names(edge.profs[i]), edge.id.BR=paste0(edge.id,'-', names(edge.profs[i])) )
  )
}




#   ____________________________________________________________________________
#   Plotting                                                                ####


my.periods <- data.frame(period = c("Quaternary", "Neogene", "Paleogene", "Cretaceous", "Jurassic", "Triassic", "Permian", "Carboniferous", "Devonian", "Silurian", "Ordovician", "Cambrian", "Ediacaran", "Cryogenian", "Tonian"),
                         max_age = c(2.588, 23.03, 66, 145, 201.3, 252.2, 298.9, 358.9, 419.2, 443.4, 485.4, 541, 635, 720, 1000),
                         min_age = c(0, 2.588, 23.03, 66, 145, 201.3, 252.2, 298.9, 358.9, 419.2, 443.4, 485.4, 541, 635, 720),
                         abbr = c("Q", "N", "Pg", "K", "J", "Tr", "P", "C", "D", "S", "O", "Cm","E","Cr","To"),
                         color = c(rgb(249, 249, 127, maxColorValue = 255),rgb(255, 230, 25, maxColorValue = 255),rgb(253, 154, 82, maxColorValue = 255),rgb(127, 198, 78, maxColorValue = 255),rgb(52, 178, 201, maxColorValue = 255),rgb(129, 43, 146, maxColorValue = 255),rgb(240, 64, 40, maxColorValue = 255),rgb(103, 165, 153, maxColorValue = 255),rgb(203, 140, 55, maxColorValue = 255),rgb(179, 225, 182, maxColorValue = 255),rgb(0, 146, 112, maxColorValue = 255),rgb(127, 160, 86, maxColorValue = 255),rgb(254, 217, 106, maxColorValue = 255),rgb(254, 204, 92, maxColorValue = 255),rgb(254, 191, 78, maxColorValue = 255)),
                         stringsAsFactors = FALSE)

my.periods <- my.periods[1:7,]


#pdf('rates_bar', width = 11, height = 8)

gg <- ggplot(data = edge.profs$cranium$loess.lambda.mean, aes(x = X-Tmax, y =  Y, group = edge.id, color=Y))+
  geom_line(alpha=1, size=.5)+
  scale_color_gradientn(colours = rev(rainbow(5, start = 0, end = .7)) )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = 'none'
  )+
  #xlim(0,282)+
  xlab('time')+ylab('rate')+
  #scale_x_continuous(limits = c(0, 300), breaks=c(seq(0, 280, 50), 280), labels = (c(seq(0, 280, 50), 280)-280) %>% abs )+
  #scale_x_continuous(limits = c(-300, 0), breaks=c(seq(-300, 0, 50)))+
  scale_x_continuous(limits = c(-299, 0), breaks=-1*c(0, my.periods$max_age[2:7])%>% round(0), 
                     labels = c(0, my.periods$max_age[2:7])%>% round(0) )+
  scale_y_continuous(limits = c(-0.02, NA))+
  coord_cartesian( expand = FALSE)
#geom_vline(xintercept = 0, color='blue')+
#geom_vline(xintercept = 281, color='blue')
gggeo_scale(gg, periods=my.periods, neg=T)

#dev.off()


gg <- ggplot(data = Ed.pr.test, aes(x = X-Tmax, y =  Y, group = edge.id, color=Y))+
  geom_line(alpha=1, size=.5)+
  scale_color_gradientn(colours = rev(rainbow(5, start = 0, end = .7)) )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(colour = "black"),
        legend.position = 'none'
  )+
  #xlim(0,282)+
  xlab('time')+ylab('rate')+
  #scale_x_continuous(limits = c(0, 300), breaks=c(seq(0, 280, 50), 280), labels = (c(seq(0, 280, 50), 280)-280) %>% abs )+
  #scale_x_continuous(limits = c(-300, 0), breaks=c(seq(-300, 0, 50)))+
  scale_x_continuous(limits = c(-299, 0), breaks=-1*c(0, my.periods$max_age[2:7])%>% round(0), 
                     labels = c(0, my.periods$max_age[2:7])%>% round(0) )+
  scale_y_continuous(limits = c(-0.02, NA))+
  coord_cartesian( expand = FALSE)
#geom_vline(xintercept = 0, color='blue')+
#geom_vline(xintercept = 281, color='blue')
gggeo_scale(gg, periods=my.periods, neg=T)

#Maps.mean.loess.norm
gg <- ggplot(data = edge.profs$cranium$Maps.mean.loess.norm, aes(x = X-Tmax, y =  Y, group = edge.id, color=Y))+
  geom_line(alpha=1, size=.5)+
  scale_color_gradientn(colours = rev(rainbow(5, start = 0, end = .7)) )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(colour = "black"),
        legend.position = 'none'
  )+
  #xlim(0,282)+
  xlab('time')+ylab('rate')+
  #scale_x_continuous(limits = c(0, 300), breaks=c(seq(0, 280, 50), 280), labels = (c(seq(0, 280, 50), 280)-280) %>% abs )+
  #scale_x_continuous(limits = c(-300, 0), breaks=c(seq(-300, 0, 50)))+
  scale_x_continuous(limits = c(-299, 0), breaks=-1*c(0, my.periods$max_age[2:7])%>% round(0), 
                     labels = c(0, my.periods$max_age[2:7])%>% round(0) )+
  scale_y_continuous(limits = c(NA, NA))+
  coord_cartesian( expand = FALSE)
#geom_vline(xintercept = 0, color='blue')+
#geom_vline(xintercept = 281, color='blue')
gggeo_scale(gg, periods=my.periods, neg=T)


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Combined edges                                                          ####

# LAMBDA
ggplot(data = edge.profs.Maps.mean.loess.norm, aes(x = X-Tmax, y =  Y, group = edge.id.BR, color=Y))+
  geom_line(alpha=1, size=.2)+
  scale_color_gradientn(colours = rev(rainbow(5, start = 0, end = .7)) )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(colour = "black")
        
  )+
  #xlim(0,282)+
  xlab('time')+ylab('rate')+
  #scale_x_continuous(limits = c(0, 300), breaks=c(seq(0, 280, 50), 280), labels = (c(seq(0, 280, 50), 280)-280) %>% abs )+
  #scale_x_continuous(limits = c(-300, 0), breaks=c(seq(-300, 0, 50)))+
  scale_x_continuous(limits = c(-299, 0), breaks=-1*c(0, my.periods$max_age[2:7])%>% round(0), 
                     labels = c(0, my.periods$max_age[2:7])%>% round(0) )+
  scale_y_continuous(limits = c(NA, NA))+
  coord_cartesian( expand = FALSE)+
  #geom_vline(xintercept = 0, color='blue')+
  #geom_vline(xintercept = 281, color='blue')+
  facet_wrap(~ BR, ncol=2)

ggplot(data = edge.profs.Maps.mean.loess.norm, aes(x = X-Tmax, y =  Y, group = edge.id.BR, color=BR))+
  geom_line(alpha=.7, size=.8)+
  #scale_color_gradientn(colours = rev(rainbow(5, start = 0, end = .7)) )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(colour = "black")
        
  )+
  #xlim(0,282)+
  xlab('time')+ylab('rate')+
  #scale_x_continuous(limits = c(0, 300), breaks=c(seq(0, 280, 50), 280), labels = (c(seq(0, 280, 50), 280)-280) %>% abs )+
  #scale_x_continuous(limits = c(-300, 0), breaks=c(seq(-300, 0, 50)))+
  scale_x_continuous(limits = c(-299, 0), breaks=-1*c(0, my.periods$max_age[2:7])%>% round(0), 
                     labels = c(0, my.periods$max_age[2:7])%>% round(0) )+
  scale_y_continuous(limits = c(NA, NA))+
  coord_cartesian( expand = FALSE)
#geom_vline(xintercept = 0, color='blue')+
#geom_vline(xintercept = 281, color='blue')+



# DERIVATIVE
ggplot(data = edge.profs.comb.lambda.deriv, aes(x = X-Tmax, y =  Y, group = edge.id.BR, color=Y))+
  geom_line(alpha=1, size=.5)+
  scale_color_gradientn(colours = rev(rainbow(5, start = 0, end = .7)) )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA), 
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(colour = "black")
        
  )+
  #xlim(0,282)+
  xlab('time')+ylab('rate')+
  #scale_x_continuous(limits = c(0, 300), breaks=c(seq(0, 280, 50), 280), labels = (c(seq(0, 280, 50), 280)-280) %>% abs )+
  #scale_x_continuous(limits = c(-300, 0), breaks=c(seq(-300, 0, 50)))+
  scale_x_continuous(limits = c(-299, 0), breaks=-1*c(0, my.periods$max_age[2:7])%>% round(0), 
                     labels = c(0, my.periods$max_age[2:7])%>% round(0) )+
  scale_y_continuous(limits = c(-0.02, NA))+
  coord_cartesian( expand = FALSE)+
  #geom_vline(xintercept = 0, color='blue')+
  #geom_vline(xintercept = 281, color='blue')+
  facet_wrap(~ BR, ncol=2)


#   ____________________________________________________________________________
#   Aligned Plotting RATES                                                    ####


library(gridBase)
library(ggplot2)
library("grid")

setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR1/Profiles_and_Tree_Rates")
#--- Test
# setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR1")
# tree.test <-lapply(c(1:1), function(x) rds_from_zip(zip.name='female_genitalia.zip', rds.prefix='female_genitalia_', rds.id=x)[[1]])
# closeAllConnections()
# tree.test <-merge_tree_cat(tree.test[[1]])
#---

i <- 1
for (i in 1:length(Edge.KDE)){
  
  BR.now <- names(Edge.KDE[i])
  
  #png("test_rates.png", width = 774, height = 513, res = NA) 
  pdf(paste0(BR.now, '_rate_lambda.pdf'), width = 7, height = 5)
  
  #par(mfrow = c(2,1), mar=c(0,0,0,0), oma=c(0,0,0,0))
  layout(matrix(c(1,2),ncol=1), heights=c(2,1))
  
  plot.contMap(NHPP.Map[[BR.now]]$loess.lambda.mean, lwd=2, outline=F, legend=F, plot=F, ftype="off", mar=c(0.1, 3.45 ,0.1, .35))
  #plotSimmap(tree.test, lwd=2, pts=F, ftype="off", ylim=c(0,100), mar=c(0.1, 3.45 ,0.1, .35)) # t ,l, b, r
  
  
  Pl <- ggplot(data = edge.profs[[BR.now]]$loess.lambda.mean, aes(x = X-Tmax, y =  Y, group = edge.id, color=Y))+
    geom_line(alpha=1, size=.3)+
    scale_color_gradientn(colours = rev(rainbow(5, start = 0, end = .7)) )+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA), 
          plot.background = element_rect(fill = "transparent",colour = NA),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(2.3,.87,.1,0.1), 'cm'), # t ,r, b, l
          legend.position = 'none'
          
    )+
    #xlim(0,282)+
    xlab('time')+ylab('rate')+
    #scale_x_continuous(limits = c(0, 300), breaks=c(seq(0, 280, 50), 280), labels = (c(seq(0, 280, 50), 280)-280) %>% abs )+
    #scale_x_continuous(limits = c(-300, 0), breaks=c(seq(-300, 0, 50)))+
    scale_x_continuous(limits = c(-299, 0), breaks=-1*c(0, my.periods$max_age[2:7])%>% round(0), 
                       labels = c(0, my.periods$max_age[2:7])%>% round(0) )+
    scale_y_continuous(limits = c(-0.02, NA))+
    coord_cartesian( expand = FALSE)
  #geom_vline(xintercept = 0, color='blue')+
  #geom_vline(xintercept = -281, color='blue')
  #Pl <- gggeo_scale(Pl, periods=my.periods, neg=T)
  
  vp <- viewport(height = unit(.5,"npc"), width=unit(1, "npc"),
                 just = c("left",'top'),
                 y = .5, x = 0)
  print(Pl, vp = vp)
  
  dev.off()
  
}



#   ____________________________________________________________________________
#   Aligned Plotting Deriv                                                 ####

setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR1/Profiles_and_Tree_Derivs")

i <- 1
for (i in 1:length(Edge.KDE)){
  
  BR.now <- names(Edge.KDE[i])
  
  #png("test_rates.png", width = 774, height = 513, res = NA) 
  pdf(paste0(BR.now, '_deriv_lambda.pdf'), width = 7, height = 5)
  
  #par(mfrow = c(2,1), mar=c(0,0,0,0), oma=c(0,0,0,0))
  layout(matrix(c(1,2),ncol=1), heights=c(2,1))
  
  plot.contMap(NHPP.Map[[BR.now]]$loess.lambda.mean.deriv, lwd=2, outline=F, legend=F, plot=F, ftype="off", mar=c(0.1, 3.45 ,0.1, .35))
  #plotSimmap(tree.test, lwd=2, pts=F, ftype="off", ylim=c(0,100), mar=c(0.1, 3.45 ,0.1, .35)) # t ,l, b, r
  
  
  Pl <- ggplot(data = edge.profs[[BR.now]]$loess.lambda.mean.deriv, aes(x = X-Tmax, y =  Y, group = edge.id, color=Y))+
    geom_line(alpha=1, size=.3)+
    scale_color_gradientn(colours = rev(rainbow(5, start = 0, end = .7)) )+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA), 
          plot.background = element_rect(fill = "transparent",colour = NA),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(2.3,.87,.1,0.1), 'cm'), # t ,r, b, l
          legend.position = 'none'
          
    )+
    #xlim(0,282)+
    xlab('time')+ylab('rate')+
    #scale_x_continuous(limits = c(0, 300), breaks=c(seq(0, 280, 50), 280), labels = (c(seq(0, 280, 50), 280)-280) %>% abs )+
    #scale_x_continuous(limits = c(-300, 0), breaks=c(seq(-300, 0, 50)))+
    scale_x_continuous(limits = c(-299, 0), breaks=-1*c(0, my.periods$max_age[2:7])%>% round(0), 
                       labels = c(0, my.periods$max_age[2:7])%>% round(0) )+
    scale_y_continuous(limits = c(NA, NA))+
    coord_cartesian( expand = FALSE)
  #geom_vline(xintercept = 0, color='blue')+
  #geom_vline(xintercept = -281, color='blue')
  #Pl <- gggeo_scale(Pl, periods=my.periods, neg=T)
  
  vp <- viewport(height = unit(.5,"npc"), width=unit(1, "npc"),
                 just = c("left",'top'),
                 y = .5, x = 0)
  print(Pl, vp = vp)
  
  dev.off()
  
}


#   ____________________________________________________________________________
#   Edge Contribution plots                                                 ####

#NHPP.Map.contr[[i]]$Maps.mean.loess.norm 


setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR1/Edge_contribution_plots")

i <- 1
for (i in 1:length(Edge.KDE)){
  
  BR.now <- names(Edge.KDE[i])
  
  #png("test_rates.png", width = 774, height = 513, res = NA) 
  pdf(paste0(BR.now, '_edge_contr_lambda.pdf'), width = 7, height = 5)
  
  #par(mfrow = c(2,1), mar=c(0,0,0,0), oma=c(0,0,0,0))
  layout(matrix(c(1,2),ncol=1), heights=c(2,1))
  
  plot.contMap(NHPP.Map.contr[[BR.now]]$Maps.mean.loess.norm, lwd=2, outline=F, legend=F, plot=F, ftype="off", mar=c(0.1, 3.45 ,0.1, .35))
  #plotSimmap(tree.test, lwd=2, pts=F, ftype="off", ylim=c(0,100), mar=c(0.1, 3.45 ,0.1, .35)) # t ,l, b, r
  
  #edge.profs.contr[[i]]$Maps.mean.loess.norm
  Pl <- ggplot(data = edge.profs.contr[[BR.now]]$Maps.mean.loess.norm, aes(x = X-Tmax, y =  Y, group = edge.id, color=Y))+
    geom_line(alpha=1, size=.3)+
    scale_color_gradientn(colours = rev(rainbow(5, start = 0, end = .7)) )+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA), 
          plot.background = element_rect(fill = "transparent",colour = NA),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(2.3,.87,.1,0.1), 'cm'), # t ,r, b, l
          legend.position = 'none'
          
    )+
    #xlim(0,282)+
    xlab('time')+ylab('rate')+
    #scale_x_continuous(limits = c(0, 300), breaks=c(seq(0, 280, 50), 280), labels = (c(seq(0, 280, 50), 280)-280) %>% abs )+
    #scale_x_continuous(limits = c(-300, 0), breaks=c(seq(-300, 0, 50)))+
    scale_x_continuous(limits = c(-299, 0), breaks=-1*c(0, my.periods$max_age[2:7])%>% round(0), 
                       labels = c(0, my.periods$max_age[2:7])%>% round(0) )+
    scale_y_continuous(limits = c(NA, NA))+
    coord_cartesian( expand = FALSE)
  #geom_vline(xintercept = 0, color='blue')+
  #geom_vline(xintercept = -281, color='blue')
  #Pl <- gggeo_scale(Pl, periods=my.periods, neg=T)
  
  vp <- viewport(height = unit(.5,"npc"), width=unit(1, "npc"),
                 just = c("left",'top'),
                 y = .5, x = 0)
  print(Pl, vp = vp)
  
  dev.off()
  
}

