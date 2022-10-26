library("grImport")
library(ggplot2)
library(RColorBrewer)
library(colorRamps)
library("plotrix")
library(tibble)
library("dplyr")
library(ontologyIndex)

#setwd("~/Documents/My_papers/onto_phylo/data/working/R/June2019")
source('R/Plot_Images_functions.R')
rootdir <- getwd()


#   ____________________________________________________________________________
#   Read in precooked picture                                               ####
#setwd("~/Documents/My_papers/onto_phylo/data/working/Graphics")
setwd("data/wasp_graphics")
# prosessing ps file, this automatically creates "hym-04.ps.xml"
PostScriptTrace("hym-04.ps")
Wasp.img <- readPicture("hym-04.ps.xml")
grid.picture(Wasp.img)
graphics.off()
setwd(rootdir)
#   ____________________________________________________________________________
#   Get Objects                                                             ####




####### GL obj
#setwd("~/Documents/My_papers/onto_phylo/data/working/dataset/Data_final")
#setwd("data/ontology_queries")
C<-read.csv("data/ontology_queries/Terms4Graphics.csv", header=T,  stringsAsFactors = FALSE, na.strings = c("", NA))
C<-as_tibble(C)
GL<-vector("list", nrow(C))
names(GL)<-C$ID

for (i in 1:nrow(C)){
  GL[[i]]$name<-C$Term[i]
  GL[[i]]$layer<-strsplit(C$pic_id[i], ", ")[[1]]
}

GL

#---

# get focal.BR that is two columns: BR name, and BR IRI
# focal.BR for BR1 is obtained by running scripts in NHPP_KDE_BR1.r
rstudioapi::navigateToFile("NHPP_KDE_BR1.r")

focal.BR

#ONT=get_OBO("HAO.obo", extract_tags="everything", propagate_relationships = c("BFO:0000050", "is_a"))
ONT= get_OBO("data/ontology_queries/HAO.obo", extract_tags="everything", propagate_relationships = c("BFO:0000050", "is_a"))

Layers <- lapply(focal.BR$ Level_1_ids, function(x) get_vector_ids_per_term(term =x, ONT, GL))
names(Layers) <- focal.BR$Level_1_names
# make vector
Layers <- setNames(unlist(Layers, use.names=F),rep(names(Layers), lengths(Layers)))


#   ____________________________________________________________________________
#   Scale const (sums of all edges)                                         ####

### Read backbone tree                                                      ####
library(ape)
Tmax <- 280.8679
tree.H <- read.tree('data/timetree/Hymenoptera_br_resolved.tre')
#tree.discr <- discr_Simmap(tree.H, res=2000)
time.sum <- sum(tree.discr$edge.length)


#   ____________________________________________________________________________
#   Get Statistics                                                          ####


# Statistics

# get lambda.post that is rate statisitcs for each body region
# lambda.post for BR1 is obtained by running scripts in NHPP_KDE_BR1.r
rstudioapi::navigateToFile("NHPP_KDE_BR1.r")

Stat <- lapply(lambda.post, function(x) x$Mean) %>% unlist



#   ____________________________________________________________________________
#   Make colors                                                             ####

cols.maps <-make_colors(Stat, palette=rainbow(100, start = 0, end = 0.7) %>% rev)

hm.palette <- colorRampPalette(brewer.pal(9, 'Spectral') %>% rev, space='Lab')
hm.palette <- colorRampPalette(brewer.pal(9, 'Accent') %>% rev, space='Lab')
hm.palette <- colorRampPalette(brewer.pal(10, 'Paired') %>% rev, space='Lab')
hm.palette <- colorRampPalette(brewer.pal(9, 'Blues'), space='Lab')

cols.maps <-make_colors(Stat, palette=hm.palette(100))

#   ____________________________________________________________________________
#   Plot Picture                                                            ####

#---
#Stat
#cols
Layers
cols.maps


new.pic <- make_pic(Wasp.img, Layers, cols.maps)
grid.picture(new.pic)
graphics.off()

# color bar
layout(matrix(c(1:4),ncol=2, nrow=2), heights=c(.7, 1), widths=c(.15, 1))
color.bar(hm.palette(100), min(Stat), max(Stat), ticks=c(min(Stat), max(Stat)/2, max(Stat)) %>% round(.,1),
          title="Rate")

# plot pic
new.pic <- make_pic(Wasp.img, Layers, cols.maps)
grid.picture(new.pic)
graphics.off()






#   ____________________________________________________________________________
#   Plot Different Pics                                                     ####



#   Get Statistics

Stat <- lapply(lambda.post, function(x) x$Mean) %>% unlist

# scale given time
#Stat <- Stat/time.sum

#----
plotSimmap(Trees.focal.BRs[[1]][[1]], get_rough_state_cols(Trees.focal.BRs[[1]][[1]]),  
           lwd=3, pts=F,ftype="off", ylim=c(0,100))
edgelabels()

plotSimmap(Trees.focal.BRs[[1]][[1]], get_rough_state_cols(Trees.focal.BRs[[1]][[1]]),  
           lwd=3, pts=F,fsize=.5, ylim=c(0,100))

# Get stat from focal edges


focal.edge <- 7
Stat <- lapply(Edge.KDE, function(x) x$loess.lambda.mean[[focal.edge]] %>% mean )
Stat <- unlist(Stat)

#--
#focal.edge <- c(2,7,8,9,48, 172, 92, 13)
focal.edge <- c(1:172)
Stat.list <- vector('list', length = length(focal.edge))
names(Stat.list) <- focal.edge

for (i in 1:length(focal.edge)){
  Stat.list[[i]] <- lapply(Edge.KDE, function(x) x$loess.lambda.mean[[focal.edge[i]]] %>% mean ) %>% unlist
}

#   ____________________________________________________________________________
#   Make colors                                                             ####

cols.maps <-make_colors(Stat, palette=rainbow(100, start = 0, end = 0.7) %>% rev)

hm.palette <- colorRampPalette(brewer.pal(9, 'Spectral') %>% rev, space='Lab')
# hm.palette <- colorRampPalette(brewer.pal(9, 'Accent') %>% rev, space='Lab')
# hm.palette <- colorRampPalette(brewer.pal(10, 'Paired') %>% rev, space='Lab')
# hm.palette <- colorRampPalette(brewer.pal(9, 'Blues'), space='Lab')

cols.maps <-make_colors(Stat, palette=hm.palette(100))

#   ____________________________________________________________________________
#   Plot Picture                                                            ####

#---
#Stat
#cols
Layers
cols.maps

#cols.maps <-make_colors(Stat, palette=hm.palette(100))
Max <-  max(Stat.list%>%unlist)
Min <-  min(Stat.list%>%unlist)

#cols.maps <-make_colors_relative_scale(Stat, palette=hm.palette(100), lims=c(Min,Max))

#setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR1/Rate_BRs_Images")
setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR1/Rates_BRs_Images_all_edges")

focal.edge
for (i in 1:length(focal.edge)){
  
  pdf(paste0('BRs_Images_edge_', focal.edge[i],'.pdf'), width=8, height=11)
  # color bar
  layout(matrix(c(1:4),ncol=2, nrow=2), heights=c(.7, 1), widths=c(.2, 1))
  
  # plot pic
  cols.maps <-make_colors_relative_scale(Stat.list[[i]], palette=hm.palette(100), lims=c(Min,Max))
  new.pic <- make_pic(Wasp.img, Layers, cols.maps)
  grid.picture(new.pic)
  dev.off()
  
}


# Plot one color bar


pdf(paste0('BRs_Images_edge_BAR_', focal.edge[i],'.pdf'), width=8, height=11)

# color bar
layout(matrix(c(1:4),ncol=2, nrow=2), heights=c(.7, 1), widths=c(.2, 1))

color.bar(hm.palette(100), min(Stat), max(Stat), 
          ticks=c(min(Stat), max(Stat)/2, max(Stat))%>% round(.,2),
          title="")

# plot pic
cols.maps <-make_colors_relative_scale(Stat.list[[i]], palette=hm.palette(100), lims=c(Min,Max))
new.pic <- make_pic(Wasp.img, Layers, cols.maps)
grid.picture(new.pic)
dev.off()

#graphics.off()