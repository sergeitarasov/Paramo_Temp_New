
#   ____________________________________________________________________________
#   MultiDim Scaling throu time                                             ####

#setwd("~/Documents/My_papers/onto_phylo/data/working/R/June2019")
source('R/Multidim_Scaling_Functions.R')
source("R/R_PARAMO/Functions_Stack_maps.R")
library("ggplot2")
library("stringdist")
library("RColorBrewer")
library("gganimate")
library(dplyr)
library(purrr)

rootdir <- getwd()

##  ............................................................................
##  Test using simulated tree                                               ####

# setwd("~/Documents/My_papers/onto_phylo/data/working/R/June2019")
# #saveRDS(tree.test, file = 'morphospace_tree_test.rds')
# tree.test <- readRDS(file = 'morphospace_tree_test.rds')
# tree.merge <- tree.test
# 
# plot(tree.merge, get_rough_state_cols(tree.merge))
# 
# Md <- MD_morpho(tree.merge)
# add.noise <- c(.1,.1)
# Md.n <- add_noise(Md, add.noise)
# 
# # add groups to points
# Md.n$Points <- mutate(Md.n$Points, tip.id=c(1:nrow(Md.n$Points)))
# 
# gg <- ggplot(Md.n$Points)+
#   geom_segment(data=Md.n$Lines, aes(x = start.V1, y = start.V2, 
#                                      xend = end.V1, yend = end.V2), 
#                colour = "black", size=.5, linetype=1, alpha=.5)+
#   geom_point( aes(x=V1, y=V2, color=time, group=tip.id), alpha=1, size=3) +
#   scale_colour_gradientn(colours = terrain.colors(10)[1:8])
# 
# 
# gg.anim <- gg + transition_reveal( time, keep_last = TRUE) 
# gg.anim



############ For Entire Phenotype

setwd("data/Run1_all/stm_R500_BR/EP")
tree.list <-lapply(c(1:2), function(x) rds_from_zip(zip.name='whole_organism.zip', rds.prefix='whole_organism_', rds.id=x)[[1]])
tree.list <-lapply(tree.list, function(x) merge_tree_cat(x) )
tree.merge <- tree.list[[1]]
setwd(rootdir)

Md <- MD_morpho(tree.merge)
add.noise <- c(.3,.3)
Md.n <- add_noise(Md, add.noise)
#Md.n <- Md

Md.n$Points <- mutate(Md.n$Points, tip.id=c(1:nrow(Md.n$Points)))

gg <- 
  ggplot(Md.n$Points)+
  geom_segment(data=Md.n$Lines, aes(x = start.V1, y = start.V2, 
                                    xend = end.V1, yend = end.V2), 
               colour = "red", size=.3, linetype=1, alpha=.7) + 
  geom_point( aes(x=V1, y=V2, color=time, group=tip.id), alpha=.7, size=1) +
  #scale_colour_gradientn(colours = terrain.colors(10)[1:8])+
    scale_colour_gradient(low = "purple",  high = "white", breaks = seq(0, 280, 50), labels = seq(0, 280, 50)%>%rev)+
    theme_black()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab("Coor1") + ylab("Coor2")+
    guides(color = guide_colourbar(title ='Time Myr' ,reverse =F, draw.ulim = T, draw.llim = T))
    
  ####

  ####
gg
#ggsave("output/phenome_flux.png", units = c("mm"), width=200, height=150)

# insert ggplot code
dev.off()

gg.anim <- gg + transition_reveal( time, keep_last = TRUE)+ labs(title = "Myr: {abs(as.integer(frame_along)-281)}")
animate(gg.anim, height = 500, width =600, nframes=100, res=100)

setwd(rootdir)
anim_save("output/Phenotype_flux.gif")



