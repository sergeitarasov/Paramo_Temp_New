# This script plots something like edge profiles. Not reaally needed for PARAMO

#   ____________________________________________________________________________
#    distance change                                                      ####

#setwd("~/Documents/My_papers/onto_phylo/data/working/R/June2019")
source('R/Multidim_Scaling_Functions.R')
source('R/Hamming_Dist_Functions.R')

#setwd("~/Documents/My_papers/onto_phylo/data/working/R/PARAMO")
source("R/R_PARAMO/Functions_Stack_maps.R")

library("ggplot2")
library("stringdist")
library("RColorBrewer")
library("ape")
library("phytools")
library("phangorn")
library(tibble)
library(dplyr)

#   ____________________________________________________________________________
#   Simulated Data                                                          ####

# setwd("~/Documents/My_papers/onto_phylo/data/working/R/June2019")
# #saveRDS(tree.test, file = 'morphospace_tree_test.rds')
# tree.test <- readRDS(file = 'morphospace_tree_test.rds')
# tree.merge <- tree.test
# 
# plot(tree.merge, get_rough_state_cols(tree.merge))


##  ............................................................................
##  Test using real data                                                    ####

# setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR2")
# tree.list <-lapply(c(1:3), function(x) rds_from_zip(zip.name='wing.zip', rds.prefix='wing_', rds.id=x)[[1]])
# tree.list <-lapply(tree.list, function(x) merge_tree_cat(x) )

#setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/EP")
setwd("data/Run1_all/stm_R500_BR/EP")
tree.merge <-rds_from_zip(zip.name='whole_organism.zip', rds.prefix='whole_organism_', rds.id=10)[[1]]
tree.merge <-merge_tree_cat(tree.merge)

plot(tree.merge, get_rough_state_cols(tree.merge))

Tb.d <- path_hamming_over_tips(tree.merge)
#filter(Tb.d, tip.id==1)$States %>% unique

ggplot()+
  geom_line(data=Tb.d, aes(x=t.end, y=Ham.dist, group=tip.id), alpha=0.5)


#   ____________________________________________________________________________
#   Many trees                                                              ####

setwd("data/Run1_all/stm_R500_BR/EP")
tree.list <-lapply(c(1:1), function(x) rds_from_zip(zip.name='whole_organism.zip', rds.prefix='whole_organism_', rds.id=x)[[1]])
tree.list <-lapply(tree.list, function(x) merge_tree_cat(x) )
closeAllConnections()
#tree.merge <- tree.list[[1]]

Tb.trees <- path_hamming_over_trees(tree.list)
ggplot()+
  geom_line(data=Tb.trees, aes(x=t.end, y=Ham.dist, group=tree.tip.id), alpha=0.08)

ggplot()+
  geom_line(data=Tb.trees, aes(x=t.end, y=Pois.count, group=tree.tip.id), alpha=0.05)


ttt <- filter(Tb.trees, tip.id==87)
hist(ttt$t.start, breaks = 30)

ggplot(ttt, aes(x=t.start, group=tree.id))
  geom_density() 

ggplot(filter(Tb.trees, tip.id==6), aes(x=t.start, group=tree.tip.id 
                                        ))+
  geom_histogram(binwidth = 5)
  geom_density()
colour=alpha("red", 0.1)

ggplot(filter(Tb.trees, tip.id==87), aes(x=t.start,  
                                   ))+
  geom_histogram(binwidth = 20)
  #geom_density( colour=alpha("red", 0.2))
