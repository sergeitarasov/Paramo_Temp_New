#   ____________________________________________________________________________
#   NHPP using Kernel denisty Estimates  for Entire Phenotype               ####

#setwd("~/Documents/My_papers/onto_phylo/data/working/R/June2019")
source('R/Multidim_Scaling_Functions.R')
source('R/Hamming_Dist_Functions.R')
source('R/Hamming_Dist_Functions_4NHPP_KDE.R')

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

# Read data by parts cuz they are large
setwd("data/Run1_all/stm_R500_BR/EP")
tree.list1 <-lapply(c(1:100), function(x) rds_from_zip(zip.name='whole_organism.zip', rds.prefix='whole_organism_', rds.id=x)[[1]])
closeAllConnections()
tree.list2 <-lapply(c(101:200), function(x) rds_from_zip(zip.name='whole_organism.zip', rds.prefix='whole_organism_', rds.id=x)[[1]])
closeAllConnections()
tree.list3 <-lapply(c(201:300), function(x) rds_from_zip(zip.name='whole_organism.zip', rds.prefix='whole_organism_', rds.id=x)[[1]])
closeAllConnections()
tree.list4 <-lapply(c(301:400), function(x) rds_from_zip(zip.name='whole_organism.zip', rds.prefix='whole_organism_', rds.id=x)[[1]])
closeAllConnections()
tree.list5 <-lapply(c(401:500), function(x) rds_from_zip(zip.name='whole_organism.zip', rds.prefix='whole_organism_', rds.id=x)[[1]])
closeAllConnections()

setwd(rootdir)

tree.list <-c(tree.list1, tree.list2, tree.list3, tree.list4, tree.list5)
length(tree.list)
tree.list <-lapply(tree.list, function(x) merge_tree_cat(x) )
closeAllConnections()

plotSimmap(tree.list[[1]], get_rough_state_cols(tree.list[[1]]),  lwd=3, pts=F,ftype="off", ylim=c(0,100))
edgelabels()
nodelabels()
#Ancestors(tree.list[[1]], node=96, type = c("all"))
#tree.list[[1]]$edge

#------ !!!! This is slow!!!!
#Tb.trees <- path_hamming_over_trees(tree.list)
Tb.trees <-path_hamming_over_trees_KDE(tree.list)
#Tb.trees$Edge.id %>% unique()
#saveRDS(Tb.trees, "Tb.trees_500_EP.rds")
Tb.trees <- readRDS("data/merged_trees/Tb.trees_500_EP.rds")
#------ !!!! This is slow!!!!

# Make data fro all tips
#Path.data <- make_data_NHPP(Tb.trees)
Path.data <-make_data_NHPP_KDE_Markov_kernel(Tb.trees)
#saveRDS(Path.data, "Path.data_500_EP.rds")

Tmax <- 280.8679
tree.H <- read.tree('data/timetree/Hymenoptera_br_resolved.tre')
tree.discr <- discr_Simmap(tree.H, res=2000)
tree.discr$edge

#---------------
# ADD PSEUDODATA
# ST: as far as I remeber pseudodata (unigorm) are edded either to pre-root or post-tips to make tha part without changes 
root.node <- 88
root.edges <- which(tree.discr$edge[,1]==root.node)
Edge.goups <- list(c(1:171), 172)

# For edge group 1
dt1 <- Path.data[[1]]
hist(dt1)
psd1 <- runif(-60,0,n=round((60/3.1)*length(dt1),0) )
#psd1 <- lapply(cumsum(c(0, rep(3.1, 20)) ), function(x) x+dt1) %>% unlist()*-1
hist(c(psd1, dt1), breaks = 50)

Pseudo.data <- list(
  psd1,
  rep((-1*Path.data[[172]]))
)

Pseudo.path.data <-add_pseudodata(Edge.goups, Pseudo.data, Path.data)
hist(Pseudo.path.data[[14]], breaks=50)

#--------------

# Bin width
h.nrd <- Estimate_band_W(tree.discr, Pseudo.path.data, band.width='bw.nrd' )
h.nrd0<- Estimate_band_W(tree.discr, Pseudo.path.data, band.width='bw.nrd0' )
h.ucv <- Estimate_band_W(tree.discr, Pseudo.path.data, band.width='bw.ucv' )

h <- mean(h.nrd)
h
# 8.851224

# normalized KDE
Edge.KDE <- estimate_edge_KDE(tree.discr, Path.data=Pseudo.path.data , h=h)
#Edge.KDE <- estimate_edge_KDE(tree.discr, root.taxon=87, Path.data, band.width='bw.bcv')
#Edge.KDE <- estimate_edge_KDE(tree.discr, root.taxon=87, Path.data, band.width='bw.ucv')
#saveRDS(Edge.KDE, "data/merged_trees/Edge.KDE_EP_new.rds")
# Edge.KDE <-readRDS("data/merged_trees/Edge.KDE_EP_new.rds")

Edge.KDE$Maps.mean.norm
#Edge.KDE$Band.width

# check if integral over all edges equals 1
integrate_edge_KDE(tree.discr, Edge.KDE$Maps.mean.norm)

#-------- Loess
Edge.KDE$Maps.mean.loess <- Loess_smooting_KDE(tree.discr, Edge.KDE)
Edge.KDE$Maps.mean.loess.norm <-normilize_KDE(tree.discr, Edge.KDE$Maps.mean.loess)
integrate_edge_KDE(tree.discr, Edge.KDE$Maps.mean.loess.norm)
#-----


# analytical posterior for lambda
lambda.post <- posterior_lambda_KDE(tree.list)


#   ____________________________________________________________________________
#   Plotting on trees                                                       ####

# Get POSTERIORS for KDEs and plot them
Edge.KDE$lambda.mean <- make_postPois_KDE(Edge.KDE$Maps.mean.norm, lambda.post, lambda.post.stat='Mean')
Edge.KDE$loess.lambda.mean <- make_postPois_KDE(Edge.KDE$Maps.mean.loess.norm, lambda.post, lambda.post.stat='Mean')
#saveRDS(Edge.KDE, "Edge.KDE_500_EP_1000bins.rds")
#saveRDS(Edge.KDE, "Edge.KDE_500_EP_2000bins.rds")
# PLOT
NHPP.Map.loess <- make_contMap_KDE(tree.discr, Edge.KDE$loess.lambda.mean)
plot.contMap(NHPP.Map.loess, lwd=2, outline=F, legend=F, ftype="off", outline=F)


# Get DERIVATIVES for KDEs and plot them
Edge.KDE$loess.lambda.mean.deriv <-derivative_KDE(tree.discr, Edge.KDE$loess.lambda.mean)
NHPP.Map.loess.deriv <- make_contMap_KDE(tree.discr, 
                                         lapply(Edge.KDE$loess.lambda.mean.deriv, function(x) x^3)
)
plot.contMap(NHPP.Map.loess.deriv, outline=F, lwd=2, ftype="off", outline=F)



#   ____________________________________________________________________________
#   Plotting Edge Profiles                                                  ####

# This gives error due to edge_profiles4plotting() func. Probably earlier version of tibble is needed?

edge.profs$lambda.mean <- edge_profiles4plotting(tree.discr, Edge.KDE$lambda.mean)
edge.profs$loess.lambda.mean <- edge_profiles4plotting(tree.discr, Edge.KDE$loess.lambda.mean)
edge.profs$loess.lambda.mean.deriv <- edge_profiles4plotting(tree.discr, Edge.KDE$loess.lambda.mean.deriv)
edge.profs$loess.lambda.mean.deriv3 <- edge_profiles4plotting(tree.discr, lapply(Edge.KDE$loess.lambda.mean.deriv, function(x) x^3))


#   ____________________________________________________________________________
#   Plotting                                                                ####

ggplot(data = edge.profs$lambda.mean, aes(x = X, y =  Y, group = edge.id, color=Y))+
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
  scale_x_continuous(limits = c(0, 285), breaks=c(seq(0, 280, 50), 280))+
  scale_y_continuous(limits = c(NA, NA))
#geom_vline(xintercept = 0, color='blue')+
#geom_vline(xintercept = 281, color='blue')
                   

ggplot(data = edge.profs$loess.lambda.mean, aes(x = X, y =  Y, group = edge.id, color=Y))+
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
  scale_x_continuous(limits = c(0, 285), breaks=c(seq(0, 280, 50), 280))+
  scale_y_continuous(limits = c(NA, NA))
#geom_vline(xintercept = 0, color='blue')+
#geom_vline(xintercept = 281, color='blue')

ggplot(data = edge.profs$loess.lambda.mean.deriv, aes(x = X, y =  Y, group = edge.id, color=Y))+
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
  scale_x_continuous(limits = c(0, 285), breaks=c(seq(0, 280, 50), 280))+
  scale_y_continuous(limits = c(NA, NA))
#geom_vline(xintercept = 0, color='blue')+
#geom_vline(xintercept = 281, color='blue')

ggplot(data = edge.profs$loess.lambda.mean.deriv3, aes(x = X, y =  Y, group = edge.id, color=Y))+
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
  scale_x_continuous(limits = c(0, 285), breaks=c(seq(0, 280, 50), 280))+
  scale_y_continuous(limits = c(NA, NA))
#geom_vline(xintercept = 0, color='blue')+
#geom_vline(xintercept = 281, color='blue')



#   ____________________________________________________________________________
#   Aligned Plotting                                                        ####


library(gridBase)
library(ggplot2)

#--- Test
# setwd("~/Documents/My_papers/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_BR/BR1")
# tree.test <-lapply(c(1:1), function(x) rds_from_zip(zip.name='female_genitalia.zip', rds.prefix='female_genitalia_', rds.id=x)[[1]])
# closeAllConnections()
# tree.test <-merge_tree_cat(tree.test[[1]])
#---

#png("test_rates.png", width = 774, height = 513, res = NA) 
pdf("test_rates.pdf", width = 7, height = 5)

par(mfrow = c(2,1), mar=c(0,0,0,0), oma=c(0,0,0,0))
plot.contMap(NHPP.Map.loess, lwd=2, outline=F, legend=F, plot=F, ftype="off", mar=c(0.1, 2.85 ,0.1, 1.25))
#plotSimmap(tree.test, lwd=2, pts=F, ftype="off", ylim=c(0,100), mar=c(0.1, 2.85 ,0.1, 1.25))


Pl <- ggplot(data = edge.profs, aes(x = X, y =  Y, group = edge.id, color=Y))+
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
  scale_x_continuous(limits = c(0, 285), breaks=c(seq(0, 280, 50), 280))+
  scale_y_continuous(limits = c(NA, NA))
  #geom_vline(xintercept = 0, color='blue')+
  #geom_vline(xintercept = 281, color='blue')

vp <- viewport(height = unit(.5,"npc"), width=unit(1, "npc"),
               just = c("left",'top'),
               y = 0.5, x = 0)
print(Pl, vp = vp)
dev.off()



