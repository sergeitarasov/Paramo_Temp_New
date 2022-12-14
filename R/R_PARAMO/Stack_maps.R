
#######################
# Stacking archivi
# sdsds
# read descritized trees from archive
setwd("~/Documents/Recon-Anc_Anat/Supplementary_materials/STEP_5/Discr_maps")

level2
#rm(F)
F1<-F.obj
# mouthparts
cc<-sub("CHAR:", "C",
       get_descendants_chars(ONT, annotations="manual", level2[3] )  )

ntrees<-100
tr<-vector("list", ntrees)

setwd("~/Documents/Recon-Anc_Anat/Supplementary_materials/STEP_5/Discr_maps")
i=1
for (i in 1:ntrees){
  
  fl<-paste0(cc, "_", i, ".rds")  

  stack.L<-vector("list", length(fl))
  j=1
  for (j in 1:length(fl)){
    
    #con<-unz("CHAR-7.zip", filename="CHAR-7_1.rds")
    print(paste0("Reading ", paste0(cc[j], ".zip"), " and ", fl[j]))
    con<-unz(paste0(cc[j], ".zip"), filename=paste0(dirW, fl[j]) )
    con2 <- gzcon(con)
    stack.L[[j]] <- readRDS(con2)
    close(con)
  }
  
  tr[[i]]<- stack_stm(stack.L)
}


pdf(file='test.pdf')
svg(filename='test.svg')
png(filename='test.png')
jpeg(filename='test.jpg')
plot(tr[[1]], pts=F,ftype="off")
dev.off()

library("png", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("jpeg", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
img <- readPNG(system.file("img", "test.png", package="png"))
img <- readJPEG(system.file("img", "test.jpg", package="jpeg"))
img<-readJPEG( "test.jpg")
rasterImage(img)
plot(img)

library("RColorBrewer", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
z<-tr[[1]]
lapply(z$maps, names) %>% unlist %>% unique->states
length(states)

hm.palette <- colorRampPalette(brewer.pal(9, 'Spectral') %>% rev, space='Lab')
hm.palette <- colorRampPalette(brewer.pal(3, 'Blues'), space='Lab')
hm.palette <- colorRampPalette(brewer.pal(3, 'Accent') %>% rev, space='Lab')
hm.palette <- colorRampPalette(brewer.pal(9, 'Set1'), space='Lab')
color<-hm.palette(length(states))

plotSimmap(z, setNames(color, states),  lwd=3, pts=F,ftype="off")
plot(tr[[1]])

#####################################################################
setwd("~/Documents/Recon-Anc_Anat/PARAMO/STEP_5/Discr_maps/maps_pics")

z<-L2.maps$head[[1]]
z<-L2.maps$wings[[1]]
z<-L2.maps$legs[[1]]

z<-L3.maps$whole_organism[[1]]

lapply(z$maps, names) %>% unlist %>% unique->states
hm.palette <- colorRampPalette(brewer.pal(9, 'Set1'), space='Lab')
color<-hm.palette(length(states))

plotSimmap(z, setNames(color, states),  lwd=3, pts=F,ftype="off")

svg(filename='head.svg')
plotSimmap(z, setNames(color, states),  lwd=3, pts=F,ftype="off")
dev.off()

svg(filename='wings.svg')
plotSimmap(z, setNames(color, states),  lwd=3, pts=F,ftype="off")
dev.off()

svg(filename='legs.svg')
plotSimmap(z, setNames(color, states),  lwd=3, pts=F,ftype="off")
dev.off()



lapply(z$maps, names) %>% unlist %>% unique->states
hm.palette <- colorRampPalette(brewer.pal(8, 'Set1'), space='Lab')
color<-hm.palette(length(states))

svg(filename='all_body.svg')
plotSimmap(z, setNames(color, states),  lwd=3, pts=F,ftype="off")
dev.off()

# intial characters
setwd("~/Documents/Recon-Anc_Anat/PARAMO/STEP_5/Discr_maps")
hh<-phytools::read.simmap(file="C3-2.stmr", format="phylip")
plot(hh[[4]])


svg(filename='inti1.svg')
plotSimmap(hh[[5]],  lwd=3, pts=F,ftype="off")
dev.off()

svg(filename='inti2.svg')
plotSimmap(hh[[2]],  lwd=3, pts=F,ftype="off")
dev.off()

svg(filename='inti3.svg')
plotSimmap(hh[[3]],  lwd=3, pts=F,ftype="off")
dev.off()


svg(filename='inti4.svg')
plotSimmap(hh[[1]],  lwd=3, pts=F,ftype="off")
dev.off()
