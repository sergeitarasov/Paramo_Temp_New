

#' Title Make color pallette for Image plotting
#'
#' @param Stat vecor where values is stat, and names are BRs labels
#' @param palette color pallete palette<- rainbow(100, start = 0, end = 0.7) %>% rev
#'
#' @return
#' @export
#'
#' @examples
make_colors <- function(Stat, palette){
  # Colors
  ncols <- length(palette)
  lims <- c(min(Stat), max(Stat))
  
  #cols <- rainbow(ncols, start = 0, end = 0.7) %>% rev
  cols <-palette
  names(cols) <- 1:ncols
  
  # remap, scale colors
  odr <- (1+((ncols - 1)/(lims[2] - lims[1]))*(Stat-lims[1])) %>% round(.,0)
  cols <-cols[odr]
  
  cols.maps <- setNames(cols, names(Stat))
  
  return(cols.maps)
}


make_colors_relative_scale <- function(Stat, palette, lims){
  # Colors
  ncols <- length(palette)
  #lims <- c(min(Stat), max(Stat))
  
  #cols <- rainbow(ncols, start = 0, end = 0.7) %>% rev
  cols <-palette
  names(cols) <- 1:ncols
  
  # remap, scale colors
  odr <- (1+((ncols - 1)/(lims[2] - lims[1]))*(Stat-lims[1])) %>% round(.,0)
  cols <-cols[odr]
  
  cols.maps <- setNames(cols, names(Stat))
  
  return(cols.maps)
}

#' Title Assign colors to Picture layers
#'
#' @param petal Picture obj
#' @param Layers vecotr where values is layer id, and name is BR label
#' @param cols.maps vecor where values is color id, and name is BR label
#'
#' @return
#' @export
#'
#' @examples
make_pic<-function(petal, Layers, cols.maps){
  
  i <- 2
  for (i in 1:length(Layers)){
    petal@paths[Layers[i]]$path@rgb <- cols.maps[names(cols.maps)==names(Layers[i])]
  }
  
  return(petal)
}


get_vector_ids_per_term <- function(term ='HAO:0000349', ONT, GL){
  
  all.ids <-names(GL)
  
  des <- get_descendants(ONT, roots=term, exclude_roots = FALSE)
  pp <- all.ids[(all.ids %in% des)]
  selected.ids <- lapply(GL[pp], function(x) x$layer) %>% unlist %>% as.numeric()
  selected.ids <-selected.ids[!is.na(selected.ids)]
  return(selected.ids)
}


##########################################
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  #dev.new(width=1.75, height=5)
  #par(mar=c(bottom=3, left=3, top=19, right=35))
  #plot(0:4,  type="n", axes=FALSE, xlab = "", ylab = "", bty="n")
  plot(c(0,1), c(min,max),  type="n",  bty='n', xlab='',  ylab='', main=title, xaxt='n',  yaxt='n')
  #plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }	
}
#########