

################ Distrubution under different null models 
plotDistributionSES <- function(outputSES, breaks, labels, mainTitle, values, colors, max.y){
  pools <- outputSES
  pools.melt <- melt(pools, id=c(colnames(pools[[1]])))
  pools.melt <- na.omit(pools.melt)
  pools.melt$L1 <- factor(pools.melt$L1)
  ggplot(pools.melt, aes_q(x = as.name(names(pools.melt)[6]), colour = as.name(names(pools.melt)[9]), 
                           linetype=as.name(names(pools.melt)[9]))) + 
    geom_density() +
    scale_color_manual(name="Source Pool", breaks=c(as.character(breaks)), labels=c(as.character(labels)), values = colors) +
    scale_linetype_manual(name="Source Pool", values=values,
                          breaks=c(as.character(breaks)), labels=c(as.character(labels))) +
    scale_x_continuous(limits=c(-5,5)) +
    scale_y_continuous(limits=c(0,max.y)) +
    ggtitle(label = mainTitle)
  
}


############# Metrics and significance across communities 
plotSESdispersion <- function(mntd, mpd, mainTitle){
  
  mntd[,"sig"] <- ifelse(mntd$mntd.obs.p <= 0.05, "TRUE", "FALSE")
  mpd[,"sig1"] <- ifelse(mpd$mpd.obs.p <= 0.05, "TRUE", "FALSE")
  tmp <- cbind(mntd[c(1,6,9)], mpd[c(1,6,9)], id=rownames(mntd))
  dispersion.metrics.long <- reshape(tmp, direction = "long", idvar = c("id", "ntaxa"), 
                                     varying = list(c("mntd.obs.z", "mpd.obs.z"), c("sig", "sig1")), 
                                     times = c("mntd.obs.z", "mpd.obs.z"), v.names = c("metric", "value"))
  colnames (dispersion.metrics.long) <- c("ntaxa", "ntaxa.1", "summit", "metric", "value", "sig")
  head(dispersion.metrics.long)
  #11                 mont.pelvoux    17 mntd.obs.z  1.7740514
  
  dispersion.metrics.long <- arrange(dispersion.metrics.long, ntaxa)
  dispersion.metrics.long$summit <- factor(dispersion.metrics.long$summit, levels = unique(dispersion.metrics.long$summit))
  
  ggplot(dispersion.metrics.long, aes(x=summit, y=value)) +
    geom_point(aes(colour = metric, fill= factor(sig), shape=factor(metric)), size = 4) +
    scale_y_continuous("obs.z")  +
    scale_x_discrete("Summit (increasing species richness)") +
    scale_fill_manual(values=c("FALSE"= "white", "TRUE" = "lightgrey"), guide="none") +
    scale_color_manual(values=c("chocolate3", "cyan4")) +
    scale_shape_manual(name=" ", values=c("mntd.obs.z" = 21, "mpd.obs.z" =21),  guide="none") +
    theme_bw() +
    theme(legend.position="top") +
    theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
    ggtitle(mainTitle) +
    theme(plot.title=element_text(size=rel(1.5))) 
}


# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

