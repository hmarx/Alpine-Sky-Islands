
############################### Contemporary Source Pool : Ecrins NP ##################################### 
## Random resample from phylogeny pool (==Ecrins NP), equal probability random draw from phylogeny pool
ecrins.sesmpd.phylonull <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/ecrins.sesmpd.phylonull.csv", row.names = 1)
ecrins.sesmntd.phylonull <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/ecrins.sesmntd.phylonull.csv", row.names = 1)

############################### Contemporary Pool: Summits ###########################
## Source pool = summits, equal probability random draw from phylogeny (pruned to summits)
summits.sesmpd.phylonull <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmpd.phylonull.csv", row.names = 1)
summits.sesmntd.phylonull <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmntd.phylonull.csv", row.names = 1)

## Source pool = summits, "abundance" weighted random draw from phylogeny (pruned to summits)
summits.sesmpd.phylonull.abun <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmpd.phylonull.abun.csv", row.names = 1)
summits.sesmntd.phylonull.abun <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmntd.phylonull.abun.csv", row.names = 1)

## source pool = summits (persistent + under ice?), weighted for abundance in summits
summits.sesmpd.sourceSummits <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmpd.sourceSummits.csv", row.names = 1)
summits.sesmntd.sourceSummits <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmntd.sourceSummits.csv", row.names = 1)

############################### Contemporary Pool : Persistent ############################ 
## Source pool = Persistent, equal probability random draw from phylogeny (pruned to persistent species)
persistent.sesmpd.phylonull <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmpd.phylonull.csv", row.names = 1)
persistent.sesmntd.phylonull <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmntd.phylonull.csv", row.names = 1)

## Source pool = Persistent, "abundance" weighted draw from phylogeny (pruned to persistent species)
persistent.sesmpd.phylonull.abun <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmpd.phylonull.abun.csv", row.names = 1)
persistent.sesmntd.phylonull.abun <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmntd.phylonull.abun.csv", row.names = 1)

## source pool = persistent above glacier through LGM, weighted for abundance in persistent pool
persistent.sesmpd.sourcePersis <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmpd.sourcePersis.csv", row.names = 1)
persistent.sesmntd.sourcePersis <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmntd.sourcePersis.csv", row.names = 1)

## source pool = summits (persistent + under ice?), weighted for abundance in summits
persistent.sesmpd.sourceSummits <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmpd.sourceSummits.csv", row.names = 1)
persistent.sesmntd.sourceSummits <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmntd.sourceSummits.csv", row.names = 1)



############################### Five Clades Separately ###############################
############################### Contemporary Source Pool : Ecrins NP #################

phylogeny.poolSES <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/phylogeny.pool.SES.csv", row.names=1)
head(phylogeny.poolSES)
str(phylogeny.poolSES)

## Plot SES mntd/mpd for each summit, each clade, randomization = "phylogeny.pool"
ses.plot.function <- ggplot(phylogeny.poolSES, aes(y=clade, x=reorder(factor(summits), as.numeric(as.character(ntaxa))), fill=obs.z))#, col=c("magenta1", "green3")) y=reorder(factor(island), as.numeric(as.character(Area.m2)))
ses.plot.function <- ses.plot.function + geom_tile(colour="white", alpha=.75)
ses.plot.function <- ses.plot.function + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white")
ses.plot.function <- ses.plot.function + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
ses.plot.function <- ses.plot.function + scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
ses.plot.function <- ses.plot.function + theme_grey(base_size = 6) + labs(x = "",  y = "") 
ses.plot.function <- ses.plot.function + facet_grid(metric ~ .)#,space="free",scales="free", as.table = F)   
ses.plot.function <- ses.plot.function + scale_x_discrete(expand = c(0, 0)) 
ses.plot.function <- ses.plot.function + scale_y_discrete(expand = c(0, 0)) 
ses.plot.function <- ses.plot.function + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
ses.plot.function <- ses.plot.function + coord_fixed(ratio=1)
ses.plot.function <- ses.plot.function + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
ses.plot.function <- ses.plot.function + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/SES_RD_phylopool.pdf")
ses.plot.function
dev.off()

pdf(file="figs/SES_RD_phylopool.pdf")
ses.plot.function
dev.off()

############################### Five Clades Separately ###############################
############################### Contemporary Source Pool : Summits #################

summit.poolSES <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/summit.pool.SES.csv")
head(summit.poolSES)

## Plot SES mntd/mpd for each summit, each clade, randomization = "phylogeny.pool"
ses.summit.plot.function <- ggplot(summit.poolSES, aes(y=clade, x=reorder(factor(summits), as.numeric(as.character(ntaxa))), fill=obs.z))#, col=c("magenta1", "green3")) y=reorder(factor(island), as.numeric(as.character(Area.m2)))
ses.summit.plot.function <- ses.summit.plot.function + geom_tile(colour="white", alpha=.75)
ses.summit.plot.function <- ses.summit.plot.function + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white")
ses.summit.plot.function <- ses.summit.plot.function + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
ses.summit.plot.function <- ses.summit.plot.function + scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
ses.summit.plot.function <- ses.summit.plot.function + theme_grey(base_size = 6) + labs(x = "",  y = "") 
ses.summit.plot.function <- ses.summit.plot.function + facet_grid(metric ~ .)#,space="free",scales="free", as.table = F)   
ses.summit.plot.function <- ses.summit.plot.function + scale_x_discrete(expand = c(0, 0)) 
ses.summit.plot.function <- ses.summit.plot.function + scale_y_discrete(expand = c(0, 0)) 
ses.summit.plot.function <- ses.summit.plot.function + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
ses.summit.plot.function <- ses.summit.plot.function + coord_fixed(ratio=1)
ses.summit.plot.function <- ses.summit.plot.function + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
ses.summit.plot.function <- ses.summit.plot.function + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")

pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/SES_RD_summitpool.pdf")
ses.summit.plot.function
dev.off()

pdf(file="figs/SES_RD_summitpool.pdf")
ses.summit.plot.function
dev.off()



############################### Five Clades Separately ###############################
############################### Contemporary Source Pool : persistent #################

persistent.poolSES <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/persistent.pool.SES.csv")
head(persistent.poolSES)

## Plot SES mntd/mpd for each persistent, each clade, randomization = "phylogeny.pool"
ses.persistent.plot.function <- ggplot(persistent.poolSES, aes(y=clade, x=reorder(factor(summits), as.numeric(as.character(ntaxa))), fill=obs.z))#, col=c("magenta1", "green3")) y=reorder(factor(island), as.numeric(as.character(Area.m2)))
ses.persistent.plot.function <- ses.persistent.plot.function + geom_tile(colour="white", alpha=.75)
ses.persistent.plot.function <- ses.persistent.plot.function + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white")
ses.persistent.plot.function <- ses.persistent.plot.function + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
ses.persistent.plot.function <- ses.persistent.plot.function + scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
ses.persistent.plot.function <- ses.persistent.plot.function + theme_grey(base_size = 6) + labs(x = "",  y = "") 
ses.persistent.plot.function <- ses.persistent.plot.function + facet_grid(metric ~ .)#,space="free",scales="free", as.table = F)   
ses.persistent.plot.function <- ses.persistent.plot.function + scale_x_discrete(expand = c(0, 0)) 
ses.persistent.plot.function <- ses.persistent.plot.function + scale_y_discrete(expand = c(0, 0)) 
ses.persistent.plot.function <- ses.persistent.plot.function + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
ses.persistent.plot.function <- ses.persistent.plot.function + coord_fixed(ratio=1)
ses.persistent.plot.function <- ses.persistent.plot.function + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
ses.persistent.plot.function <- ses.persistent.plot.function + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
ses.persistent.plot.function

pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/SES_RD_persistentpool.pdf")
ses.persistent.plot.function
dev.off()

pdf(file="figs/SES_RD_persistentpool.pdf")
ses.persistent.plot.function
dev.off()


############## Master alpha SES divided by clade
master.ses.alpha <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/master.ses.static.alpha.csv", row.names=1)
head(master.ses.alpha)

master.ses.alpha <- master.ses.alpha[!master.ses.alpha$summits == "Under Ice",] #reove this...not interesting
master.ses.alpha$pool <- factor(master.ses.alpha$pool, levels = c( "Ecrins NP",  "Summits", "Persistent LGM"))

master.ses.alpha.mntd <- master.ses.alpha[master.ses.alpha$metric == "mntd",]
master.ses.alpha.mpd <- master.ses.alpha[master.ses.alpha$metric == "mpd",]

master.ses.alpha.mntd$clade <- factor(master.ses.alpha.mntd$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
plot.mntd.pools <- ggplot(master.ses.alpha.mntd, aes(y=clade, x=reorder(factor(summits), as.numeric(as.character(ntaxa))), fill=obs.z))
plot.mntd.pools <- plot.mntd.pools + geom_tile(colour="white", alpha=.75)
plot.mntd.pools <- plot.mntd.pools + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white", limits=c(-5, 5))
plot.mntd.pools <- plot.mntd.pools + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
plot.mntd.pools <- plot.mntd.pools + scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
plot.mntd.pools <- plot.mntd.pools + theme_grey(base_size = 6) + labs(x = "",  y = "") 
plot.mntd.pools <- plot.mntd.pools + facet_grid(pool ~ .)#,space="free",scales="free", as.table = F)   
plot.mntd.pools <- plot.mntd.pools + scale_x_discrete(expand = c(0, 0)) 
plot.mntd.pools <- plot.mntd.pools + scale_y_discrete(expand = c(0, 0)) 
plot.mntd.pools <- plot.mntd.pools + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
plot.mntd.pools <- plot.mntd.pools + coord_fixed(ratio=1)
plot.mntd.pools <- plot.mntd.pools + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
plot.mntd.pools <- plot.mntd.pools + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/rd_SES_mntd_tile_pools.pdf")
plot.mntd.pools
dev.off()

pdf(file="figs/rd_SES_mntd_tile_pools.pdf")
plot.mntd.pools
dev.off()

master.ses.alpha.mpd$clade <- factor(master.ses.alpha.mpd$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
plot.mpd.pools <- ggplot(master.ses.alpha.mpd, aes(y=clade, x=reorder(factor(summits), as.numeric(as.character(ntaxa))), fill=obs.z))
plot.mpd.pools <- plot.mpd.pools + geom_tile(colour="white", alpha=.75)
plot.mpd.pools <- plot.mpd.pools + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white", limits=c(-5, 5))
plot.mpd.pools <- plot.mpd.pools + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
plot.mpd.pools <- plot.mpd.pools + scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
plot.mpd.pools <- plot.mpd.pools + theme_grey(base_size = 6) + labs(x = "",  y = "") 
plot.mpd.pools <- plot.mpd.pools + facet_grid(pool ~ .)#,space="free",scales="free", as.table = F)   
plot.mpd.pools <- plot.mpd.pools + scale_x_discrete(expand = c(0, 0)) 
plot.mpd.pools <- plot.mpd.pools + scale_y_discrete(expand = c(0, 0)) 
plot.mpd.pools <- plot.mpd.pools + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
plot.mpd.pools <- plot.mpd.pools + coord_fixed(ratio=1)
plot.mpd.pools <- plot.mpd.pools + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
plot.mpd.pools <- plot.mpd.pools + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/rd_SES_mpd_tile_pools.pdf")
plot.mpd.pools
dev.off()

pdf(file="figs/rd_SES_mpd_tile_pools.pdf")
plot.mpd.pools
dev.off()


