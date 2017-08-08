#####################################################################################################################
############# Plotting phylogenetic alpha diversity within alpine summits ########################################### 
############# Random Draw (RD) Null Model ######################################################################
############# Hannah E. Marx, 25 April 2017 #########################################################################
#####################################################################################################################

############## Tile plots of alpha SES (static null model) per clade
master.ses.alpha <- read.csv(file="output/8_PhyoDiversity/alpha/static/master.ses.static.alpha.csv", row.names=1)
head(master.ses.alpha)
master.ses.alpha$summits
master.ses.alpha$pool <- factor(master.ses.alpha$pool, levels = c( "Regional",  "All Summits", "LGM"))
master.ses.alpha$clade <- factor(master.ses.alpha$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
master.ses.alpha$metric <- factor(master.ses.alpha$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))

plot.pools <- ggplot(master.ses.alpha, aes(y=clade, x=reorder(factor(summits), as.numeric(as.character(ntaxa))), fill=obs.z))
plot.pools <- plot.pools + geom_tile(colour="white", alpha=.75)
plot.pools <- plot.pools + scale_fill_gradient2(low="#44AAAA", high="#AA7744", mid="beige", na.value="white", limits=c(-5, 5)) #"cyan, red"
plot.pools <- plot.pools + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
plot.pools <- plot.pools + scale_size_manual(values=c(dot=1, no_dot=NA), guide="none")
plot.pools <- plot.pools + theme_grey(base_size = 6) + labs(x = "",  y = "") 
plot.pools <- plot.pools + facet_grid(pool ~ metric)#,space="free",scales="free", as.table = F)   
plot.pools <- plot.pools + scale_x_discrete(expand = c(0, 0), labels = c("Brevoort" = "Pointe Brevoort",
                                                                         "Pelvoux"= "Mont Pelvoux",
                                                                         "Choisy" = "Tour Choisy",
                                                                         "Occidentale Ailefroide" = "Occidentale l'Ailefroide",
                                                                         "Plat de la Selle" = "Aiguille du plat de la Selle",
                                                                         "Burlan" = "Pointes de Burlan",
                                                                         "La Meije" = "la Meije",
                                                                         "Barre des Ecrins" = "Barre des Ecrins",
                                                                         "Rouies" = "les Rouies",
                                                                         "Muraillette" = "Muraillette",
                                                                         "Olan"  = "l'Olan",
                                                                         "Sirac" = "le Sirac",
                                                                         "Lauvitel" = "Signal du Lauvitel",
                                                                         "Rateau" = "le Rateau",
                                                                         "Rocher de la Selle" = "Rocher de la Selle",
                                                                         "Persistent" = "LGM",
                                                                         "Summits" = "All Summits")) 
plot.pools <- plot.pools + scale_y_discrete(expand = c(0, 0)) 
plot.pools <- plot.pools + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = 6 * 0.8, angle = -45, hjust = 0, colour = "black"), 
  strip.text.x = element_text(size = 8, face="bold"),
  strip.text.y = element_text(size = 8, face="bold"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES", face="italic"))
plot.pools <- plot.pools + coord_fixed(ratio=1)
plot.pools <- plot.pools + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
plot.pools <- plot.pools + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
plot.pools

pdf(file="figs/Figure3_rd_alpha_SES_tile_pools.pdf")
plot.pools
dev.off()


