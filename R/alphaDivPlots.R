
############## Tile plots of alpha SES divided by clade
master.ses.alpha <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/master.ses.static.alpha.csv", row.names=1)
head(master.ses.alpha)

master.ses.alpha <- master.ses.alpha[!master.ses.alpha$summits == "Under Ice",] #reove this...not interesting
master.ses.alpha$pool <- factor(master.ses.alpha$pool, levels = c( "Ecrins NP",  "Summits", "Persistent LGM"))
master.ses.alpha$pool <- factor(master.ses.alpha$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
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
  axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"), 
  strip.text.x = element_text(size = 8, face="bold"),
  strip.text.y = element_text(size = 8, face="bold"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES", face="italic"))
plot.pools <- plot.pools + coord_fixed(ratio=1)
plot.pools <- plot.pools + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
plot.pools <- plot.pools + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
plot.pools

pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/rd_alpha_SES_tile_pools.pdf")
plot.pools
dev.off()

pdf(file="figs/Figure2_rd_alpha_SES_tile_pools.pdf")
plot.pools
dev.off()


######################### Boxplots for each clade SES by species pool ######################### 

############## Master alpha SES divided by clade
master.ses.static <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/master.ses.static.alpha.csv", row.names=1)
head(master.ses.static)
master.ses.static <- cbind(master.ses.static, model = rep("RD", nrow(master.ses.static)))

master.ses.dynamic <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/dynamic/master.ses.dynamic.alpha.csv", row.names=1)
head(master.ses.dynamic)
master.ses.dynamic.damo <- master.ses.dynamic[master.ses.dynamic$model == "DAMOCLES",]

alpha.div.master <- rbind(master.ses.static[-1], master.ses.dynamic.damo)
head(alpha.div.master)
str(alpha.div.master)
alpha.div.master$pool <- factor(alpha.div.master$pool, levels = c("Ecrins NP", "Summits", "Persistent LGM"))
alpha.div.master$pool <- factor(alpha.div.master$pool, labels = c("Ecrins NP" = "Regional", "Summits" = "All Summits", "Persistent LGM" = "LGM"))
alpha.div.master$clade <- factor(alpha.div.master$clade, levels = c("Spermatophyta", "Asterales", "Poales", "Rosales",  "Lamiales", "Caryophyllales"))
alpha.div.master$metric <- factor(alpha.div.master$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))


alpha.mntd.box <- ggplot(alpha.div.master[alpha.div.master$metric == "MNTD",], aes(x=pool, y=obs.z, fill=clade))
alpha.mntd.box <- alpha.mntd.box + geom_boxplot()
alpha.mntd.box <- alpha.mntd.box + facet_grid(model ~ .)
alpha.mntd.box <- alpha.mntd.box + labs(title = "MNTD")
alpha.mntd.box <- alpha.mntd.box + labs(y = "SES")
alpha.mntd.box
pdf(file="figs/supplemental/alpha.ses.mntd.box.pdf")
alpha.mntd.box
dev.off()


alpha.mpd.box <- ggplot(alpha.div.master[alpha.div.master$metric == "MPD",], aes(x=pool, y=obs.z, fill=clade))
alpha.mpd.box <- alpha.mpd.box + geom_boxplot()
alpha.mpd.box <- alpha.mpd.box + facet_grid(model ~ .)
alpha.mntd.box <- alpha.mntd.box + labs(title = "MPD")
alpha.mntd.box <- alpha.mntd.box + labs(y = "SES")
alpha.mpd.box
pdf(file="figs/supplemental/alpha.ses.mpd.box.pdf")
alpha.mpd.box
dev.off()


alpha.box <- ggplot(alpha.div.master, aes(x=pool, y=obs.z, fill=clade))
alpha.box <- alpha.box + geom_boxplot()
alpha.box <- alpha.box + facet_grid(model ~ metric)
#alpha.box <- alpha.box + labs(title = "Alpha Diveristy")
alpha.box <- alpha.box + labs(y = "SES")
alpha.box <- alpha.box + labs(x = "source pool")
alpha.box
ggsave("figs/supplemental/alpha.ses.box.pdf", alpha.box, width=11, height=8.5)


alpha.pools.rd.box <- ggplot(alpha.div.master[alpha.div.master$model == "RD",], aes(x=clade, y=obs.z, fill=pool)) + 
  geom_boxplot() +
  facet_grid(metric ~.) +
  labs(title = "Alpha Diveristy (RD)") +
  labs(y = "SES") 
alpha.pools.rd.box
ggsave("output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/alpha.pools.rd.box.pdf", alpha.pools.rd.box, width=8.5, height=8.5)
ggsave("figs/supplemental/alpha.pools.rd.box.pdf", alpha.pools.rd.box, width=8.5, height=8.5)


alpha.pools.damo.box <- ggplot(alpha.div.master[alpha.div.master$model == "DAMOCLES",], aes(x=clade, y=obs.z, fill=pool)) + 
  geom_boxplot() +
  facet_grid(metric ~.) +
  labs(title = "Alpha Diveristy (DAMOCLES)") +
  labs(y = "SES") 
alpha.pools.damo.box
ggsave("output/9_PhyoDiversity/Spermatophyta/dynamic/figs/alpha.pools.damo.box.pdf", alpha.pools.damo.box, width=8.5, height=8.5)
ggsave("figs/supplemental/alpha.pools.damo.box.pdf", alpha.pools.damo.box, width=8.5, height=8.5)


