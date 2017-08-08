#####################################################################################################################
############# Plotting phylogenetic alpha diversity within alpine summits ########################################### 
############# Both Random Draw (RD) and DAMOCLES Dynamic Null Models ################################################
############# Hannah E. Marx, 25 April 2017 #########################################################################
#####################################################################################################################

######################### Boxplots for each clade SES by species pool ######################### 

############## Master alpha SES divided by clade
master.ses.static <- read.csv(file="output/8_PhyoDiversity/alpha/static/Dryad_master.ses.static.alpha.csv", row.names=1)
head(master.ses.static)
master.ses.static <- cbind(master.ses.static, model = rep("RD", nrow(master.ses.static)))

master.ses.dynamic <- read.csv(file="output/8_PhyoDiversity/alpha/dynamic/Dryad_master.ses.dynamic.alpha.csv", row.names=1)
head(master.ses.dynamic)
master.ses.dynamic.damo <- master.ses.dynamic[master.ses.dynamic$model == "DAMOCLES",]

alpha.div.master <- rbind(master.ses.static, master.ses.dynamic.damo)
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
ggsave("figs/supplemental/alpha.pools.rd.box.pdf", alpha.pools.rd.box, width=8.5, height=8.5)

alpha.pools.damo.box <- ggplot(alpha.div.master[alpha.div.master$model == "DAMOCLES",], aes(x=clade, y=obs.z, fill=pool)) + 
  geom_boxplot() +
  facet_grid(metric ~.) +
  labs(title = "Alpha Diveristy (DAMOCLES)") +
  labs(y = "SES") 
alpha.pools.damo.box
ggsave("figs/supplemental/alpha.pools.damo.box.pdf", alpha.pools.damo.box, width=8.5, height=8.5)


