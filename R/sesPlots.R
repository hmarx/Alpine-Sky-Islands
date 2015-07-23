
###################################################################################### 
############################### Pool : Ecrins NP ##################################### 
###################################################################################### 

shape.divs.long <- melt(as.matrix(coef(shape.Alpes))) # shape.divs
colnames (shape.divs.long) <- c("summit", "metric", "value")
head(shape.divs.long)

colourCount = length(unique(shape.divs.long$metric))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ggplot(shape.divs.long, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric)) +
  scale_y_log10("metric")  +
  scale_x_discrete("Summit (increasing PD)") +
  #facet_wrap(~metric) + 
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Shape Metrics") +
  theme(plot.title=element_text(size=rel(1.5))) #+
#theme(plot.margin = units(c(0.5,2,0.5,0.5), "cm"))

comm.sesmntd[,"sig"] <- ifelse(comm.sesmntd$mntd.obs.p <= 0.05, "TRUE", "FALSE")
comm.sesmpd.phylonull[,"sig1"] <- ifelse(comm.sesmpd.phylonull$mpd.obs.p <= 0.05, "TRUE", "FALSE")
tmp <- cbind(comm.sesmntd[c(1,6,9)], comm.sesmpd.phylonull[c(1,6,9)], id=rownames(comm.sesmntd))
dispersion.metrics.long <- reshape(tmp, direction = "long", idvar = c("id", "ntaxa"), 
                                   varying = list(c("mntd.obs.z", "mpd.obs.z"), c("sig", "sig1")), 
                                   times = c("mntd.obs.z", "mpd.obs.z"), v.names = c("metric", "value"))
colnames (dispersion.metrics.long) <- c("ntaxa", "ntaxa.1", "summit", "metric", "value", "sig")
head(dispersion.metrics.long)
#11                 mont.pelvoux    17 mntd.obs.z  1.7740514

dispersion.metrics.long <- arrange(dispersion.metrics.long, ntaxa)
dispersion.metrics.long$summit <- factor(dispersion.metrics.long$summit, levels = unique(dispersion.metrics.long$summit))

# Manual levels
#disp_table <- table(dispersion.metrics.long$summit)
#disp_levels <- c(names(disp_table[c(5, 20, 13, 23)]), (names(disp_table)[order(disp_table)][c(-5, -20, -13, -23)]))
#dispersion.metrics.long$summit <- factor(dispersion.metrics.long$summit, as.character(disp_levels))

#pdf(file="output/9_PhyoDiversity/DispersionMetricsEcrinsPool.pdf") 
ggplot(dispersion.metrics.long, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric, fill= factor(sig), shape=factor(metric)), size = 4) +
  scale_y_continuous("obs.z")  +
  scale_x_discrete("Summit (increasing species richness)") +
  scale_fill_manual(values=c("FALSE"= "white", "TRUE" = "grey"), guide="none") +
  scale_color_manual(values=c("chocolate3", "cyan4")) +
  scale_shape_manual(name=" ", values=c("mntd.obs.z" = 21, "mpd.obs.z" =24),  guide="none") +
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Dispersion Metrics\nEcrins Species Pool: null.model = phylogeny.pool") +
  theme(plot.title=element_text(size=rel(1.5))) 
#dev.off()

###################################################################################### 
############################### Reduced Pool : Summits ############################### 
###################################################################################### 


summits.dispersion.metrics.long <- melt(as.matrix(cbind(summits.sesmpd.phylonull[6], summits.sesmntd.phylonull[6])))
colnames (summits.dispersion.metrics.long) <- c("summit", "metric", "value")
head(summits.dispersion.metrics.long)

# Manual levels
disp_table <- table(summits.dispersion.metrics.long$summit)
disp_levels <- c(names(disp_table[18]), (names(disp_table)[order(disp_table)][-18]))
summits.dispersion.metrics.long$summit <- factor(summits.dispersion.metrics.long$summit, levels = disp_levels)

#pdf(file="output/9_PhyoDiversity/DispersionMetrics_summitPool.pdf") 
ggplot(summits.dispersion.metrics.long, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric)) +
  scale_y_continuous("obs.z : null=phylogeny.pool")  +
  scale_x_discrete("Summit") +
  #facet_wrap(~metric) + 
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Dispersion Metrics\nSummit Species Pool: null.model = phylogeny.pool") +
  theme(plot.title=element_text(size=rel(1.5))) 
#dev.off()


###################################################################################### 
############################### Reduced Pool : Persistent ############################ 
###################################################################################### 

persistent.dispersion.metrics.long <- melt(as.matrix(cbind(persistent.sesmpd.phylonull[6], persistent.sesmntd.phylonull[6])))
colnames(persistent.dispersion.metrics.long) <- c("summit", "metric", "value")
head(persistent.dispersion.metrics.long)

# Manual levels
disp_table <- table(persistent.dispersion.metrics.long$summit)
disp_levels <- c(names(disp_table[c(12, 19, 22)]), (names(disp_table)[order(names(disp_table))][c(-12, -19, -22)]))
persistent.dispersion.metrics.long$summit <- factor(persistent.dispersion.metrics.long$summit, levels = disp_levels)

pdf(file="output/9_PhyoDiversity/DispersionMetrics_persistentPool.pdf") 
ggplot(persistent.dispersion.metrics.long, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric)) +
  scale_y_continuous("obs.z : null=phylogeny.pool")  +
  scale_x_discrete("Summit") +
  #facet_wrap(~metric) + 
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Persistent Pool: Dispersion Metrics") +
  theme(plot.title=element_text(size=rel(1.5))) 
dev.off()


###################################################################################### 
############################### Reduced Pool : UnderIce ##############################
###################################################################################### 

UnderIce.dispersion.metrics.long <- melt(as.matrix(cbind(UnderIce.sesmpd.phylonull[6], UnderIce.sesmntd.phylonull[6])))
colnames(UnderIce.dispersion.metrics.long) <- c("summit", "metric", "value")
head(UnderIce.dispersion.metrics.long)

# Manual levels
disp_table <- table(UnderIce.dispersion.metrics.long$summit)
disp_levels <- c(names(disp_table[c(22, 19, 12)]), (names(disp_table)[order(names(disp_table))][c(-22, -12, -19)]))
UnderIce.dispersion.metrics.long$summit <- factor(UnderIce.dispersion.metrics.long$summit, levels = disp_levels)

#pdf(file="output/9_PhyoDiversity/DispersionMetrics_UnderIcePool.pdf") 
ggplot(UnderIce.dispersion.metrics.long, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric)) +
  scale_y_continuous("obs.z : null=phylogeny.pool")  +
  scale_x_discrete("Summit") +
  #facet_wrap(~metric) + 
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("UnderIce Pool: Dispersion Metrics") +
  theme(plot.title=element_text(size=rel(1.5))) 
#dev.off()



###################################################################################### 
########################   Summarize Reduced Source pool effects ########################   
###################################################################################### 

mpd.pools <- list(comm.sesmpd.phylonull, summits.sesmpd.phylonull, persistent.sesmpd.phylonull, UnderIce.sesmpd.phylonull)
mpd.pools.melt <- melt(mpd.pools, id=c(colnames(comm.sesmpd.phylonull)))
mpd.pools.melt <- na.omit(mpd.pools.melt)
mpd.pools.melt$L1 <- factor(mpd.pools.melt$L1)

pdf(file="output/9_PhyoDiversity/mpd.sourcepools.pdf") 
p1 <- ggplot(mpd.pools.melt, aes(mpd.obs.z, colour=L1)) + 
  geom_density() +
  scale_color_discrete(name="Source Pool",
                       breaks=c("1", "2", "3", "4"),
                       labels=c("Ecrins NP", "Summits", "Peristent", "Under Ice")) +
  ggtitle("Distribution of SES mpd for different source pools") 
p1
dev.off()


mntd.pools <- list(comm.sesmntd, summits.sesmntd.phylonull, persistent.sesmntd.phylonull, UnderIce.sesmntd.phylonull)
mntd.pools.melt <- melt(mntd.pools, id=c(colnames(comm.sesmntd)))
mntd.pools.melt <- na.omit(mntd.pools.melt)
mntd.pools.melt$L1 <- factor(mntd.pools.melt$L1)

#pdf(file="output/9_PhyoDiversity/mntd.sourcepools.pdf") 
p1 <- ggplot(mntd.pools.melt, aes(mntd.obs.z, colour=L1)) + 
  geom_density() +
  scale_color_discrete(name="Source Pool",
                       breaks=c("1", "2", "3", "4"),
                       labels=c("Ecrins NP", "Summits", "Peristent", "Under Ice")) +
  ggtitle("Distribution of SES mntd for different source pools") 
p1
#dev.off()




###################################################################################### 
############################### Null Pool: Summits ###################################
###################################################################################### 

##### Summits Source Pool

comm.sesmpd.sourcePoolNull[,"sig"] <- ifelse(comm.sesmpd.sourcePoolNull$mpd.obs.p <= 0.05, "TRUE", "FALSE")
comm.sesmntd.sourcePoolNull[,"sig1"] <- ifelse(comm.sesmntd.sourcePoolNull$mntd.obs.p <= 0.05, "TRUE", "FALSE")
tmp <- cbind(comm.sesmntd.sourcePoolNull[c(1,6,9)], comm.sesmpd.sourcePoolNull[c(1,6,9)], id=rownames(comm.sesmpd.sourcePoolNull))
dispersion.sourceSummit <- reshape(tmp, direction = "long", idvar = c("id", "ntax"), 
                                   varying = list(c("mntd.obs.z", "mpd.obs.z"), c("sig1", "sig")), 
                                   times = c("mntd.obs.z", "mpd.obs.z"), v.names = c("metric", "value"))
colnames (dispersion.sourceSummit) <- c("ntaxa", "ntaxa.1", "summit", "metric", "value", "sig")
head(dispersion.sourceSummit)

dispersion.sourceSummit <- arrange(dispersion.sourceSummit, ntaxa)
dispersion.sourceSummit$summit <- factor(dispersion.sourceSummit$summit, levels = unique(dispersion.sourceSummit$summit))

#pdf(file="output/9_PhyoDiversity/DispersionMetricsEcrinsPool.pdf") 
ggplot(dispersion.sourceSummit, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric, fill= factor(sig), shape=factor(metric)), size = 4) +
  scale_y_continuous("obs.z")  +
  scale_x_discrete("Summit (increasing species richness)") +
  scale_fill_manual(values=c("FALSE"= "white", "TRUE" = "grey"), guide="none") +
  scale_color_manual(values=c("chocolate3", "cyan4")) +
  scale_shape_manual(name=" ", values=c("mntd.obs.z" = 21, "mpd.obs.z" =24),  guide="none") +
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Dispersion Metrics\nLGM Summit Source Pool: null.model = sourcePool") +
  theme(plot.title=element_text(size=rel(1.5))) 
#dev.off()


##### Persistant Source Pool

comm.sesmpd.sourcePoolNull.pers[,"sig"] <- ifelse(comm.sesmpd.sourcePoolNull.pers$mpd.obs.p <= 0.05, "TRUE", "FALSE")
comm.sesmntd.sourcePoolNull.pers[,"sig1"] <- ifelse(comm.sesmntd.sourcePoolNull.pers$mntd.obs.p <= 0.05, "TRUE", "FALSE")
tmp <- cbind(comm.sesmntd.sourcePoolNull.pers[c(1,6,9)], comm.sesmpd.sourcePoolNull.pers[c(1,6,9)], id=rownames(comm.sesmpd.sourcePoolNull.pers))
dispersion.sourcePers <- reshape(tmp, direction = "long", idvar = c("id", "ntax"), 
                                   varying = list(c("mntd.obs.z", "mpd.obs.z"), c("sig1", "sig")), 
                                   times = c("mntd.obs.z", "mpd.obs.z"), v.names = c("metric", "value"))
colnames (dispersion.sourcePers) <- c("ntaxa", "ntaxa.1", "summit", "metric", "value", "sig")
head(dispersion.sourcePers)

dispersion.sourcePers <- arrange(dispersion.sourcePers, ntaxa)
dispersion.sourcePers$summit <- factor(dispersion.sourcePers$summit, levels = unique(dispersion.sourcePers$summit))

# Manual levels
#disp_table <- table(dispersion.metrics.long$summit)
#disp_levels <- c(names(disp_table[c(5, 20, 13, 23)]), (names(disp_table)[order(disp_table)][c(-5, -20, -13, -23)]))
#dispersion.metrics.long$summit <- factor(dispersion.metrics.long$summit, as.character(disp_levels))

#pdf(file="output/9_PhyoDiversity/DispersionMetricsEcrinsPool.pdf") 
ggplot(dispersion.sourcePers, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric, fill= factor(sig), shape=factor(metric)), size = 4) +
  scale_y_continuous("obs.z")  +
  scale_x_discrete("Summit (increasing species richness)") +
  scale_fill_manual(values=c("FALSE"= "white", "TRUE" = "grey"), guide="none") +
  scale_color_manual(values=c("chocolate3", "cyan4")) +
  scale_shape_manual(name=" ", values=c("mntd.obs.z" = 21, "mpd.obs.z" =24),  guide="none") +
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Dispersion Metrics\nLGM Persistent Source Pool: null.model = sourcePool") +
  theme(plot.title=element_text(size=rel(1.5))) 
#dev.off()



########################################################################################################### 
################# Distribution of SES under reduced source pools and  sourcePool Null ##################### 
########################################################################################################### 


colnames(mntd.pools.sourcePool.melt)[1] <- "ntaxa"
levels(mntd.pools.sourcePool.melt$L1) <- c(5,6,7)
all.mntd <- rbind(mntd.pools.melt, mntd.pools.sourcePool.melt)

#pdf(file="output/9_PhyoDiversity/mntd.all.sourcepools.pdf") 
p1 <- ggplot(all.mntd, aes(mntd.obs.z, colour=L1, linetype=L1)) + 
  geom_density() +
  scale_color_discrete(name="Source Pool",
                       breaks=c("1", "2", "3", "4", "5", "6", "7"),
                       labels=c("Ecrins NP_phyNull", "Summits_phyNull", "Peristent_phyNull", "Under Ice_phyNull", 
                                "Summits_sourceNull", "Peristent_sourceNull", "Under Ice_sourceNull")) +
  scale_linetype_manual(name="Source Pool", breaks=c("1", "2", "3", "4", "5", "6", "7"), values=c(1,1,1,1, 5,5,5), 
                        labels=c("Ecrins NP_phyNull", "Summits_phyNull", "Peristent_phyNull", "Under Ice_phyNull", 
                                 "Summits_sourceNull", "Peristent_sourceNull", "Under Ice_sourceNull")) +
  ggtitle("Distribution of SES mntd for different source pools") 
p1
#dev.off()


colnames(mpd.pools.sourcePool.melt)[1] <- "ntaxa"
levels(mpd.pools.sourcePool.melt$L1) <- c(5,6,7)
all.mpd <- rbind(mpd.pools.melt, mpd.pools.sourcePool.melt)

#pdf(file="output/9_PhyoDiversity/mpd.all.sourcepools.pdf") 
p1 <- ggplot(all.mpd, aes(mpd.obs.z, colour=L1, linetype=L1)) + 
  geom_density() +
  scale_color_discrete(name="Source Pool",
                       breaks=c("1", "2", "3", "4", "5", "6", "7"),
                       labels=c("Ecrins NP_phyNull", "Summits_phyNull", "Peristent_phyNull", "Under Ice_phyNull", 
                                "Summits_sourceNull", "Peristent_sourceNull", "Under Ice_sourceNull")) +
  scale_linetype_manual(name="Source Pool", 
                        breaks=c("1", "2", "3", "4", "5", "6", "7"), values=c(1,1,1,1, 5,5,5), 
                        labels=c("Ecrins NP_phyNull", "Summits_phyNull", "Peristent_phyNull", "Under Ice_phyNull", 
                                 "Summits_sourceNull", "Peristent_sourceNull", "Under Ice_sourceNull")) +
  ggtitle("Distribution of SES mpd for different source pools") 
p1
#dev.off()



