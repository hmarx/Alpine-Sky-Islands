
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




