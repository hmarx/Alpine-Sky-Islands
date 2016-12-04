#####################################################################################################################
############# Tile heat plots of DAMOCLES SES output for each summit community  #####################################
############# Regional (Ecrins), All Summits, and LGM species pools ################################################# 
############# Hannah E. Marx, 7 Mar 2016 ############################################################################
#####################################################################################################################


################################ Ecrins Species Pool: summit and LGM communities #################################### 

## Summit community 
summit.ast.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/asterales.summary_table.csv")
summit.ast.EcrinsPool.summary
summit.ast.EcrinsPool.summary_tab <- t(as.data.frame(summit.ast.EcrinsPool.summary[,3], strings.as.factors=F))
colnames(summit.ast.EcrinsPool.summary_tab) <- t(summit.ast.EcrinsPool.summary)[2,]
rownames(summit.ast.EcrinsPool.summary_tab) <- NULL
summit.ast.EcrinsPool.summary_tab <- as.data.frame(cbind(summit = "All Summits", clade = "Asterales", summit.ast.EcrinsPool.summary_tab, pool = "Regional"))

summit.poa.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/poales.summary_table.csv")
summit.poa.EcrinsPool.summary
summit.poa.EcrinsPool.summary_tab <- t(as.data.frame(summit.poa.EcrinsPool.summary[,3], strings.as.factors=F))
colnames(summit.poa.EcrinsPool.summary_tab) <- t(summit.poa.EcrinsPool.summary)[2,]
rownames(summit.poa.EcrinsPool.summary_tab) <- NULL
summit.poa.EcrinsPool.summary_tab <- as.data.frame(cbind(summit = "All Summits", clade = "Poales", summit.poa.EcrinsPool.summary_tab, pool = "Regional"))

summit.lam.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/Lamiales.summary_table.csv")
summit.lam.EcrinsPool.summary
summit.lam.EcrinsPool.summary_tab <- t(as.data.frame(summit.lam.EcrinsPool.summary[,3], strings.as.factors=F))
colnames(summit.lam.EcrinsPool.summary_tab) <- t(summit.lam.EcrinsPool.summary)[2,]
rownames(summit.lam.EcrinsPool.summary_tab) <- NULL
summit.lam.EcrinsPool.summary_tab <- as.data.frame(cbind(summit = "All Summits", clade = "Lamiales", summit.lam.EcrinsPool.summary_tab, pool = "Regional"))

summit.cary.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/caryophyllales.summary_table.csv")
summit.cary.EcrinsPool.summary
summit.cary.EcrinsPool.summary_tab <- t(as.data.frame(summit.cary.EcrinsPool.summary[,3], strings.as.factors=F))
colnames(summit.cary.EcrinsPool.summary_tab) <- t(summit.cary.EcrinsPool.summary)[2,]
rownames(summit.cary.EcrinsPool.summary_tab) <- NULL
summit.cary.EcrinsPool.summary_tab <- as.data.frame(cbind(summit = "All Summits", clade = "Caryophyllales", summit.cary.EcrinsPool.summary_tab, pool = "Regional"))

## Persistent community
persistent.ast.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/asterales.persis.summary_table.csv")
persistent.ast.EcrinsPool.summary
persistent.ast.EcrinsPool.summary_tab <- t(as.data.frame(persistent.ast.EcrinsPool.summary[,3], strings.as.factors=F))
colnames(persistent.ast.EcrinsPool.summary_tab) <- t(persistent.ast.EcrinsPool.summary)[2,]
rownames(persistent.ast.EcrinsPool.summary_tab) <- NULL
persistent.ast.EcrinsPool.summary_tab <- as.data.frame(cbind(summit = "Persistent", clade = "Asterales", persistent.ast.EcrinsPool.summary_tab, pool="Regional"))

persistent.poa.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/poales.persis.summary_table.csv")
persistent.poa.EcrinsPool.summary
persistent.poa.EcrinsPool.summary_tab <- t(as.data.frame(persistent.poa.EcrinsPool.summary[,3], strings.as.factors=F))
colnames(persistent.poa.EcrinsPool.summary_tab) <- t(persistent.poa.EcrinsPool.summary)[2,]
rownames(persistent.poa.EcrinsPool.summary_tab) <- NULL
persistent.poa.EcrinsPool.summary_tab <- as.data.frame(cbind(summit = "Persistent", clade = "Poales", persistent.poa.EcrinsPool.summary_tab, pool="Regional"))

persistent.lam.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/lamiales.persis.summary_table.csv")
persistent.lam.EcrinsPool.summary
persistent.lam.EcrinsPool.summary_tab <- t(as.data.frame(persistent.lam.EcrinsPool.summary[,3], strings.as.factors=F))
colnames(persistent.lam.EcrinsPool.summary_tab) <- t(persistent.lam.EcrinsPool.summary)[2,]
rownames(persistent.lam.EcrinsPool.summary_tab) <- NULL
persistent.lam.EcrinsPool.summary_tab <- as.data.frame(cbind(summit = "Persistent", clade = "Lamiales", persistent.lam.EcrinsPool.summary_tab, pool="Regional"))

persistent.cary.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/caryophyllales.persis.summary_table.csv")
persistent.cary.EcrinsPool.summary
persistent.cary.EcrinsPool.summary_tab <- t(as.data.frame(persistent.cary.EcrinsPool.summary[,3], strings.as.factors=F))
colnames(persistent.cary.EcrinsPool.summary_tab) <- t(persistent.cary.EcrinsPool.summary)[2,]
rownames(persistent.cary.EcrinsPool.summary_tab) <- NULL
persistent.cary.EcrinsPool.summary_tab <- as.data.frame(cbind(summit = "Persistent", clade = "Caryophyllales", persistent.cary.EcrinsPool.summary_tab, pool="Regional"))


################################ Summits Species Pool: within each summit separately ################################ 
source(file="R/chooseClade.R")

DAMOCLES_summits_aster <- data.frame()
for (i in c(2:7, 9:16)){
  file = paste("asterales.summary_table", i, "csv", sep=".")
  summary <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Asterales/", file, sep=""))
  summary_tab <- t(as.data.frame(summary$value, strings.as.factors=F))
  colnames(summary_tab) <- t(as.data.frame(summary[,2]))
  summary_tab.df <- cbind(summit = colnames(alps.damocles.aster$data)[i], clade = "Asterales", summary_tab, pool = "All Summits")
  rownames(summary_tab.df) <- NULL
  DAMOCLES_summits_aster <- rbind(DAMOCLES_summits_aster, summary_tab.df)
  
}
head(DAMOCLES_summits_aster)

DAMOCLES_summits_caryo <- data.frame()
for (i in 2:16){
  file = paste("caryophyllales.summary_table", i, "csv", sep=".")
  summary <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Caryo/", file, sep=""))
  summary_tab <- t(as.data.frame(summary$value, strings.as.factors=F))
  colnames(summary_tab) <- t(as.data.frame(summary[,2]))
  summary_tab.df <- cbind(summit = colnames(alps.damocles.caryophyllales$data)[i], clade = "Caryophyllales", summary_tab, pool = "All Summits")
  rownames(summary_tab.df) <- NULL
  DAMOCLES_summits_caryo <- rbind(DAMOCLES_summits_caryo, summary_tab.df)
  
}
head(DAMOCLES_summits_caryo)


DAMOCLES_summits_lam <- data.frame()
for (i in 2:16){
  file = paste("lamiales.summary_table", i, "csv", sep=".")
  summary <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Lamiales/", file, sep=""))
  summary_tab <- t(as.data.frame(summary[,3], strings.as.factors=F))
  colnames(summary_tab) <- t(summary[,2])
  summary_tab.df <- cbind(summit = colnames(alps.damocles.lamiales$data)[i], clade = "Lamiales", summary_tab, pool = "All Summits")
  rownames(summary_tab.df) <- NULL
  DAMOCLES_summits_lam <- rbind(DAMOCLES_summits_lam, summary_tab.df)
  
}
head(DAMOCLES_summits_lam)


DAMOCLES_summits_poa <- data.frame()
for (i in 2:16){
  file = paste("poales.summary_table", i, "csv", sep=".")
  summary <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Poales/", file, sep=""))
  summary_tab <- t(as.data.frame(summary[,3], strings.as.factors=F))
  colnames(summary_tab) <- t(summary[,2])
  summary_tab.df <- cbind(summit = colnames(alps.damocles.poales$data)[i], clade = "Poales", summary_tab, pool = "All Summits")
  rownames(summary_tab.df) <- NULL
  DAMOCLES_summits_poa <- rbind(DAMOCLES_summits_poa, summary_tab.df)
  
}
head(DAMOCLES_summits_poa)


################################ Persistent Species Pool:  within each summit separately ################################ 

DAMOCLES_persistent_aster <- data.frame()
for (i in c(2:16)){
  file = paste("asterales.persistent.summary_table", i, "csv", sep=".")
  summary <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Asterales/", file, sep=""))
  summary_tab <- t(as.data.frame(summary$value, strings.as.factors=F))
  colnames(summary_tab) <- t(as.data.frame(summary[,2]))
  summary_tab.df <- cbind(summit = colnames(alps.damocles.aster$data)[i], clade = "Asterales", summary_tab, pool = "LGM")
  rownames(summary_tab.df) <- NULL
  DAMOCLES_persistent_aster <- rbind(DAMOCLES_persistent_aster, summary_tab.df)
  
}
head(DAMOCLES_persistent_aster)

DAMOCLES_persistent_caryo <- data.frame()
for (i in c(2:9, 11:12, 14:16)){
  file = paste("caryophyllales.persistent.summary_table", i, "csv", sep=".")
  summary <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Caryo/", file, sep=""))
  summary_tab <- t(as.data.frame(summary$value, strings.as.factors=F))
  colnames(summary_tab) <- t(as.data.frame(summary[,2]))
  summary_tab.df <- cbind(summit = colnames(alps.damocles.caryophyllales$data)[i], clade = "Caryophyllales", summary_tab, pool = "LGM")
  rownames(summary_tab.df) <- NULL
  DAMOCLES_persistent_caryo <- rbind(DAMOCLES_persistent_caryo, summary_tab.df)
  
}
head(DAMOCLES_persistent_caryo)


DAMOCLES_persistent_lam <- data.frame()
for (i in 3:16){
  file = paste("lamiales.persistent.summary_table", i, "csv", sep=".")
  summary <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Lamiales/", file, sep=""))
  summary_tab <- t(as.data.frame(summary[,3], strings.as.factors=F))
  colnames(summary_tab) <- t(summary[,2])
  summary_tab.df <- cbind(summit = colnames(alps.damocles.lamiales$data)[i], clade = "Lamiales", summary_tab, pool = "LGM")
  rownames(summary_tab.df) <- NULL
  DAMOCLES_persistent_lam <- rbind(DAMOCLES_persistent_lam, summary_tab.df)
  
}
head(DAMOCLES_persistent_lam)


DAMOCLES_persistent_poa <- data.frame()
for (i in 2:16){
  file = paste("poales.persistent.summary_table", i, "csv", sep=".")
  summary <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Poales/", file, sep=""))
  summary_tab <- t(as.data.frame(summary[,3], strings.as.factors=F))
  colnames(summary_tab) <- t(summary[,2])
  summary_tab.df <- cbind(summit = colnames(alps.damocles.poales$data)[i], clade = "Poales", summary_tab, pool = "LGM")
  rownames(summary_tab.df) <- NULL
  DAMOCLES_persistent_poa <- rbind(DAMOCLES_persistent_poa, summary_tab.df)
  
}
head(DAMOCLES_persistent_poa)



################################ Combine all output summaries ################################ 
SES_all <- rbind(summit.ast.EcrinsPool.summary_tab, summit.poa.EcrinsPool.summary_tab, summit.lam.EcrinsPool.summary_tab, summit.cary.EcrinsPool.summary_tab,
                 persistent.ast.EcrinsPool.summary_tab, persistent.poa.EcrinsPool.summary_tab, persistent.lam.EcrinsPool.summary_tab, persistent.cary.EcrinsPool.summary_tab,
                 DAMOCLES_summits_aster, DAMOCLES_summits_poa, DAMOCLES_summits_lam, DAMOCLES_summits_caryo,
                 DAMOCLES_persistent_aster, DAMOCLES_persistent_poa, DAMOCLES_persistent_lam, DAMOCLES_persistent_caryo)


SES_RD_mntd <- cbind(SES_all[,c("n.obs", "mntd.obs", "mntd.mean.RD", "mntd.sd.RD", "mntd.obs.rank.RD","mntd.obs.z.RD", "mntd.obs.q.RD", 
                  "runs", "summit")], metric =  rep("mntd", times = nrow(SES_all)), SES_all[,c("clade")],
                  sig=ifelse(!is.na(as.character(SES_all[,"mntd.obs.z.RD"])) & as.numeric(as.character(SES_all[,"mntd.obs.q.RD"])) <= 0.05 | as.numeric(as.character(SES_all[,"mntd.obs.q.RD"])) > 0.95, 1,0), SES_all[,c("pool")],
                  model = rep("RD", times = nrow(SES_all)))

SES_RD_mpd <- cbind(SES_all[,c("n.obs", "mpd.obs", "mpd.mean.RD", "mpd.sd.RD", "mpd.obs.rank.RD","mpd.obs.z.RD", "mpd.obs.q.RD",
                                      "runs", "summit")], metric =  rep("mpd", times = nrow(SES_all)), SES_all[,c("clade")],
                           sig=ifelse(!is.na((SES_all[,"mpd.obs.z.RD"])) & as.numeric(as.character(SES_all[,"mpd.obs.q.RD"])) <= 0.05 | as.numeric(as.character(SES_all[,"mpd.obs.q.RD"])) > 0.95, 1,0), SES_all[,c("pool")],
                    model = rep("RD", times = nrow(SES_all)))

SES_DAMOCLES_mntd <- cbind(SES_all[,c("n.obs", "mntd.obs", "mntd.mean.DAMOCLES", "mntd.sd.DAMOCLES", "mntd.obs.rank.DAMOCLES","mntd.obs.z.DAMOCLES", "mntd.obs.q.DAMOCLES",
                                      "runs", "summit")], metric =  rep("mntd", times = nrow(SES_all)), SES_all[,c("clade")],
                           sig=ifelse(!is.na((SES_all[,"mntd.obs.z.DAMOCLES"])) & as.numeric(as.character(SES_all[,"mntd.obs.q.DAMOCLES"])) <= 0.05 | as.numeric(as.character(SES_all[,"mntd.obs.q.DAMOCLES"])) > 0.95, 1,0), SES_all[,c("pool")], 
                           model = rep("DAMOCLES", times = nrow(SES_all)))

SES_DAMOCLES_mpd <- cbind(SES_all[,c("n.obs", "mpd.obs", "mpd.mean.DAMOCLES", "mpd.sd.DAMOCLES", "mpd.obs.rank.DAMOCLES","mpd.obs.z.DAMOCLES","mpd.obs.q.DAMOCLES",
                               "runs", "summit")], metric =  rep("mpd", times = nrow(SES_all)), SES_all[,c("clade")],
                    sig=ifelse(!is.na((SES_all[,"mpd.obs.z.DAMOCLES"])) & as.numeric(as.character(SES_all[,"mpd.obs.q.DAMOCLES"])) <= 0.05 | as.numeric(as.character(SES_all[,"mpd.obs.q.DAMOCLES"])) > 0.95, 1,0), SES_all[,c("pool")],
                    model = rep("DAMOCLES", times = nrow(SES_all)))

tmpcol <- c("ntaxa", "obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p" , "runs", "summits", "metric", "clade", "sig", "pool", "model" )
colnames(SES_RD_mntd) <- tmpcol
colnames(SES_RD_mpd) <- tmpcol
colnames(SES_DAMOCLES_mntd) <- tmpcol
colnames(SES_DAMOCLES_mpd) <- tmpcol

ses.master.tmp <- rbind(SES_RD_mntd, SES_RD_mpd, SES_DAMOCLES_mntd, SES_DAMOCLES_mpd)
ses.master.tmp$obs.z <- -1*as.numeric(as.character(ses.master.tmp$obs.z))

ses.master <- ses.master.tmp #phylogeny.poolSES, 
str(ses.master)
ses.master$clade <- factor(ses.master$clade, levels = c( "Caryophyllales",  "Lamiales","Rosales", "Poales", "Asterales","Spermatophyta"))
ses.master$metric <- factor(ses.master$metric, levels = c("mntd", "mpd"))
levels(ses.master$metric)
levels(ses.master$metric) <- c("MNTD", "MPD")
head(ses.master)
#write.csv(ses.master, file="output/8_PhyoDiversity/alpha/dynamic/Dryad_master.ses.dynamic.alpha")

## Plot SES mntd/mpd for each summit, each clade, randomization = "phylogeny.pool"
ses.rd.damo <- ggplot(ses.master, aes(y=clade, x=reorder(factor(summits), as.numeric(as.character(ntaxa))), fill=as.numeric(as.character(obs.z))))
ses.rd.damo <- ses.rd.damo + geom_tile(colour="white", alpha=.75)
ses.rd.damo <- ses.rd.damo + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white", limits=c(-5,5))
ses.rd.damo <- ses.rd.damo + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
ses.rd.damo <- ses.rd.damo + scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
ses.rd.damo <- ses.rd.damo + theme_grey(base_size = 6) + labs(x = "",  y = "") 
ses.rd.damo <- ses.rd.damo + facet_grid(metric + model + pool ~ .)#,space="free",scales="free", as.table = F)   
ses.rd.damo <- ses.rd.damo + scale_x_discrete(expand = c(0, 0)) 
ses.rd.damo <- ses.rd.damo + scale_y_discrete(expand = c(0, 0)) 
ses.rd.damo <- ses.rd.damo + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
ses.rd.damo <- ses.rd.damo + coord_fixed(ratio=1)
ses.rd.damo <- ses.rd.damo + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
ses.rd.damo <- ses.rd.damo + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
ses.rd.damo

pdf("output/8_PhyoDiversity/alpha/dynamic/figs/damocles_SES_tile.pdf_all_pools.pdf")
ses.rd.damo
dev.off()

ses.master.damo <- rbind(ses.master[ses.master$model == "DAMOCLES",])
## Plot SES mntd/mpd for each summit, each clade, randomization = "phylogeny.pool"
ses.damo <- ggplot(ses.master.damo, aes(y=clade, x=reorder(factor(summits), as.numeric(as.character(ntaxa))), fill=as.numeric(as.character(obs.z))))
ses.damo <- ses.damo + geom_tile(colour="white", alpha=.75)
ses.damo <- ses.damo + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white", limits=c(-5,5))
ses.damo <- ses.damo + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
ses.damo <- ses.damo + scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")
ses.damo <- ses.damo + theme_grey(base_size = 6) + labs(x = "",  y = "") 
ses.damo <- ses.damo + facet_grid(metric + pool ~ .)#,space="free",scales="free", as.table = F)   
ses.damo <- ses.damo + scale_x_discrete(expand = c(0, 0)) 
ses.damo <- ses.damo + scale_y_discrete(expand = c(0, 0)) 
ses.damo <- ses.damo + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
ses.damo <- ses.damo + coord_fixed(ratio=1)
ses.damo <- ses.damo + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
ses.damo <- ses.damo + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
ses.damo

pdf("output/8_PhyoDiversity/alpha/dynamic/figs/damocles_SES_tile_pools.pdf")
ses.damo
dev.off()

pdf("figs/supplemental/damocles_SES_tile_pools.pdf")
ses.damo
dev.off()


##### Table with rates estimated for each clade | source pool
head(SES_all)
SES_all_table <- SES_all[-c(3, 8:ncol(SES_all)-1)]
head(SES_all_table)

rate_summary <- SES_all_table[, c("clade", "summit", "pool", "mu", "gamma_0", "loglik.obs")]
colnames(rate_summary)[2] <- "communtiy"
write.csv(rate_summary, file="output/8_PhyoDiversity/alpha/dynamic/clade_summit_summary.csv")




