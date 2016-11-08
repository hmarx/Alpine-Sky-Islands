
## Box and density plots for summit community, Ecrins Source pool
# Density distributions of MNTD and MPD standardized effect sizes for SESDAMOCLES (blue) compared to the SESRD_RD (red) within the for four clades. 

######### Asterales
################################ Ecrins Source Pool ################################ 
## Summit community 
summit.ast.EcrinsPool.summary <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/asterales.summary_table.csv")
summit.ast.EcrinsPool.summary #summit.ast.EcrinsPool.summary

summit.ast.EcrinsPool.com <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/asterales.null_community_data.csv")
head(summit.ast.EcrinsPool.com) #summit.ast.EcrinsPool.com

comDataSESdamocles.ast.mntd <- cbind(summit.ast.EcrinsPool.com %>% 
                                       mutate(RD = (mntd.RD - as.numeric(as.character(summit.ast.EcrinsPool.summary[11,3])))/ as.numeric(as.character(summit.ast.EcrinsPool.summary[12,3]))) %>% 
                                       mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(summit.ast.EcrinsPool.summary[11,3])))/ as.numeric(as.character(summit.ast.EcrinsPool.summary[12,3])))), 
                                     metric = "mntd",
                                     obs.RD = as.numeric(as.character(summit.ast.EcrinsPool.summary[13, 3])), #mntd.obs.RD.z.RD
                                     obs.DAMO = as.numeric(as.character(summit.ast.EcrinsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                     clade = "Asterales")

comDataSESdamocles.ast.mpd <- cbind(summit.ast.EcrinsPool.com %>%  
                                      mutate(RD = ((mpd.RD - as.numeric(as.character(summit.ast.EcrinsPool.summary[16,3])))/ as.numeric(as.character(summit.ast.EcrinsPool.summary[17,3])))) %>% 
                                      mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(summit.ast.EcrinsPool.summary[16,3])))/ as.numeric(as.character(summit.ast.EcrinsPool.summary[17,3])))), 
                                    metric = "mpd",
                                    obs.RD = as.numeric(as.character(summit.ast.EcrinsPool.summary[18, 3])),
                                    obs.DAMO = as.numeric(as.character(summit.ast.EcrinsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                    clade = "Asterales")

comDataSESdamocles.ast.mntd.melt <- melt(comDataSESdamocles.ast.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.ast.mpd.melt <- melt(comDataSESdamocles.ast.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.ast.melt <- rbind(comDataSESdamocles.ast.mntd.melt, comDataSESdamocles.ast.mpd.melt)
tail(comDataSESdamocles.ast.melt)

pdf("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/SES.1000.asterales.pdf")
ggplot(comDataSESdamocles.ast.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Asterales") 
dev.off()



######### Poales
summit.poa.EcrinsPool.summary <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/poales.summary_table.csv")
summit.poa.EcrinsPool.summary #summit.poa.EcrinsPool.summary

summit.ast.EcrinsPool.com <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/poales.null_community_data.csv")
head(summit.ast.EcrinsPool.com)  # #summit.ast.EcrinsPool.com

comDataSESdamocles.poa.mntd <- cbind(summit.ast.EcrinsPool.com %>% 
                                       mutate(RD = (mntd.RD - as.numeric(as.character(summit.poa.EcrinsPool.summary[11,3])))/ as.numeric(as.character(summit.poa.EcrinsPool.summary[12,3]))) %>% 
                                       mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(summit.poa.EcrinsPool.summary[11,3])))/ as.numeric(as.character(summit.poa.EcrinsPool.summary[12,3])))), 
                                     metric = "mntd",
                                     obs.RD = as.numeric(as.character(summit.poa.EcrinsPool.summary[13, 3])), 
                                     obs.DAMO = as.numeric(as.character(summit.poa.EcrinsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                     clade = "Poales")

comDataSESdamocles.poa.mpd <- cbind(summit.ast.EcrinsPool.com %>%  
                                      mutate(RD = ((mpd.RD - as.numeric(as.character(summit.poa.EcrinsPool.summary[16,3])))/ as.numeric(as.character(summit.poa.EcrinsPool.summary[17,3])))) %>% 
                                      mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(summit.poa.EcrinsPool.summary[16,3])))/ as.numeric(as.character(summit.poa.EcrinsPool.summary[17,3])))), 
                                    metric = "mpd",
                                    obs.RD = as.numeric(as.character(summit.poa.EcrinsPool.summary[18, 3])),
                                    obs.DAMO = as.numeric(as.character(summit.poa.EcrinsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                    clade = "Poales")

comDataSESdamocles.poa.mntd.melt <- melt(comDataSESdamocles.poa.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.poa.mpd.melt <- melt(comDataSESdamocles.poa.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.poa.melt <- rbind(comDataSESdamocles.poa.mntd.melt, comDataSESdamocles.poa.mpd.melt)
tail(comDataSESdamocles.poa.melt)

pdf("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/SES.1000.poales.pdf")
ggplot(comDataSESdamocles.poa.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Poales") 
dev.off()



######### Lamiales
summit.lam.EcrinsPool.com <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/Lamiales.null_community_data.csv")
head(summit.lam.EcrinsPool.com)

summit.lam.EcrinsPool.summary <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/Lamiales.summary_table.csv")
summit.lam.EcrinsPool.summary #summit.lam.EcrinsPool.summary

comDataSESdamocles.lam.mntd <- cbind(summit.lam.EcrinsPool.com %>% 
                                        mutate(RD = (mntd.RD - as.numeric(as.character(summit.lam.EcrinsPool.summary[11,3])))/ as.numeric(as.character(summit.lam.EcrinsPool.summary[12,3]))) %>% 
                                        mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(summit.lam.EcrinsPool.summary[11,3])))/ as.numeric(as.character(summit.lam.EcrinsPool.summary[12,3])))), 
                                      metric = "mntd",
                                      obs.RD = as.numeric(as.character(summit.lam.EcrinsPool.summary[13, 3])),
                                     obs.DAMO = as.numeric(as.character(summit.lam.EcrinsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                     clade = "Lamiales")

comDataSESdamocles.lam.mpd <- cbind(summit.lam.EcrinsPool.com %>%  
                                       mutate(RD = ((mpd.RD - as.numeric(as.character(summit.lam.EcrinsPool.summary[16,3])))/ as.numeric(as.character(summit.lam.EcrinsPool.summary[17,3])))) %>% 
                                       mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(summit.lam.EcrinsPool.summary[16,3])))/ as.numeric(as.character(summit.lam.EcrinsPool.summary[17,3])))), 
                                     metric = "mpd",
                                     obs.RD = as.numeric(as.character(summit.lam.EcrinsPool.summary[18, 3])),
                                    obs.DAMO = as.numeric(as.character(summit.lam.EcrinsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                    clade = "Lamiales")

comDataSESdamocles.lam.mntd.melt <- melt(comDataSESdamocles.lam.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.lam.mpd.melt <- melt(comDataSESdamocles.lam.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.lam.melt <- rbind(comDataSESdamocles.lam.mntd.melt, comDataSESdamocles.lam.mpd.melt)
tail(comDataSESdamocles.lam.melt)

pdf("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/SES.1000.Lamiales.pdf")
ggplot(comDataSESdamocles.lam.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Lamiales") 
dev.off()


######### Caryophyllales
summit.cary.EcrinsPool.com <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/caryophyllales.null_community_data.csv")
head(summit.cary.EcrinsPool.com)

summit.cary.EcrinsPool.summary <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/caryophyllales.summary_table.csv")
summit.cary.EcrinsPool.summary #summit.cary.EcrinsPool.summary

comDataSESdamocles.cary.mntd <- cbind(summit.cary.EcrinsPool.com %>% 
                                        mutate(RD = (mntd.RD - as.numeric(as.character(summit.cary.EcrinsPool.summary[11,3])))/ as.numeric(as.character(summit.cary.EcrinsPool.summary[12,3]))) %>% 
                                        mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(summit.cary.EcrinsPool.summary[11,3])))/ as.numeric(as.character(summit.cary.EcrinsPool.summary[12,3])))), 
                                      metric = "mntd",
                                      obs.RD = as.numeric(as.character(summit.cary.EcrinsPool.summary[13, 3])),
                                      obs.DAMO = as.numeric(as.character(summit.cary.EcrinsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                      clade = "Caryophyllales")

comDataSESdamocles.cary.mpd <- cbind(summit.cary.EcrinsPool.com %>%  
                                       mutate(RD = ((mpd.RD - as.numeric(as.character(summit.cary.EcrinsPool.summary[16,3])))/ as.numeric(as.character(summit.cary.EcrinsPool.summary[17,3])))) %>% 
                                       mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(summit.cary.EcrinsPool.summary[16,3])))/ as.numeric(as.character(summit.cary.EcrinsPool.summary[17,3])))), 
                                     metric = "mpd",
                                     obs.RD = as.numeric(as.character(summit.cary.EcrinsPool.summary[18, 3])),
                                     obs.DAMO = as.numeric(as.character(summit.cary.EcrinsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                     clade = "Caryophyllales")

comDataSESdamocles.cary.mntd.melt <- melt(comDataSESdamocles.cary.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.cary.mpd.melt <- melt(comDataSESdamocles.cary.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.cary.melt <- rbind(comDataSESdamocles.cary.mntd.melt, comDataSESdamocles.cary.mpd.melt)
head(comDataSESdamocles.cary.melt)

pdf("output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/SES.1000.Caryophyl.pdf")
ggplot(comDataSESdamocles.cary.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Caryophyllales") 
dev.off()




##### ALL

master.ses <- rbind(comDataSESdamocles.ast.mntd.melt, comDataSESdamocles.ast.mpd.melt,
                    comDataSESdamocles.poa.mntd.melt, comDataSESdamocles.poa.mpd.melt,
                    comDataSESdamocles.lam.mntd.melt, comDataSESdamocles.lam.mpd.melt,
                    comDataSESdamocles.cary.mntd.melt, comDataSESdamocles.cary.mpd.melt)
head(master.ses)
master.ses$value <- -1*master.ses$value
master.ses$obs.RD <- -1*master.ses$obs.RD
master.ses$obs.DAMO <- -1*master.ses$obs.DAMO
summary(master.ses)

ggplot(master.ses, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(clade ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show.legend = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "") 

pdf(file="output/9_PhyoDiversity/Spermatophyta/dynamic/cladesEachSummit_EcrinsPool/summits/SES.all.density.pdf")
ggplot(master.ses, aes(x = value, fill = variable)) + 
  geom_histogram(aes(y=0.5*..density..), alpha=0.7,position='identity',binwidth=0.5) +
  facet_grid(clade ~ metric) +
  geom_vline(aes(xintercept=obs.RD), color = "black", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  #geom_vline(aes(xintercept=obs.DAMO), color = "green", show.legend = F, lty=2) + #
  labs(x = "SES", y = "Density", title = "") +
  scale_fill_discrete(guide = guide_legend(title = "Null Models"))
dev.off()

