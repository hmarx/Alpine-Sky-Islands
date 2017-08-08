#####################################################################################################################
############# Box and Density Plots of observed SES compared to expected under DAMOCLES null model ##################
############# Regional (Ecrins) Species pool, Summits & LGM communities #############################################
############# Hannah E. Marx, 7 Mar 2016 ############################################################################
#####################################################################################################################


################################ Ecrins Species Pool: summit community ################################ 

######### Asterales ######### 
summit.ast.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/asterales.summary_table.csv")
summit.ast.EcrinsPool.com <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/asterales.null_community_data.csv")
head(summit.ast.EcrinsPool.summary)
head(summit.ast.EcrinsPool.com)

comDataSESdamocles.ast.mntd <- cbind(summit.ast.EcrinsPool.com %>% 
                                       mutate(RD = (mntd.RD - as.numeric(as.character(summit.ast.EcrinsPool.summary[11,3])))/ as.numeric(as.character(summit.ast.EcrinsPool.summary[12,3]))) %>%  #mntd.mean.RD
                                       mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(summit.ast.EcrinsPool.summary[11,3])))/ as.numeric(as.character(summit.ast.EcrinsPool.summary[12,3])))), #mntd.sd.RD
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

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/SES.1000.asterales.pdf")
ggplot(comDataSESdamocles.ast.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Asterales") 
dev.off()



######### Poales ######### 
summit.poa.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/poales.summary_table.csv")
summit.ast.EcrinsPool.com <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/poales.null_community_data.csv")

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

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/SES.1000.poales.pdf")
ggplot(comDataSESdamocles.poa.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Poales") 
dev.off()



######### Lamiales ############
summit.lam.EcrinsPool.com <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/Lamiales.null_community_data.csv")
summit.lam.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/Lamiales.summary_table.csv")

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

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/SES.1000.Lamiales.pdf")
ggplot(comDataSESdamocles.lam.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Lamiales") 
dev.off()


######### Caryophyllales #############
summit.cary.EcrinsPool.com <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/caryophyllales.null_community_data.csv")
summit.cary.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/caryophyllales.summary_table.csv")

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

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/SES.1000.Caryophyl.pdf")
ggplot(comDataSESdamocles.cary.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Caryophyllales") 
dev.off()


##### ALL clades together ######

master.ses <- rbind(comDataSESdamocles.ast.mntd.melt, comDataSESdamocles.ast.mpd.melt,
                    comDataSESdamocles.poa.mntd.melt, comDataSESdamocles.poa.mpd.melt,
                    comDataSESdamocles.lam.mntd.melt, comDataSESdamocles.lam.mpd.melt,
                    comDataSESdamocles.cary.mntd.melt, comDataSESdamocles.cary.mpd.melt)
head(master.ses)
master.ses$value <- -1*master.ses$value
master.ses$obs.RD <- -1*master.ses$obs.RD
master.ses$obs.DAMO <- -1*master.ses$obs.DAMO
summary(master.ses)
master.ses$metric <- factor(master.ses$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))

ggplot(master.ses, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(clade ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show.legend = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "") 

ggplot(filter(master.ses, clade !="Poales"), aes(x = value, fill = variable)) + 
  geom_histogram(aes(y=0.5*..density.., alpha=variable),position='identity',binwidth=0.2) +
  facet_grid(clade ~ metric) +
  geom_vline(aes(xintercept=obs.RD), color = "black", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  #geom_vline(aes(xintercept=obs.DAMO), color = "green", show.legend = F, lty=2) + #
  labs(x = "SES", y = "Density", title = "") +
  scale_fill_manual(guide = guide_legend(title = "Null Models"), values = c("RD" = "#4477AA", "DAMOCLES" = "#AAAA44")) + 
  ggtitle("Ecrins Species Pool, Summit Community") + 
  scale_alpha_discrete(range = c("RD" = 1, "DAMOCLES" = .7), guide=FALSE) +
  scale_x_continuous(limits = c(-5,5)) + 
  scale_y_continuous(limits = c(0, .25)) +
  theme(strip.text.x = element_text(size = 12, face="bold"))
ggsave(filename="output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/SES.all.density.summits.pdf")
ggsave(filename="figs/Figure4a_SES.all.density.summits.pdf")


pdf(file="output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/AllSummits/SES.all.density.green.pdf")
ggplot(master.ses, aes(x = value, fill = variable)) + 
  geom_histogram(aes(y=0.5*..density..), alpha=0.7,position='identity',binwidth=0.5) +
  facet_grid(clade ~ metric) +
  geom_vline(aes(xintercept=obs.RD), color = "black", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  geom_vline(aes(xintercept=obs.DAMO), color = "green", show.legend = F, lty=2) + #
  labs(x = "SES", y = "Density", title = "") +
  scale_fill_discrete(guide = guide_legend(title = "Null Models")) + 
  ggtitle("Ecrins Species Pool, Summit Community") + 
  scale_x_continuous(limits = c(-5,5)) + 
  scale_y_continuous(limits = c(0, .25))
dev.off()


################################ Ecrins Species Pool: LGM community ################################ 

######### Asterales #############
persist.ast.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/asterales.persis.summary_table.csv")
persist.ast.EcrinsPool.com <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/asterales.persis.null_community_data.csv")

comDataSESdamocles.ast.mntd <- cbind(persist.ast.EcrinsPool.com %>% 
                                       mutate(RD = (mntd.RD - as.numeric(as.character(persist.ast.EcrinsPool.summary[11,3])))/ as.numeric(as.character(persist.ast.EcrinsPool.summary[12,3]))) %>% 
                                       mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(persist.ast.EcrinsPool.summary[11,3])))/ as.numeric(as.character(persist.ast.EcrinsPool.summary[12,3])))), 
                                     metric = "mntd",
                                     obs.RD = as.numeric(as.character(persist.ast.EcrinsPool.summary[13, 3])), #mntd.obs.RD.z.RD
                                     obs.DAMO = as.numeric(as.character(persist.ast.EcrinsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                     clade = "Asterales")

comDataSESdamocles.ast.mpd <- cbind(persist.ast.EcrinsPool.com %>%  
                                      mutate(RD = ((mpd.RD - as.numeric(as.character(persist.ast.EcrinsPool.summary[16,3])))/ as.numeric(as.character(persist.ast.EcrinsPool.summary[17,3])))) %>% 
                                      mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(persist.ast.EcrinsPool.summary[16,3])))/ as.numeric(as.character(persist.ast.EcrinsPool.summary[17,3])))), 
                                    metric = "mpd",
                                    obs.RD = as.numeric(as.character(persist.ast.EcrinsPool.summary[18, 3])),
                                    obs.DAMO = as.numeric(as.character(persist.ast.EcrinsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                    clade = "Asterales")

comDataSESdamocles.ast.mntd.melt <- melt(comDataSESdamocles.ast.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.ast.mpd.melt <- melt(comDataSESdamocles.ast.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.ast.melt <- rbind(comDataSESdamocles.ast.mntd.melt, comDataSESdamocles.ast.mpd.melt)
tail(comDataSESdamocles.ast.melt)

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/SES.1000.asterales.pdf")
ggplot(comDataSESdamocles.ast.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Asterales") 
dev.off()



######### Poales #############
persist.poa.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/poales.persis.summary_table.csv")
persist.ast.EcrinsPool.com <- read.csv("output/output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/poales.persis.null_community_data.csv")

comDataSESdamocles.poa.mntd <- cbind(persist.ast.EcrinsPool.com %>% 
                                       mutate(RD = (mntd.RD - as.numeric(as.character(persist.poa.EcrinsPool.summary[11,3])))/ as.numeric(as.character(persist.poa.EcrinsPool.summary[12,3]))) %>% 
                                       mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(persist.poa.EcrinsPool.summary[11,3])))/ as.numeric(as.character(persist.poa.EcrinsPool.summary[12,3])))), 
                                     metric = "mntd",
                                     obs.RD = as.numeric(as.character(persist.poa.EcrinsPool.summary[13, 3])), 
                                     obs.DAMO = as.numeric(as.character(persist.poa.EcrinsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                     clade = "Poales")

comDataSESdamocles.poa.mpd <- cbind(persist.ast.EcrinsPool.com %>%  
                                      mutate(RD = ((mpd.RD - as.numeric(as.character(persist.poa.EcrinsPool.summary[16,3])))/ as.numeric(as.character(persist.poa.EcrinsPool.summary[17,3])))) %>% 
                                      mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(persist.poa.EcrinsPool.summary[16,3])))/ as.numeric(as.character(persist.poa.EcrinsPool.summary[17,3])))), 
                                    metric = "mpd",
                                    obs.RD = as.numeric(as.character(persist.poa.EcrinsPool.summary[18, 3])),
                                    obs.DAMO = as.numeric(as.character(persist.poa.EcrinsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                    clade = "Poales")

comDataSESdamocles.poa.mntd.melt <- melt(comDataSESdamocles.poa.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.poa.mpd.melt <- melt(comDataSESdamocles.poa.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.poa.melt <- rbind(comDataSESdamocles.poa.mntd.melt, comDataSESdamocles.poa.mpd.melt)
tail(comDataSESdamocles.poa.melt)

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/SES.1000.poales.pdf")
ggplot(comDataSESdamocles.poa.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Poales") 
dev.off()



######### Lamiales #############
persist.lam.EcrinsPool.com <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/Lamiales.persis.null_community_data.csv")
persist.lam.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/Lamiales.persis.summary_table.csv")

comDataSESdamocles.lam.mntd <- cbind(persist.lam.EcrinsPool.com %>% 
                                       mutate(RD = (mntd.RD - as.numeric(as.character(persist.lam.EcrinsPool.summary[11,3])))/ as.numeric(as.character(persist.lam.EcrinsPool.summary[12,3]))) %>% 
                                       mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(persist.lam.EcrinsPool.summary[11,3])))/ as.numeric(as.character(persist.lam.EcrinsPool.summary[12,3])))), 
                                     metric = "mntd",
                                     obs.RD = as.numeric(as.character(persist.lam.EcrinsPool.summary[13, 3])),
                                     obs.DAMO = as.numeric(as.character(persist.lam.EcrinsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                     clade = "Lamiales")

comDataSESdamocles.lam.mpd <- cbind(persist.lam.EcrinsPool.com %>%  
                                      mutate(RD = ((mpd.RD - as.numeric(as.character(persist.lam.EcrinsPool.summary[16,3])))/ as.numeric(as.character(persist.lam.EcrinsPool.summary[17,3])))) %>% 
                                      mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(persist.lam.EcrinsPool.summary[16,3])))/ as.numeric(as.character(persist.lam.EcrinsPool.summary[17,3])))), 
                                    metric = "mpd",
                                    obs.RD = as.numeric(as.character(persist.lam.EcrinsPool.summary[18, 3])),
                                    obs.DAMO = as.numeric(as.character(persist.lam.EcrinsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                    clade = "Lamiales")

comDataSESdamocles.lam.mntd.melt <- melt(comDataSESdamocles.lam.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.lam.mpd.melt <- melt(comDataSESdamocles.lam.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.lam.melt <- rbind(comDataSESdamocles.lam.mntd.melt, comDataSESdamocles.lam.mpd.melt)
tail(comDataSESdamocles.lam.melt)

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/SES.1000.Lamiales.pdf")
ggplot(comDataSESdamocles.lam.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Lamiales") 
dev.off()


######### Caryophyllales ############
persist.cary.EcrinsPool.com <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/caryophyllales.persis.null_community_data.csv")
persist.cary.EcrinsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/caryophyllales.persis.summary_table.csv")

comDataSESdamocles.cary.mntd <- cbind(persist.cary.EcrinsPool.com %>% 
                                        mutate(RD = (mntd.RD - as.numeric(as.character(persist.cary.EcrinsPool.summary[11,3])))/ as.numeric(as.character(persist.cary.EcrinsPool.summary[12,3]))) %>% 
                                        mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(persist.cary.EcrinsPool.summary[11,3])))/ as.numeric(as.character(persist.cary.EcrinsPool.summary[12,3])))), 
                                      metric = "mntd",
                                      obs.RD = as.numeric(as.character(persist.cary.EcrinsPool.summary[13, 3])),
                                      obs.DAMO = as.numeric(as.character(persist.cary.EcrinsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                      clade = "Caryophyllales")

comDataSESdamocles.cary.mpd <- cbind(persist.cary.EcrinsPool.com %>%  
                                       mutate(RD = ((mpd.RD - as.numeric(as.character(persist.cary.EcrinsPool.summary[16,3])))/ as.numeric(as.character(persist.cary.EcrinsPool.summary[17,3])))) %>% 
                                       mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(persist.cary.EcrinsPool.summary[16,3])))/ as.numeric(as.character(persist.cary.EcrinsPool.summary[17,3])))), 
                                     metric = "mpd",
                                     obs.RD = as.numeric(as.character(persist.cary.EcrinsPool.summary[18, 3])),
                                     obs.DAMO = as.numeric(as.character(persist.cary.EcrinsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                     clade = "Caryophyllales")

comDataSESdamocles.cary.mntd.melt <- melt(comDataSESdamocles.cary.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.cary.mpd.melt <- melt(comDataSESdamocles.cary.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.cary.melt <- rbind(comDataSESdamocles.cary.mntd.melt, comDataSESdamocles.cary.mpd.melt)
head(comDataSESdamocles.cary.melt)

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/SES.1000.Caryophyl.pdf")
ggplot(comDataSESdamocles.cary.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Caryophyllales") 
dev.off()

 
##### ALL clades together ########

master.ses.persis <- rbind(comDataSESdamocles.ast.mntd.melt, comDataSESdamocles.ast.mpd.melt,
                           comDataSESdamocles.poa.mntd.melt, comDataSESdamocles.poa.mpd.melt,
                           comDataSESdamocles.lam.mntd.melt, comDataSESdamocles.lam.mpd.melt,
                           comDataSESdamocles.cary.mntd.melt, comDataSESdamocles.cary.mpd.melt)
head(master.ses.persis)
master.ses.persis$value <- -1*master.ses.persis$value
master.ses.persis$obs.RD <- -1*master.ses.persis$obs.RD
master.ses.persis$obs.DAMO <- -1*master.ses.persis$obs.DAMO
summary(master.ses.persis)
master.ses.persis$metric <- factor(master.ses.persis$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))

ggplot(master.ses.persis, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(clade ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show.legend = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "") 

ggplot(filter(master.ses.persis, clade != "Poales"), aes(x = value, fill = variable)) + 
  geom_histogram(aes(y=0.5*..density.., alpha=variable), position='identity',binwidth=0.2) +
  facet_grid(clade ~ metric) +
  geom_vline(aes(xintercept=obs.RD), color = "black", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  #geom_vline(aes(xintercept=obs.DAMO), color = "green", show.legend = F, lty=2) + #
  labs(x = "SES", y = "Density", title = "") +
  scale_fill_manual(guide = guide_legend(title = "Null Models"), values = c("RD" = "#4477AA", "DAMOCLES" = "#AAAA44")) + 
  ggtitle("Ecrins Species Pool, LGM Community") + 
  scale_alpha_discrete(range = c("RD" = 1, "DAMOCLES" = .7), guide=FALSE) +
  scale_x_continuous(limits = c(-5,5)) + 
  scale_y_continuous(limits = c(0, .25)) +
  theme(strip.text.x = element_text(size = 12, face="bold"))
ggsave(filename="output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/SES.all.density.LGM.pdf")
ggsave(filename="figs/Figure4b_SES.all.density.LGM.pdf")

pdf(file="output/8_PhyoDiversity/alpha/dynamic/speciesPools/1_regionalEcrins/LGM/SES.all.density.green.pdf")
ggplot(master.ses.persis, aes(x = value, fill = variable)) + 
  geom_histogram(aes(y=0.5*..density..), alpha=0.7,position='identity',binwidth=0.5) +
  facet_grid(clade ~ metric) +
  geom_vline(aes(xintercept=obs.RD), color = "black", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  geom_vline(aes(xintercept=obs.DAMO), color = "green", show.legend = F, lty=2) + #
  labs(x = "SES", y = "Density", title = "") +
  scale_fill_discrete(guide = guide_legend(title = "Null Models")) + 
  ggtitle("Ecrins Species Pool, LGM Community") + 
  scale_x_continuous(limits = c(-5,5)) + 
  scale_y_continuous(limits = c(0, .25))
dev.off()


################################ All Summits Species Pool: LGM community ################################ 

######### Asterales #############
persist.ast.SummitsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/asterales.persis.summary_table.csv")
persist.ast.SummitsPool.com <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/asterales.persis.null_community_data.csv")

comDataSESdamocles.ast.mntd <- cbind(persist.ast.SummitsPool.com %>% 
                                       mutate(RD = (mntd.RD - as.numeric(as.character(persist.ast.SummitsPool.summary[11,3])))/ as.numeric(as.character(persist.ast.SummitsPool.summary[12,3]))) %>% 
                                       mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(persist.ast.SummitsPool.summary[11,3])))/ as.numeric(as.character(persist.ast.SummitsPool.summary[12,3])))), 
                                     metric = "mntd",
                                     obs.RD = as.numeric(as.character(persist.ast.SummitsPool.summary[13, 3])), #mntd.obs.RD.z.RD
                                     obs.DAMO = as.numeric(as.character(persist.ast.SummitsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                     clade = "Asterales")

comDataSESdamocles.ast.mpd <- cbind(persist.ast.SummitsPool.com %>%  
                                      mutate(RD = ((mpd.RD - as.numeric(as.character(persist.ast.SummitsPool.summary[16,3])))/ as.numeric(as.character(persist.ast.SummitsPool.summary[17,3])))) %>% 
                                      mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(persist.ast.SummitsPool.summary[16,3])))/ as.numeric(as.character(persist.ast.SummitsPool.summary[17,3])))), 
                                    metric = "mpd",
                                    obs.RD = as.numeric(as.character(persist.ast.SummitsPool.summary[18, 3])),
                                    obs.DAMO = as.numeric(as.character(persist.ast.SummitsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                    clade = "Asterales")

comDataSESdamocles.ast.mntd.melt <- melt(comDataSESdamocles.ast.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.ast.mpd.melt <- melt(comDataSESdamocles.ast.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.ast.melt <- rbind(comDataSESdamocles.ast.mntd.melt, comDataSESdamocles.ast.mpd.melt)
tail(comDataSESdamocles.ast.melt)

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/SES.1000.asterales.pdf")
ggplot(comDataSESdamocles.ast.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Asterales") 
dev.off()



######### Poales #############
persist.poa.SummitsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/poales.persis.summary_table.csv")
persist.ast.SummitsPool.com <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/poales.persis.null_community_data.csv")

comDataSESdamocles.poa.mntd <- cbind(persist.ast.SummitsPool.com %>% 
                                       mutate(RD = (mntd.RD - as.numeric(as.character(persist.poa.SummitsPool.summary[11,3])))/ as.numeric(as.character(persist.poa.SummitsPool.summary[12,3]))) %>% 
                                       mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(persist.poa.SummitsPool.summary[11,3])))/ as.numeric(as.character(persist.poa.SummitsPool.summary[12,3])))), 
                                     metric = "mntd",
                                     obs.RD = as.numeric(as.character(persist.poa.SummitsPool.summary[13, 3])), 
                                     obs.DAMO = as.numeric(as.character(persist.poa.SummitsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                     clade = "Poales")

comDataSESdamocles.poa.mpd <- cbind(persist.ast.SummitsPool.com %>%  
                                      mutate(RD = ((mpd.RD - as.numeric(as.character(persist.poa.SummitsPool.summary[16,3])))/ as.numeric(as.character(persist.poa.SummitsPool.summary[17,3])))) %>% 
                                      mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(persist.poa.SummitsPool.summary[16,3])))/ as.numeric(as.character(persist.poa.SummitsPool.summary[17,3])))), 
                                    metric = "mpd",
                                    obs.RD = as.numeric(as.character(persist.poa.SummitsPool.summary[18, 3])),
                                    obs.DAMO = as.numeric(as.character(persist.poa.SummitsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                    clade = "Poales")

comDataSESdamocles.poa.mntd.melt <- melt(comDataSESdamocles.poa.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.poa.mpd.melt <- melt(comDataSESdamocles.poa.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.poa.melt <- rbind(comDataSESdamocles.poa.mntd.melt, comDataSESdamocles.poa.mpd.melt)
tail(comDataSESdamocles.poa.melt)

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/SES.1000.poales.pdf")
ggplot(comDataSESdamocles.poa.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Poales") 
dev.off()



######### Lamiales #############
persist.lam.SummitsPool.com <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/Lamiales.persis.null_community_data.csv")
persist.lam.SummitsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/Lamiales.persis.summary_table.csv")

comDataSESdamocles.lam.mntd <- cbind(persist.lam.SummitsPool.com %>% 
                                       mutate(RD = (mntd.RD - as.numeric(as.character(persist.lam.SummitsPool.summary[11,3])))/ as.numeric(as.character(persist.lam.SummitsPool.summary[12,3]))) %>% 
                                       mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(persist.lam.SummitsPool.summary[11,3])))/ as.numeric(as.character(persist.lam.SummitsPool.summary[12,3])))), 
                                     metric = "mntd",
                                     obs.RD = as.numeric(as.character(persist.lam.SummitsPool.summary[13, 3])),
                                     obs.DAMO = as.numeric(as.character(persist.lam.SummitsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                     clade = "Lamiales")

comDataSESdamocles.lam.mpd <- cbind(persist.lam.SummitsPool.com %>%  
                                      mutate(RD = ((mpd.RD - as.numeric(as.character(persist.lam.SummitsPool.summary[16,3])))/ as.numeric(as.character(persist.lam.SummitsPool.summary[17,3])))) %>% 
                                      mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(persist.lam.SummitsPool.summary[16,3])))/ as.numeric(as.character(persist.lam.SummitsPool.summary[17,3])))), 
                                    metric = "mpd",
                                    obs.RD = as.numeric(as.character(persist.lam.SummitsPool.summary[18, 3])),
                                    obs.DAMO = as.numeric(as.character(persist.lam.SummitsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                    clade = "Lamiales")

comDataSESdamocles.lam.mntd.melt <- melt(comDataSESdamocles.lam.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.lam.mpd.melt <- melt(comDataSESdamocles.lam.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.lam.melt <- rbind(comDataSESdamocles.lam.mntd.melt, comDataSESdamocles.lam.mpd.melt)
tail(comDataSESdamocles.lam.melt)

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/SES.1000.Lamiales.pdf")
ggplot(comDataSESdamocles.lam.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Lamiales") 
dev.off()


######### Caryophyllales ############
persist.cary.SummitsPool.com <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/caryophyllales.persis.null_community_data.csv")
persist.cary.SummitsPool.summary <- read.csv("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/caryophyllales.persis.summary_table.csv")

comDataSESdamocles.cary.mntd <- cbind(persist.cary.SummitsPool.com %>% 
                                        mutate(RD = (mntd.RD - as.numeric(as.character(persist.cary.SummitsPool.summary[11,3])))/ as.numeric(as.character(persist.cary.SummitsPool.summary[12,3]))) %>% 
                                        mutate(DAMOCLES = ((mntd.DAMOCLES - as.numeric(as.character(persist.cary.SummitsPool.summary[11,3])))/ as.numeric(as.character(persist.cary.SummitsPool.summary[12,3])))), 
                                      metric = "mntd",
                                      obs.RD = as.numeric(as.character(persist.cary.SummitsPool.summary[13, 3])),
                                      obs.DAMO = as.numeric(as.character(persist.cary.SummitsPool.summary[24, 3])), #mntd.obs.z.DAMOCLES
                                      clade = "Caryophyllales")

comDataSESdamocles.cary.mpd <- cbind(persist.cary.SummitsPool.com %>%  
                                       mutate(RD = ((mpd.RD - as.numeric(as.character(persist.cary.SummitsPool.summary[16,3])))/ as.numeric(as.character(persist.cary.SummitsPool.summary[17,3])))) %>% 
                                       mutate(DAMOCLES = ((mpd.DAMOCLES - as.numeric(as.character(persist.cary.SummitsPool.summary[16,3])))/ as.numeric(as.character(persist.cary.SummitsPool.summary[17,3])))), 
                                     metric = "mpd",
                                     obs.RD = as.numeric(as.character(persist.cary.SummitsPool.summary[18, 3])),
                                     obs.DAMO = as.numeric(as.character(persist.cary.SummitsPool.summary[29, 3])), #mpd.obs.z.DAMOCLES
                                     clade = "Caryophyllales")

comDataSESdamocles.cary.mntd.melt <- melt(comDataSESdamocles.cary.mntd, measure.vars = c("RD", "DAMOCLES"))
comDataSESdamocles.cary.mpd.melt <- melt(comDataSESdamocles.cary.mpd, measure.vars = c("RD", "DAMOCLES"))

comDataSESdamocles.cary.melt <- rbind(comDataSESdamocles.cary.mntd.melt, comDataSESdamocles.cary.mpd.melt)
head(comDataSESdamocles.cary.melt)

pdf("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/SES.1000.Caryophyl.pdf")
ggplot(comDataSESdamocles.cary.melt, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(. ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show_guide = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "Caryophyllales") 
dev.off()


##### ALL clades together ########

master.ses.persis <- rbind(comDataSESdamocles.ast.mntd.melt, comDataSESdamocles.ast.mpd.melt,
                           comDataSESdamocles.poa.mntd.melt, comDataSESdamocles.poa.mpd.melt,
                           comDataSESdamocles.lam.mntd.melt, comDataSESdamocles.lam.mpd.melt,
                           comDataSESdamocles.cary.mntd.melt, comDataSESdamocles.cary.mpd.melt)
head(master.ses.persis)
master.ses.persis$value <- -1*master.ses.persis$value
master.ses.persis$obs.RD <- -1*master.ses.persis$obs.RD
master.ses.persis$obs.DAMO <- -1*master.ses.persis$obs.DAMO
summary(master.ses.persis)
master.ses.persis$metric <- factor(master.ses.persis$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))

ggplot(master.ses.persis, aes(x = variable, y=value, fill = variable)) + 
  geom_boxplot() +
  facet_grid(clade ~ metric) +
  geom_hline(aes(yintercept=obs.RD), color = "black", show.legend = T, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  labs(x = "Null Model", y = "SES", title = "") 

ggplot(filter(master.ses.persis, clade != "Poales"), aes(x = value, fill = variable)) + 
  geom_histogram(aes(y=0.5*..density.., alpha=variable), position='identity',binwidth=0.2) +
  facet_grid(clade ~ metric) +
  geom_vline(aes(xintercept=obs.RD), color = "black", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  #geom_vline(aes(xintercept=obs.DAMO), color = "green", show.legend = F, lty=2) + #
  labs(x = "SES", y = "Density", title = "") +
  scale_fill_manual(guide = guide_legend(title = "Null Models"), values = c("RD" = "#4477AA", "DAMOCLES" = "#AAAA44")) + 
  ggtitle("All Summits Species Pool, LGM Community") + 
  scale_alpha_discrete(range = c("RD" = 1, "DAMOCLES" = .7), guide=FALSE) +
  scale_x_continuous(limits = c(-5,5)) + 
  scale_y_continuous(limits = c(0, .25)) +
  theme(strip.text.x = element_text(size = 12, face="bold"))
ggsave(filename="output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/SES.all.density.LGM.pdf")
ggsave(filename="figs/Figure4c_SES.all.density.summitsLGM.pdf")

pdf(file="output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/LGM/output/SES.all.density.green.pdf")
ggplot(master.ses.persis, aes(x = value, fill = variable)) + 
  geom_histogram(aes(y=0.5*..density..), alpha=0.7,position='identity',binwidth=0.5) +
  facet_grid(clade ~ metric) +
  geom_vline(aes(xintercept=obs.RD), color = "black", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  geom_vline(aes(xintercept=obs.DAMO), color = "green", show.legend = F, lty=2) + #
  labs(x = "SES", y = "Density", title = "") +
  scale_fill_discrete(guide = guide_legend(title = "Null Models")) + 
  ggtitle("All Summits Species Pool, LGM Community") + 
  scale_x_continuous(limits = c(-5,5)) + 
  scale_y_continuous(limits = c(0, .25))
dev.off()










