#####################################################################################################################
############# Effect of species pool reduction on phylogenetic alpha diversity ######################################
############# Both Static & Dynamic Null Models #####################################################################
############# Hannah E. Marx, 6 June 2016 ###########################################################################
#####################################################################################################################

##################################  SES alpha diveristy of species pool reduction ##################################  
### Difference between SES metric for paired comparisons between source pools: 
#Paired t-test (use shapiro to test normal distribution of residuals)

master.ses.alpha <- read.csv(file="output/8_PhyoDiversity/alpha/static/Dryad_master.ses.static.alpha.csv", row.names=1)
head(master.ses.alpha)
master.ses.alpha.end <- master.ses.alpha[!master.ses.alpha$summits %in% c("Summits", "Persistent","Under Ice"),] #reove these...not interesting
master.ses.static.summtis <- master.ses.alpha.end

shapiro.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                         master.ses.static.summtis$pool == "Ecrins NP" & 
                                         master.ses.static.summtis$metric == "mntd", "obs.z"]) #This p-value tells you what the chances are that the sample comes from a normal distribution. The lower this value, the smaller the chance.

t.test.ecrin.summit.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                               master.ses.static.summtis$pool == "Ecrins NP" & 
                                                               master.ses.static.summtis$metric == "mntd", "obs.z"],
                                   master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                               master.ses.static.summtis$pool == "Summits" &
                                                               master.ses.static.summtis$metric == "mntd", "obs.z"], paired = T)

t.test.ecrin.summit.mntd.out <- cbind(test = "Pools", pool = "Ecrins_Summits", clade = "Spermatophyta", metric="mntd", mean.of.diff = t.test.ecrin.summit.mntd$estimate, t.value = t.test.ecrin.summit.mntd$statistic, df = t.test.ecrin.summit.mntd$parameter,
                                      p.value= t.test.ecrin.summit.mntd$p.value, conf.low = t.test.ecrin.summit.mntd$conf.int[1], conf.high =t.test.ecrin.summit.mntd$conf.int[2])


t.test.ecrin.summit.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                              master.ses.static.summtis$pool == "Ecrins NP" & 
                                                              master.ses.static.summtis$metric == "mpd", "obs.z"],
                                  master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                              master.ses.static.summtis$pool == "Summits" &
                                                              master.ses.static.summtis$metric == "mpd", "obs.z"], paired = T)

t.test.ecrin.summit.mpd.out <- cbind(test = "Pools", pool = "Ecrins_Summits", clade = "Spermatophyta", metric="mpd", mean.of.diff = t.test.ecrin.summit.mpd$estimate, t.value = t.test.ecrin.summit.mpd$statistic, df = t.test.ecrin.summit.mpd$parameter,
                                     p.value= t.test.ecrin.summit.mpd$p.value, conf.low = t.test.ecrin.summit.mpd$conf.int[1], conf.high =t.test.ecrin.summit.mpd$conf.int[2])


t.test.ecrin.LGM.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                            master.ses.static.summtis$pool == "Ecrins NP" & 
                                                            master.ses.static.summtis$metric == "mntd", "obs.z"],
                                master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                            master.ses.static.summtis$pool == "Persistent LGM" &
                                                            master.ses.static.summtis$metric == "mntd", "obs.z"], paired = T)
t.test.ecrin.LGM.mntd.out <- cbind(test = "Pools", pool = "Ecrins_LGM",clade = "Spermatophyta", metric="mntd", mean.of.diff = t.test.ecrin.LGM.mntd$estimate, t.value = t.test.ecrin.LGM.mntd$statistic, df = t.test.ecrin.LGM.mntd$parameter,
                                   p.value= t.test.ecrin.LGM.mntd$p.value, conf.low = t.test.ecrin.LGM.mntd$conf.int[1], conf.high =t.test.ecrin.LGM.mntd$conf.int[2])


t.test.ecrin.LGM.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                           master.ses.static.summtis$pool == "Ecrins NP" & 
                                                           master.ses.static.summtis$metric == "mpd", "obs.z"],
                               master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                           master.ses.static.summtis$pool == "Persistent LGM" &
                                                           master.ses.static.summtis$metric == "mpd", "obs.z"], paired = T)
t.test.ecrin.LGM.mpd.out <- cbind(test = "Pools", pool = "Ecrins_LGM", clade = "Spermatophyta", metric="mpd",  mean.of.diff = t.test.ecrin.LGM.mpd$estimate, t.value = t.test.ecrin.LGM.mpd$statistic, df = t.test.ecrin.LGM.mpd$parameter,
                                  p.value= t.test.ecrin.LGM.mpd$p.value, conf.low = t.test.ecrin.LGM.mpd$conf.int[1], conf.high =t.test.ecrin.LGM.mpd$conf.int[2])


##############

t.test.ecrin.clades.summit.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                      master.ses.static.summtis$pool == "Ecrins NP" & 
                                                                      master.ses.static.summtis$metric == "mntd", "obs.z"],
                                          master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                      master.ses.static.summtis$pool == "Summits" &
                                                                      master.ses.static.summtis$metric == "mntd", "obs.z"], paired = T)

t.test.ecrin.clades.summit.mntd.out <- cbind(test = "Pools", pool = "Ecrins_Summits", clade = "Clades", metric="mntd", mean.of.diff = t.test.ecrin.clades.summit.mntd$estimate, t.value = t.test.ecrin.clades.summit.mntd$statistic, df = t.test.ecrin.clades.summit.mntd$parameter,
                                             p.value= t.test.ecrin.clades.summit.mntd$p.value, conf.low = t.test.ecrin.clades.summit.mntd$conf.int[1], conf.high =t.test.ecrin.clades.summit.mntd$conf.int[2])


t.test.ecrin.clades.summit.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                     master.ses.static.summtis$pool == "Ecrins NP" & 
                                                                     master.ses.static.summtis$metric == "mpd", "obs.z"],
                                         master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                     master.ses.static.summtis$pool == "Summits" &
                                                                     master.ses.static.summtis$metric == "mpd", "obs.z"], paired = T)

t.test.ecrin.clades.summit.mpd.out <- cbind(test = "Pools", pool = "Ecrins_Summits",clade = "Clades",  metric="mpd", mean.of.diff = t.test.ecrin.clades.summit.mpd$estimate, t.value = t.test.ecrin.clades.summit.mpd$statistic, df = t.test.ecrin.clades.summit.mpd$parameter,
                                            p.value= t.test.ecrin.clades.summit.mpd$p.value, conf.low = t.test.ecrin.clades.summit.mpd$conf.int[1], conf.high =t.test.ecrin.clades.summit.mpd$conf.int[2])


t.test.ecrin.clades.LGM.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                   master.ses.static.summtis$pool == "Ecrins NP" & 
                                                                   master.ses.static.summtis$metric == "mntd", "obs.z"],
                                       master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                   master.ses.static.summtis$pool == "Persistent LGM" &
                                                                   master.ses.static.summtis$metric == "mntd", "obs.z"], paired = T)
t.test.ecrin.clades.LGM.mntd.out <- cbind(test = "Pools", pool = "Ecrins_LGM", clade = "Clades",metric="mntd", mean.of.diff = t.test.ecrin.clades.LGM.mntd$estimate, t.value = t.test.ecrin.clades.LGM.mntd$statistic, df = t.test.ecrin.clades.LGM.mntd$parameter,
                                          p.value= t.test.ecrin.clades.LGM.mntd$p.value, conf.low = t.test.ecrin.clades.LGM.mntd$conf.int[1], conf.high =t.test.ecrin.clades.LGM.mntd$conf.int[2])


t.test.ecrin.clades.LGM.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                  master.ses.static.summtis$pool == "Ecrins NP" & 
                                                                  master.ses.static.summtis$metric == "mpd", "obs.z"],
                                      master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                  master.ses.static.summtis$pool == "Persistent LGM" &
                                                                  master.ses.static.summtis$metric == "mpd", "obs.z"], paired = T)
t.test.ecrin.clades.LGM.mpd.out <- cbind(test = "Pools", pool = "Ecrins_LGM", clade = "Clades", metric="mpd", mean.of.diff.of.diff = t.test.ecrin.clades.LGM.mpd$estimate, t.value = t.test.ecrin.clades.LGM.mpd$statistic, df = t.test.ecrin.clades.LGM.mpd$parameter,
                                         p.value= t.test.ecrin.clades.LGM.mpd$p.value, conf.low = t.test.ecrin.clades.LGM.mpd$conf.int[1], conf.high =t.test.ecrin.clades.LGM.mpd$conf.int[2])

pools.stats <- as.data.frame(rbind(t.test.ecrin.summit.mntd.out, t.test.ecrin.summit.mpd.out, t.test.ecrin.LGM.mntd.out, t.test.ecrin.LGM.mpd.out, 
                    t.test.ecrin.clades.summit.mntd.out, t.test.ecrin.clades.summit.mpd.out, t.test.ecrin.clades.LGM.mntd.out, t.test.ecrin.clades.LGM.mpd.out))
rownames(pools.stats) <- NULL
head(pools.stats)
str(pools.stats)

pools.stats[5:9] <- apply(pools.stats[5:9], 2, as.character)
pools.stats[5:9] <- apply(pools.stats[5:9], 2, as.numeric)

pools.stats$p.value <- round(pools.stats$p.value, 9)

#write.csv(pools.stats, file="output/8_PhyoDiversity/alpha/static/sourcePools/diff.SES.pools.csv")

