
##################################  SES alpha diveristy of species pool reduction ##################################  
### Difference between SES metric for paired comparisons between source pools: 
#Paired t-test (use shapiro to test normal distribution of residuals)

master.ses.alpha <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/master.ses.static.alpha.csv", row.names=1)
head(master.ses.alpha)
master.ses.alpha.end <- master.ses.alpha[!master.ses.alpha$summits %in% c("Summits", "Persistent","Under Ice"),] #reove this...not interesting
master.ses.static.summtis <- master.ses.alpha.end

head(master.ses.static)
#master.ses.static.summtis <- master.ses.static[master.ses.static$summits %in% rownames(alps.env.sprich.summits),]

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





##################################  SES alpha diveristy with and without endemics ##################################  

master.ses.static.NOendemics <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/Endemics/master.SES.NOendemic.csv", row.names=1)
head(master.ses.static.NOendemics)

master.ses.static.NOendemics <- master.ses.static.NOendemics[master.ses.static.NOendemics$summits %in% rownames(alps.env.sprich.summits),]
unique(master.ses.static.NOendemics$summits)
unique(master.ses.static.summtis$summits)


t.test.sperm.endemic.no.ecrins.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                               master.ses.static.summtis$pool == "Ecrins NP" & 
                                                               master.ses.static.summtis$metric == "mntd", "obs.z"],
                                   master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                  master.ses.static.NOendemics$pool == "Ecrins NP" &
                                                                  master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.sperm.endemic.no.ecrins.mntd.out <- cbind(test = "Endemics_No", pool = "Ecrins", clade = "Spermatophyta", metric="mntd", mean.of.diff = t.test.sperm.endemic.no.ecrins.mntd$estimate, t.value = t.test.sperm.endemic.no.ecrins.mntd$statistic, df = t.test.sperm.endemic.no.ecrins.mntd$parameter,
                                      p.value= t.test.sperm.endemic.no.ecrins.mntd$p.value, conf.low = t.test.sperm.endemic.no.ecrins.mntd$conf.int[1], conf.high =t.test.sperm.endemic.no.ecrins.mntd$conf.int[2])


t.test.sperm.endemic.no.ecrins.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                                          master.ses.static.summtis$pool == "Ecrins NP" & 
                                                                          master.ses.static.summtis$metric == "mpd", "obs.z"],
                                              master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                             master.ses.static.NOendemics$pool == "Ecrins NP" &
                                                                             master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.sperm.endemic.no.ecrins.mpd.out <- cbind(test = "Endemics_No", pool = "Ecrins", clade = "Spermatophyta", metric="mpd", mean.of.diff = t.test.sperm.endemic.no.ecrins.mpd$estimate, t.value = t.test.sperm.endemic.no.ecrins.mpd$statistic, df = t.test.sperm.endemic.no.ecrins.mpd$parameter,
                                                 p.value= t.test.sperm.endemic.no.ecrins.mpd$p.value, conf.low = t.test.sperm.endemic.no.ecrins.mpd$conf.int[1], conf.high =t.test.sperm.endemic.no.ecrins.mpd$conf.int[2])




t.test.sperm.endemic.no.summit.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                                          master.ses.static.summtis$pool == "Summits" & 
                                                                          master.ses.static.summtis$metric == "mntd", "obs.z"],
                                              master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                             master.ses.static.NOendemics$pool == "Summits" &
                                                                             master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.sperm.endemic.no.summit.mntd.out <- cbind(test = "Endemics_No", pool = "summit", clade = "Spermatophyta", metric="mntd", mean.of.diff = t.test.sperm.endemic.no.summit.mntd$estimate, t.value = t.test.sperm.endemic.no.summit.mntd$statistic, df = t.test.sperm.endemic.no.summit.mntd$parameter,
                                                 p.value= t.test.sperm.endemic.no.summit.mntd$p.value, conf.low = t.test.sperm.endemic.no.summit.mntd$conf.int[1], conf.high =t.test.sperm.endemic.no.summit.mntd$conf.int[2])


t.test.sperm.endemic.no.summit.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                                         master.ses.static.summtis$pool == "Summits" & 
                                                                         master.ses.static.summtis$metric == "mpd", "obs.z"],
                                             master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                            master.ses.static.NOendemics$pool == "Summits" &
                                                                            master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.sperm.endemic.no.summit.mpd.out <- cbind(test = "Endemics_No", pool = "summit", clade = "Spermatophyta", metric="mpd", mean.of.diff = t.test.sperm.endemic.no.summit.mpd$estimate, t.value = t.test.sperm.endemic.no.summit.mpd$statistic, df = t.test.sperm.endemic.no.summit.mpd$parameter,
                                                p.value= t.test.sperm.endemic.no.summit.mpd$p.value, conf.low = t.test.sperm.endemic.no.summit.mpd$conf.int[1], conf.high =t.test.sperm.endemic.no.summit.mpd$conf.int[2])




t.test.sperm.endemic.no.persistent.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                                          master.ses.static.summtis$pool == "Persistent LGM" & 
                                                                          master.ses.static.summtis$metric == "mntd", "obs.z"],
                                              master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                             master.ses.static.NOendemics$pool == "Persistent LGM" &
                                                                             master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.sperm.endemic.no.persistent.mntd.out <- cbind(test = "Endemics_No", pool = "Persistent LGM", clade = "Spermatophyta", metric="mntd", mean.of.diff = t.test.sperm.endemic.no.persistent.mntd$estimate, t.value = t.test.sperm.endemic.no.persistent.mntd$statistic, df = t.test.sperm.endemic.no.persistent.mntd$parameter,
                                                 p.value= t.test.sperm.endemic.no.persistent.mntd$p.value, conf.low = t.test.sperm.endemic.no.persistent.mntd$conf.int[1], conf.high =t.test.sperm.endemic.no.persistent.mntd$conf.int[2])


t.test.sperm.endemic.no.persistent.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                                         master.ses.static.summtis$pool == "Persistent LGM" & 
                                                                         master.ses.static.summtis$metric == "mpd", "obs.z"],
                                             master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                            master.ses.static.NOendemics$pool == "Persistent LGM" &
                                                                            master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.sperm.endemic.no.persistent.mpd.out <- cbind(test = "Endemics_No", pool = "Persistent LGM", clade = "Spermatophyta", metric="mpd", mean.of.diff = t.test.sperm.endemic.no.persistent.mpd$estimate, t.value = t.test.sperm.endemic.no.persistent.mpd$statistic, df = t.test.sperm.endemic.no.persistent.mpd$parameter,
                                                p.value= t.test.sperm.endemic.no.persistent.mpd$p.value, conf.low = t.test.sperm.endemic.no.persistent.mpd$conf.int[1], conf.high =t.test.sperm.endemic.no.persistent.mpd$conf.int[2])


###########

t.test.clades.endemic.no.ecrins.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                            master.ses.static.summtis$pool == "Ecrins NP" & 
                                                                            master.ses.static.summtis$metric == "mntd", "obs.z"],
                                                master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                               master.ses.static.NOendemics$pool == "Ecrins NP" &
                                                                               master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.clades.endemic.no.ecrins.mntd.out <- cbind(test = "Endemics_No", pool = "Ecrins NP", clade = "Clades", metric="mntd", mean.of.diff = t.test.clades.endemic.no.ecrins.mntd$estimate, t.value = t.test.clades.endemic.no.ecrins.mntd$statistic, df = t.test.clades.endemic.no.ecrins.mntd$parameter,
                                                   p.value= t.test.clades.endemic.no.ecrins.mntd$p.value, conf.low = t.test.clades.endemic.no.ecrins.mntd$conf.int[1], conf.high =t.test.clades.endemic.no.ecrins.mntd$conf.int[2])


t.test.clades.endemic.no.ecrins.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                           master.ses.static.summtis$pool == "Ecrins NP" & 
                                                                           master.ses.static.summtis$metric == "mpd", "obs.z"],
                                               master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                              master.ses.static.NOendemics$pool == "Ecrins NP" &
                                                                              master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.clades.endemic.no.ecrins.mpd.out <- cbind(test = "Endemics_No", pool = "Ecrins NP", clade = "Clades", metric="mpd", mean.of.diff = t.test.clades.endemic.no.ecrins.mpd$estimate, t.value = t.test.clades.endemic.no.ecrins.mpd$statistic, df = t.test.clades.endemic.no.ecrins.mpd$parameter,
                                                  p.value= t.test.clades.endemic.no.ecrins.mpd$p.value, conf.low = t.test.clades.endemic.no.ecrins.mpd$conf.int[1], conf.high =t.test.clades.endemic.no.ecrins.mpd$conf.int[2])






t.test.clades.endemic.no.summits.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                               master.ses.static.summtis$pool == "Summits" & 
                                                                               master.ses.static.summtis$metric == "mntd", "obs.z"],
                                                   master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                                  master.ses.static.NOendemics$pool == "Summits" &
                                                                                  master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.clades.endemic.no.summits.mntd.out <- cbind(test = "Endemics_No", pool = "Summits", clade = "Clades", metric="mntd", mean.of.diff = t.test.clades.endemic.no.summits.mntd$estimate, t.value = t.test.clades.endemic.no.summits.mntd$statistic, df = t.test.clades.endemic.no.summits.mntd$parameter,
                                                      p.value= t.test.clades.endemic.no.summits.mntd$p.value, conf.low = t.test.clades.endemic.no.summits.mntd$conf.int[1], conf.high =t.test.clades.endemic.no.summits.mntd$conf.int[2])


t.test.clades.endemic.no.summits.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                              master.ses.static.summtis$pool == "Summits" & 
                                                                              master.ses.static.summtis$metric == "mpd", "obs.z"],
                                                  master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                                 master.ses.static.NOendemics$pool == "Summits" &
                                                                                 master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.clades.endemic.no.summits.mpd.out <- cbind(test = "Endemics_No", pool = "Summits", clade = "Clades", metric="mpd", mean.of.diff = t.test.clades.endemic.no.summits.mpd$estimate, t.value = t.test.clades.endemic.no.summits.mpd$statistic, df = t.test.clades.endemic.no.summits.mpd$parameter,
                                                     p.value= t.test.clades.endemic.no.summits.mpd$p.value, conf.low = t.test.clades.endemic.no.summits.mpd$conf.int[1], conf.high =t.test.clades.endemic.no.summits.mpd$conf.int[2])








t.test.clades.endemic.no.persistent.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                              master.ses.static.summtis$pool == "Persistent LGM" & 
                                                                              master.ses.static.summtis$metric == "mntd", "obs.z"],
                                                  master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                                 master.ses.static.NOendemics$pool == "Persistent LGM" &
                                                                                 master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.clades.endemic.no.persistent.mntd.out <- cbind(test = "Endemics_No", pool = "Persistent LGM", clade = "Clades", metric="mntd", mean.of.diff = t.test.clades.endemic.no.persistent.mntd$estimate, t.value = t.test.clades.endemic.no.persistent.mntd$statistic, df = t.test.clades.endemic.no.persistent.mntd$parameter,
                                                     p.value= t.test.clades.endemic.no.persistent.mntd$p.value, conf.low = t.test.clades.endemic.no.persistent.mntd$conf.int[1], conf.high =t.test.clades.endemic.no.persistent.mntd$conf.int[2])


t.test.clades.endemic.no.persistent.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                             master.ses.static.summtis$pool == "Persistent LGM" & 
                                                                             master.ses.static.summtis$metric == "mpd", "obs.z"],
                                                 master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                                master.ses.static.NOendemics$pool == "Persistent LGM" &
                                                                                master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.clades.endemic.no.persistent.mpd.out <- cbind(test = "Endemics_No", pool = "Persistent LGM", clade = "Clades", metric="mpd", mean.of.diff = t.test.clades.endemic.no.persistent.mpd$estimate, t.value = t.test.clades.endemic.no.persistent.mpd$statistic, df = t.test.clades.endemic.no.persistent.mpd$parameter,
                                                    p.value= t.test.clades.endemic.no.persistent.mpd$p.value, conf.low = t.test.clades.endemic.no.persistent.mpd$conf.int[1], conf.high =t.test.clades.endemic.no.persistent.mpd$conf.int[2])


################
#The sign of a t-value tells us the direction of the difference in sample means

biotic.corr.stats <- as.data.frame(rbind(t.test.ecrin.summit.mntd.out, t.test.ecrin.summit.mpd.out, t.test.ecrin.LGM.mntd.out, t.test.ecrin.LGM.mpd.out, 
      t.test.ecrin.clades.summit.mntd.out, t.test.ecrin.clades.summit.mpd.out, t.test.ecrin.clades.LGM.mntd.out, t.test.ecrin.clades.LGM.mpd.out,
      t.test.sperm.endemic.no.ecrins.mntd.out, t.test.sperm.endemic.no.ecrins.mpd.out, t.test.sperm.endemic.no.summit.mntd.out,
      t.test.sperm.endemic.no.summit.mpd.out, t.test.sperm.endemic.no.persistent.mntd.out, t.test.sperm.endemic.no.persistent.mpd.out,
      t.test.clades.endemic.no.ecrins.mntd.out, t.test.clades.endemic.no.ecrins.mpd.out, t.test.clades.endemic.no.summits.mntd.out,
      t.test.clades.endemic.no.summits.mpd.out, t.test.clades.endemic.no.persistent.mntd.out, t.test.clades.endemic.no.persistent.mpd.out))
rownames(biotic.corr.stats) <- NULL
head(biotic.corr.stats)
str(biotic.corr.stats)

biotic.corr.stats[5:9] <- apply(biotic.corr.stats[5:9], 2, as.character)
biotic.corr.stats[5:9] <- apply(biotic.corr.stats[5:9], 2, as.numeric)

biotic.corr.stats$p.value <- round(biotic.corr.stats$p.value, 9)

#write.csv(biotic.corr.stats, file="output/9_PhyoDiversity/Spermatophyta/diff.SES.pools.endemics.csv")


