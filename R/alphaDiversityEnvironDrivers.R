#####################################################################################################################
############# Environmental drivers of phylogenetic alpha diversity patterns ########################################
############# Environmental variables ~ alpha PD ####################################################################
############# Only Random Draw Null Model (too few observations for DAMOCLES) #######################################
############# Hannah E. Marx, 22 Mar 2016 ###########################################################################
#####################################################################################################################

source("R/envCategory.R") ## calculates PC of BioClim variables
source("R/plotAlphaEnv.R") ## plot regression

## Alpha diveristy output
master.ses.static <- read.csv(file="output/8_PhyoDiversity/alpha/static/Dryad_master.ses.static.alpha.csv", row.names=1)
head(master.ses.static)
master.ses.static <- cbind(master.ses.static, model = rep("RD", nrow(master.ses.static)))

master.ses.dynamic <- read.csv(file="output/8_PhyoDiversity/alpha/dynamic/Dryad_master.ses.dynamic.alpha.csv", row.names=1)
head(master.ses.dynamic)
master.ses.dynamic.damo <- master.ses.dynamic[master.ses.dynamic$model == "DAMOCLES",]

## combine ouput for both null models
alpha.div.master <- rbind(master.ses.static[-1], master.ses.dynamic.damo)
head(alpha.div.master)
str(alpha.div.master)

## geographic (spatial) distance
spatial.dist.data = na.omit(cbind(summit=(rownames(pezAlpes$env)), pezAlpes$env[c("X_WGS84", "Y_WGS84")]))
spatial.dist.tmp <- spatial.dist.data
coordinates(spatial.dist.tmp) <- ~X_WGS84 + Y_WGS84
spatial.dist <- gDistance(spatial.dist.tmp, byid = T)
min.d <- apply(spatial.dist, 1, function(x) order(x, decreasing=F)[2])
newdata <- cbind(spatial.dist.data, spatial.dist.data[min.d,], apply(spatial.dist, 1, function(x) sort(x, decreasing=F)[2]))
colnames(newdata) <- c(colnames(spatial.dist.data), 'neighbor', 'n.lat', 'n.long', 'distance')

spatial.dist <- newdata
alpha.div.master.spatial.dist <- (merge(alpha.div.master, spatial.dist,by.x=9, by.y=0))
head(alpha.div.master.spatial.dist)

## Raw environmental variables
alpha.div.master.env <- merge(alpha.div.master, alps.env.sprich.summits, by.x=9, by.y=0)
head(alpha.div.master.env)
dim(alpha.div.master.env)
names(alpha.div.master.env)
alpha.div.master.env$area <- log(alpha.div.master.env$area)
alpha.div.master.env$S.area <- log(alpha.div.master.env$S.area)
alpha.div.master.env$elevation <- log(alpha.div.master.env$elevation)


## PCA categorical variables
alpha.div.master.env.pca <- merge(alpha.div.master, env.pca, by.x=9, by.y=0)
head(alpha.div.master.env.pca)
names(alpha.div.master.env.pca)
str(alpha.div.master.env.pca$pc.availener.score)

## PCA all presnet and past
alpha.div.master.env.pca.recon <- merge(alpha.div.master, pca.lg, by.x=9, by.y=0)
head(alpha.div.master.env.pca.recon)
names(alpha.div.master.env.pca.recon)

## PCA differnece present- past
alpha.div.master.pca.lg.diff <- merge(alpha.div.master, env.pca.past.present, by.x=9, by.y=0)
head(alpha.div.master.pca.lg.diff)
names(alpha.div.master.pca.lg.diff)


#####################################################################################################################
#################################### Correltation ###################################################################
######################### Spearman's rank correlation rho

corfun<-function(x, y) {
  corr=(cor.test(alpha.div.master.env$obs.z, alpha.div.master.env$S.area, alternative="two.sided", method="spearman")) 
}
#The Spearman correlation is less sensitive than the Pearson correlation to strong outliers that are in the tails of both samples.
#The Pearson correlation is can deal with ties
#http://rstudio-pubs-static.s3.amazonaws.com/2292_412ffe4b0e024655a1ab971c173c0d8b.html
#https://www.reddit.com/r/statistics/comments/3alyyq/r_cortest_output/

#tests of no correlation
#estimate the association between paired samples and compute a test of the value being zero.
#alternative hypothesis: true rho is not equal to 0 
#http://www.r-bloggers.com/non-parametric-methods-for-the-study-of-the-correlation-spearmans-rank-correlation-coefficient-and-kendall-rho-rank-correlation-coefficient/
#The statistical test gives us as a result rho = 0.115, which indicates a low correlation (not parametric) between the two sets of values.
#The p-value > 0.05(??) allows us to accept the value of rho calculated, being statistically significant.

############################################### Geographic : Total summit area ############################################### 
geog_summitArea = ddply(alpha.div.master.env %>% filter(model == "RD"), 
                        .(metric, clade, pool), summarise, z=corfun(obs.z, S.area)$statistic,
                        pval=corfun(obs.z, S.area)$p.value,
                        rho.est=corfun(obs.z, S.area)$estimate,
                        alt=corfun(obs.z, S.area)$alternative)
#      CI.low=corfun(obs.z, S.area)$conf.int[1], 
#CI.high=corfun(obs.z, S.area)$conf.int[2], 

geog_summitArea <- cbind(geog_summitArea, hyp = rep("Geographic", times= nrow(geog_summitArea)), measure = rep("Summit Area", times= nrow(geog_summitArea)))


################################################  Geographic: Refugial area ###############################################################  
geog_area = ddply(alpha.div.master.env %>% filter(model == "RD"), 
                  .(metric, clade, pool), summarise, z=corfun(obs.z, area)$statistic,
                  pval=corfun(obs.z, area)$p.value,
                  rho.est=corfun(obs.z, area)$estimate,
                  #CI.low=corfun(obs.z, S.area)$conf.int[1], 
                  #CI.high=corfun(obs.z, S.area)$conf.int[2], 
                  alt=corfun(obs.z, area)$alternative) 
geog_area <- cbind(geog_area, hyp = rep("Geographic", times= nrow(geog_area)), measure = rep("Refugial Area", times= nrow(geog_area)))


################################################  Geographic: Maximum slope  ###############################################################  
geog_slope = ddply(alpha.div.master.env %>% filter(model == "RD"), 
                   .(metric, clade, pool), summarise, z=corfun(obs.z, max_slope)$statistic,
                   pval=corfun(obs.z, max_slope)$p.value,
                   rho.est=corfun(obs.z, max_slope)$estimate,
                   alt=corfun(obs.z, max_slope)$alternative) 

geog_slope <- cbind(geog_slope, hyp = rep("Geographic", times= nrow(geog_slope)), measure = rep("Maximum Slope", times= nrow(geog_slope)))


################################################  Geographic: Maximum elevation  ###############################################################  
geog_elevation = ddply(alpha.div.master.env %>% filter(model == "RD"), 
                       .(metric, clade, pool), summarise, z=corfun(obs.z, elevation)$statistic,
                       pval=corfun(obs.z, elevation)$p.value,
                       rho.est=corfun(obs.z, elevation)$estimate,
                       alt=corfun(obs.z, elevation)$alternative) 

geog_elevation <- cbind(geog_elevation, hyp = rep("Geographic", times= nrow(geog_elevation)), measure = rep("Maximum Elevation", times= nrow(geog_elevation)))


################################################  Geographic: Heterogeneity Topography ###############################################################  
geog_hetero =  ddply(alpha.div.master.env.pca %>% filter(model == "RD"), 
                     .(metric, clade, pool), summarise, z=corfun(obs.z, pc.hetero.score)$statistic,
                     pval=corfun(obs.z, pc.hetero.score)$p.value,
                     rho.est=corfun(obs.z, pc.hetero.score)$estimate,
                     alt=corfun(obs.z, pc.hetero.score)$alternative)

geog_hetero <- cbind(geog_hetero, hyp = rep("Geographic", times= nrow(geog_hetero)), measure = rep("Heterogeneity Topography", times= nrow(geog_hetero)))

################################################  Geographic: Heterogeneity of lithology  ###############################################################  
geog_hetero_lith =  ddply(alpha.div.master.env %>% filter(model == "RD"), 
                          .(metric, clade, pool), summarise, z=corfun(obs.z, simpson.d.geol)$statistic,
                          pval=corfun(obs.z, simpson.d.geol)$p.value,
                          rho.est=corfun(obs.z, simpson.d.geol)$estimate,
                          alt=corfun(obs.z, simpson.d.geol)$alternative)

geog_hetero_lith <- cbind(geog_hetero_lith, hyp = rep("Geographic", times= nrow(geog_hetero_lith)), measure = rep("Heterogeneity Lithology", times= nrow(geog_hetero_lith)))


############################################# Environment: Available Energy ############################################# 
env_energy =  ddply(alpha.div.master.env.pca %>% filter(model == "RD"), 
                    .(metric, clade, pool), summarise, z=corfun(obs.z, pc.availener.score)$statistic,
                    pval=corfun(obs.z, pc.availener.score)$p.value,
                    rho.est=corfun(obs.z, pc.availener.score)$estimate,
                    alt=corfun(obs.z, pc.availener.score)$alternative)
env_energy <- cbind(env_energy, hyp = rep("Environment", times= nrow(env_energy)), measure = rep("Available Energy", times= nrow(env_energy)))

############################################# Environment: Stress #############################################
env_stress = ddply(alpha.div.master.env.pca %>% filter(model == "RD"), 
                   .(metric, clade, pool), summarise, z=corfun(obs.z, pc.stress.score)$statistic,
                   pval=corfun(obs.z, pc.stress.score)$p.value,
                   rho.est=corfun(obs.z, pc.stress.score)$estimate,
                   alt=corfun(obs.z, pc.stress.score)$alternative)
env_stress <- cbind(env_stress, hyp = rep("Environment", times= nrow(env_stress)), measure = rep("Stress", times= nrow(env_stress)))

############################################# Environment: Present Stability #############################################
env_stabil = ddply(alpha.div.master.env.pca %>% filter(model == "RD"), 
                   .(metric, clade, pool), summarise, z=corfun(obs.z, pc.stabil.score)$statistic,
                   pval=corfun(obs.z, pc.stabil.score)$p.value,
                   rho.est=corfun(obs.z, pc.stabil.score)$estimate,
                   alt=corfun(obs.z, pc.stabil.score)$alternative)
env_stabil <- cbind(env_stabil, hyp = rep("Environment", times= nrow(env_stabil)), measure = rep("Stability", times= nrow(env_stabil)))


############################################# Environment: Past Stability (Present - Past Climate) ############################################# 
### Available Energy
env_stabilityPastEnergy =ddply(alpha.div.master.env.pca %>% filter(model == "RD"), 
                               .(metric, clade, pool), summarise, z=corfun(obs.z, present.cclgmbi.diff.availener.score)$statistic,
                               pval=corfun(obs.z, present.cclgmbi.diff.availener.score)$p.value,
                               rho.est=corfun(obs.z, present.cclgmbi.diff.availener.score)$estimate,
                               alt=corfun(obs.z, present.cclgmbi.diff.availener.score)$alternative)
env_stabilityPastEnergy <- cbind(env_stabilityPastEnergy, hyp = rep("Environment", times= nrow(env_stabilityPastEnergy)), measure = rep("Past Stability Energy", times= nrow(env_stabilityPastEnergy)))

### Stress
env_stabilityPastStress =ddply(alpha.div.master.env.pca %>% filter(model == "RD"), 
                               .(metric, clade, pool), summarise, z=corfun(obs.z, present.cclgmbi.diff.stress.score)$statistic,
                               pval=corfun(obs.z, present.cclgmbi.diff.stress.score)$p.value,
                               rho.est=corfun(obs.z, present.cclgmbi.diff.stress.score)$estimate,
                               alt=corfun(obs.z, present.cclgmbi.diff.stress.score)$alternative)
env_stabilityPastStress <- cbind(env_stabilityPastStress, hyp = rep("Environment", times= nrow(env_stabilityPastStress)), measure = rep("Past Stability Stress", times= nrow(env_stabilityPastStress)))

### Stability
env_stabilityPastStabil = ddply(alpha.div.master.env.pca %>% filter(model == "RD"), 
                                .(metric, clade, pool), summarise, z=corfun(obs.z, present.cclgmbi.diff.stabil.score)$statistic,
                                pval=corfun(obs.z, present.cclgmbi.diff.stabil.score)$p.value,
                                rho.est=corfun(obs.z, present.cclgmbi.diff.stabil.score)$estimate,
                                alt=corfun(obs.z, present.cclgmbi.diff.stabil.score)$alternative)
env_stabilityPastStabil <- cbind(env_stabilityPastStabil, hyp = rep("Environment", times= nrow(env_stabilityPastStabil)), measure = rep("Past Stability", times= nrow(env_stabilityPastStabil)))

################# 
env_correlate <- rbind(geog_summitArea, geog_area, geog_slope, geog_elevation,geog_hetero, geog_hetero_lith, env_energy, env_stress, env_stabil) #env_stabilityPastEnergy, env_stabilityPastStress, env_stabilityPastStabil

#write.csv(env_correlate, file="output/9_Environment/alpha/correlation/alpha_hyp_correlates.csv")
head(env_correlate)
str(env_correlate)

env_correlate$clade <- factor(env_correlate$clade, levels = c( "Spermatophyta", "Asterales", "Poales", "Rosales", "Lamiales", "Caryophyllales"))
#env_correlate$hyp <- factor(env_correlate$hyp, levels = c( "Geographic", "Environment"))
env_correlate$metric <- factor(env_correlate$metric, labels = c( "mntd" = "MNTD", "mpd"="MPD"))
env_correlate$pool <- factor(env_correlate$pool, levels = c( "Ecrins NP", "Summits", "Persistent LGM"))
env_correlate$pool <- factor(env_correlate$pool, labels = c( "Ecrins NP" = "Regional", "Summits" = "All Summits", "Persistent LGM" = "LGM"))
env_correlate$measure = with(env_correlate, factor(measure, levels = rev(levels(measure))))
env_correlate$measure <- factor(env_correlate$measure, levels = rev(c("Maximum Elevation", "Summit Area", "Refugial Area", 
                                                                      "Maximum Slope", "Heterogeneity Topography", "Heterogeneity Lithology",
                                                                      "Available Energy", "Stress", "Stability")))

env_correlate_plot <- ggplot(env_correlate, aes(y=measure, x=clade, fill=rho.est)) +
  geom_tile(alpha=0.75) +
  scale_fill_gradient2(low="#44AA77", high="#AA4488", mid="beige", na.value="white", limits=c(-.8, .8)) + # darkgreen, violet
  geom_point(aes(size=ifelse(pval <= 0.05, "dot", "no_dot"))) +
  scale_size_manual(values=c(dot=2, no_dot=NA), guide="none") +
  facet_grid(pool ~ metric) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 12, face="bold"),
    strip.text.y = element_text(size = 12, face="bold"),
    axis.text.x = element_text(angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "rho")) +
  coord_fixed(ratio=1)
ggsave("output/9_Environment/alpha/correlation/alpha_env_correlates.pdf", env_correlate_plot, width = 8.5, height = 11)
ggsave("figs/Figure4_alpha_env_correlates.pdf", env_correlate_plot, width = 8.5, height = 11)



#####################################################################################################################
#################################### Regression ###################################################################

############################################### Geographic : Total summit area ############################################### 
geog_summitArea = alpha.div.master.env %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~S.area, data=.)  # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

geog_summitArea.df <- as.data.frame(geog_summitArea %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_summitArea.df, file="output/9_Environment/alpha/regression/geog_summitArea.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env), indepVar="S.area", xx = 12.5, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/alpha.RD.AllClades.geog_summitArea.pdf")



################################################  Geographic: Refugial area ###############################################################  
geog_area = alpha.div.master.env %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~area, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

geog_area.df <- as.data.frame(geog_area %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_area.df, file="output/9_Environment/alpha/regression/geog_area.csv")
## column area is slope (coefficient of variable area)

plotAlphaEnv(df=as.data.frame(alpha.div.master.env), indepVar="area", xx = 12.5, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/alpha.RD.AllClades.geog_area.pdf")


################################################  Geographic: Maximum slope  ###############################################################  
geog_slope = alpha.div.master.env %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~max_slope, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

geog_slope.df <- as.data.frame(geog_slope %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_slope.df, file="output/9_Environment/alpha/regression/geog_slope.csv")
## column area is slope (coefficient of variable area)

plotAlphaEnv(df=as.data.frame(alpha.div.master.env), indepVar="max_slope", xx = 68, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/alpha.RD.AllClades.geog_max_slope.pdf")

################################################  Geographic: Heterogeneity  ###############################################################  
geog_hetero = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.hetero.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

geog_hetero.df <- as.data.frame(geog_hetero %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_hetero.df, file="output/9_Environment/alpha/regression/geog_hetero.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.hetero.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/alpha.RD.AllClades.geog_hetero.pdf")

############################################# Environment: Available Energy ############################################# 
env_energy = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.availener.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_energy.df <- as.data.frame(env_energy %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_energy.df, file="output/9_Environment/alpha/regression/env_energy.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.availener.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/alpha.RD.AllClades.env_energy.pdf")

############################################# Environment: Stress #############################################
env_stress = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.stress.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_stress.df <- as.data.frame(env_stress %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stress.df, file="output/9_Environment/alpha/regression/env_stress.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.stress.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/alpha.RD.AllClades.env_stress.pdf")
plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.stress.score", xx = 0, dot=F, eqn=F, filename="output/9_Environment/alpha/regression/alpha.RD.AllClades.env_stress_notxt.pdf")


############################################# Environment: Present Stability #############################################
env_stabil = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.stabil.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_stabil.df <- as.data.frame(env_stabil %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabil.df, file="output/9_Environment/alpha/regression/env_stabil.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="pc.stabil.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/alpha.RD.AllClades.env_stability.pdf")

############################################# Environment: Past Stability (Present - Past Climate) ############################################# 
### Available Energy
env_stabilityPastEnergy = alpha.div.master.pca.lg.diff %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~present.cclgmbi.diff.availener.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_stabilityPastEnergy.df <- as.data.frame(env_stabilityPastEnergy %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabilityPastEnergy.df, file="output/9_Environment/alpha/regression/env_stabilityPastEnergy.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="present.cclgmbi.diff.availener.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/alpha.RD.AllClades.env_stabilityPastEnergy.pdf")

### Stress
env_stabilityPastStress = alpha.div.master.pca.lg.diff %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~present.cclgmbi.diff.stress.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_stabilityPastStress.df <- as.data.frame(env_stabilityPastStress %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabilityPastStress.df, file="output/9_Environment/alpha/regression/env_stabilityPastStress.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="present.cclgmbi.diff.stress.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/alpha.RD.AllClades.env_stabilityPastStress.pdf")

### Stability
env_stabilityPastStabil = alpha.div.master.pca.lg.diff %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~present.cclgmbi.diff.stabil.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_stabilityPastStabil.df <- as.data.frame(env_stabilityPastStabil %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabilityPastStabil.df, file="output/9_Environment/alpha/regression/env_stabilityPastStabil.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="present.cclgmbi.diff.stabil.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/alpha.RD.AllClades.env_stabilityPastStabil.pdf")






