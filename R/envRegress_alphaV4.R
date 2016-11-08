### Environmental variables ~ alpha PD)
### Only Random Draw Null Model (too few observations for DAMOCLES)

source("R/envCategory.R")
source("R/plotAlphaEnv.R")

## Alpha diveristy output
master.ses.static <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/master.ses.static.alpha.csv", row.names=1)
head(master.ses.static)
master.ses.static <- cbind(master.ses.static, model = rep("RD", nrow(master.ses.static)))

master.ses.dynamic <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/dynamic/master.ses.dynamic.alpha.csv", row.names=1)
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
alpha.div.master.env <- merge(alpha.div.master, alps.env.sprich, by.x=9, by.y=0)
head(alpha.div.master.env)
dim(alpha.div.master.env)
names(alpha.div.master.env)
alpha.div.master.env$area <- log(alpha.div.master.env$area)
alpha.div.master.env$S.area <- log(alpha.div.master.env$S.area)

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


############################################### Geographic : Total summit area ############################################### 
geog_summitArea = alpha.div.master.env %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~S.area, data=.)  # create your model
      data.frame(tidy(model),              # get coefficient info
                 glance(model))})

geog_summitArea.df <- as.data.frame(geog_summitArea %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_summitArea.df, file="output/10_Environment/alpha/geog_summitArea.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env), indepVar="S.area", xx = 12.5, dot=T, eqn=T, filename="output/10_Environment/alpha/alpha.RD.AllClades.geog_summitArea.pdf")



################################################  Geographic: Refugial area ###############################################################  
geog_area = alpha.div.master.env %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~area, data=.)    # create your model
      data.frame(tidy(model),              # get coefficient info
                 glance(model))})

geog_area.df <- as.data.frame(geog_area %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_area.df, file="output/10_Environment/alpha/geog_area.csv")
## column area is slope (coefficient of variable area)

plotAlphaEnv(df=as.data.frame(alpha.div.master.env), indepVar="area", xx = 12.5, dot=T, eqn=T, filename="output/10_Environment/alpha/alpha.RD.AllClades.geog_area.pdf")


################################################  Geographic: Maximum slope  ###############################################################  
geog_slope = alpha.div.master.env %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~max_slope, data=.)    # create your model
      data.frame(tidy(model),              # get coefficient info
                 glance(model))})

geog_slope.df <- as.data.frame(geog_slope %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_slope.df, file="output/10_Environment/alpha/geog_slope.csv")
## column area is slope (coefficient of variable area)

plotAlphaEnv(df=as.data.frame(alpha.div.master.env), indepVar="max_slope", xx = 68, dot=T, eqn=T, filename="output/10_Environment/alpha/alpha.RD.AllClades.geog_max_slope.pdf")

################################################  Geographic: Heterogeneity  ###############################################################  
geog_hetero = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.hetero.score, data=.)    # create your model
      data.frame(tidy(model),              # get coefficient info
                 glance(model))})

geog_hetero.df <- as.data.frame(geog_hetero %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_hetero.df, file="output/10_Environment/alpha/geog_hetero.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.hetero.score", xx = 0, dot=T, eqn=T, filename="output/10_Environment/alpha/alpha.RD.AllClades.geog_hetero.pdf")

############################################# Environment: Available Energy ############################################# 
env_energy = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.availener.score, data=.)    # create your model
      data.frame(tidy(model),              # get coefficient info
                 glance(model))})

env_energy.df <- as.data.frame(env_energy %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_energy.df, file="output/10_Environment/alpha/env_energy.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.availener.score", xx = 0, dot=T, eqn=T, filename="output/10_Environment/alpha/alpha.RD.AllClades.env_energy.pdf")

############################################# Environment: Stress #############################################
env_stress = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.stress.score, data=.)    # create your model
      data.frame(tidy(model),              # get coefficient info
                 glance(model))})

env_stress.df <- as.data.frame(env_stress %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stress.df, file="output/10_Environment/alpha/env_stress.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.stress.score", xx = 0, dot=T, eqn=T, filename="output/10_Environment/alpha/alpha.RD.AllClades.env_stress.pdf")
plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.stress.score", xx = 0, dot=F, eqn=F, filename="output/10_Environment/alpha/alpha.RD.AllClades.env_stress_notxt.pdf")


############################################# Environment: Present Stability #############################################
env_stabil = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.stabil.score, data=.)    # create your model
      data.frame(tidy(model),              # get coefficient info
                 glance(model))})

env_stabil.df <- as.data.frame(env_stabil %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabil.df, file="output/10_Environment/alpha/env_stabil.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="pc.stabil.score", xx = 0, dot=T, eqn=T, filename="output/10_Environment/alpha/alpha.RD.AllClades.env_stability.pdf")

############################################# Environment: Past Stability (Present - Past Climate) ############################################# 
### Available Energy
env_stabilityPastEnergy = alpha.div.master.pca.lg.diff %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~present.cclgmbi.diff.availener.score, data=.)    # create your model
      data.frame(tidy(model),              # get coefficient info
                 glance(model))})

env_stabilityPastEnergy.df <- as.data.frame(env_stabilityPastEnergy %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabilityPastEnergy.df, file="output/10_Environment/alpha/env_stabilityPastEnergy.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="present.cclgmbi.diff.availener.score", xx = 0, dot=T, eqn=T, filename="output/10_Environment/alpha/alpha.RD.AllClades.env_stabilityPastEnergy.pdf")

### Stress
env_stabilityPastStress = alpha.div.master.pca.lg.diff %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~present.cclgmbi.diff.stress.score, data=.)    # create your model
      data.frame(tidy(model),              # get coefficient info
                 glance(model))})

env_stabilityPastStress.df <- as.data.frame(env_stabilityPastStress %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabilityPastStress.df, file="output/10_Environment/alpha/env_stabilityPastStress.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="present.cclgmbi.diff.stress.score", xx = 0, dot=T, eqn=T, filename="output/10_Environment/alpha/alpha.RD.AllClades.env_stabilityPastStress.pdf")

### Stability
env_stabilityPastStabil = alpha.div.master.pca.lg.diff %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~present.cclgmbi.diff.stabil.score, data=.)    # create your model
      data.frame(tidy(model),              # get coefficient info
                 glance(model))})

env_stabilityPastStabil.df <- as.data.frame(env_stabilityPastStabil %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabilityPastStabil.df, file="output/10_Environment/alpha/env_stabilityPastStabil.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="present.cclgmbi.diff.stabil.score", xx = 0, dot=T, eqn=T, filename="output/10_Environment/alpha/alpha.RD.AllClades.env_stabilityPastStabil.pdf")






