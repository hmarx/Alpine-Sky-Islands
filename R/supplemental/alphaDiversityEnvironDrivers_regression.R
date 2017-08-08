#####################################################################################################################
############# Environmental predictors of phylogenetic diversity  ###################################################
############# Hannah E. Marx, 25 April 2017 #########################################################################
#####################################################################################################################

############# Regression of SES on environmental variables ##########################################################
############# Observed SES on each summit (15 points) ~ environmental variables

source("R/alphaDiversityEnvironDrivers_correlations.R") ## generates data objects 
source("R/plotAlphaEnv.R") ## plot regression

############################################### Geographic : Total summit area ############################################### 
geog_summitArea = alpha.div.master.env %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~S.area, data=.)  # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

geog_summitArea.df <- as.data.frame(geog_summitArea %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_summitArea.df, file="output/9_Environment/alpha/regression/observed/geog_summitArea.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env), indepVar="S.area", xx = 12.5, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/observed/alpha.RD.AllClades.geog_summitArea.pdf")



################################################  Geographic: Refugial area ###############################################################  
geog_area = alpha.div.master.env %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~area, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

geog_area.df <- as.data.frame(geog_area %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_area.df, file="output/9_Environment/alpha/regression/observed/geog_area.csv")
## column area is slope (coefficient of variable area)

plotAlphaEnv(df=as.data.frame(alpha.div.master.env), indepVar="area", xx = 12.5, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/observed/alpha.RD.AllClades.geog_area.pdf")


################################################  Geographic: Maximum slope  ###############################################################  
geog_slope = alpha.div.master.env %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~max_slope, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

geog_slope.df <- as.data.frame(geog_slope %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_slope.df, file="output/9_Environment/alpha/regression/observed/geog_slope.csv")
## column area is slope (coefficient of variable area)

plotAlphaEnv(df=as.data.frame(alpha.div.master.env), indepVar="max_slope", xx = 68, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/observed/alpha.RD.AllClades.geog_max_slope.pdf")

################################################  Geographic: Heterogeneity  ###############################################################  
geog_hetero = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.hetero.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

geog_hetero.df <- as.data.frame(geog_hetero %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(geog_hetero.df, file="output/9_Environment/alpha/regression/observed/geog_hetero.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.hetero.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/observed/alpha.RD.AllClades.geog_hetero.pdf")

############################################# Environment: Available Energy ############################################# 
env_energy = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.availener.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_energy.df <- as.data.frame(env_energy %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_energy.df, file="output/9_Environment/alpha/regression/observed/env_energy.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.availener.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/observed/alpha.RD.AllClades.env_energy.pdf")

############################################# Environment: Stress #############################################
env_stress = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.stress.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_stress.df <- as.data.frame(env_stress %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stress.df, file="output/9_Environment/alpha/regression/observed/env_stress.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.stress.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/observed/alpha.RD.AllClades.env_stress.pdf")
plotAlphaEnv(df=as.data.frame(alpha.div.master.env.pca), indepVar="pc.stress.score", xx = 0, dot=F, eqn=F, filename="output/9_Environment/alpha/regression/observed/alpha.RD.AllClades.env_stress_notxt.pdf")


############################################# Environment: Present Stability #############################################
env_stabil = alpha.div.master.env.pca %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~pc.stabil.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_stabil.df <- as.data.frame(env_stabil %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabil.df, file="output/9_Environment/alpha/regression/observed/env_stabil.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="pc.stabil.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/observed/alpha.RD.AllClades.env_stability.pdf")

############################################# Environment: Past Stability (Present - Past Climate) ############################################# 
### Available Energy
env_stabilityPastEnergy = alpha.div.master.pca.lg.diff %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~present.cclgmbi.diff.availener.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_stabilityPastEnergy.df <- as.data.frame(env_stabilityPastEnergy %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabilityPastEnergy.df, file="output/9_Environment/alpha/regression/observed/env_stabilityPastEnergy.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="present.cclgmbi.diff.availener.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/observed/alpha.RD.AllClades.env_stabilityPastEnergy.pdf")

### Stress
env_stabilityPastStress = alpha.div.master.pca.lg.diff %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~present.cclgmbi.diff.stress.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_stabilityPastStress.df <- as.data.frame(env_stabilityPastStress %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabilityPastStress.df, file="output/9_Environment/alpha/regression/observed/env_stabilityPastStress.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="present.cclgmbi.diff.stress.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/observed/alpha.RD.AllClades.env_stabilityPastStress.pdf")

### Stability
env_stabilityPastStabil = alpha.div.master.pca.lg.diff %>% filter(model == "RD") %>%
  group_by(metric, clade, pool) %>%
  do({model = lm(obs.z~present.cclgmbi.diff.stabil.score, data=.)    # create your model
  data.frame(tidy(model),              # get coefficient info
             glance(model))})

env_stabilityPastStabil.df <- as.data.frame(env_stabilityPastStabil %>% dplyr::select(metric, clade, pool, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
write.csv(env_stabilityPastStabil.df, file="output/9_Environment/alpha/regression/observed/env_stabilityPastStabil.csv")

plotAlphaEnv(df=as.data.frame(alpha.div.master.pca.lg.diff), indepVar="present.cclgmbi.diff.stabil.score", xx = 0, dot=T, eqn=T, filename="output/9_Environment/alpha/regression/observed/alpha.RD.AllClades.env_stabilityPastStabil.pdf")






