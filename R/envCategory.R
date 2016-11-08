
### Multiple regression on beta diversity distance 

##### http://pubs.usgs.gov/ds/691/ds691.pdf
# BIO1 = Annual Mean Temperature : The annual mean temperature approximates the total energy inputs for an ecosystem.
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)) : This index can help provide information pertaining to the relevance of temperature fluctuation for different species
# BIO3 = Isothermality (BIO2/BIO7) (* 100) : is a quantification of how large the day-to-night temperature oscillation is in comparison to the summer-to-winter oscillation. A value of 100 would represent a site where the diurnal temperature range is equal to the annual temperature range. A value of 50 would indicate a location where the diurnal temperature range is half of the annual temperature range.
# BIO4 = Temperature Seasonality (standard deviation *100) :Temperature seasonality is a measure of temperature change over the course of the year. 
# BIO5 = Max Temperature of Warmest Month : This information is useful when examining whether species distributions are affected by warm temperature anomalies throughout the year
# BIO6 = Min Temperature of Coldest Month : This information is useful when examining whether species distributions are affected by cold temperature anomalies throughout the year
# BIO7 = Temperature Annual Range (BIO5-BIO6) : whether species distributions are affected by ranges of extreme temperature conditions
# BIO8 = Mean Temperature of Wettest Quarter : This index provides mean temperatures during the wettest three months of the year, which can be useful for examining how such environmental fac- tors may affect species seasonal distributions
# BIO9 = Mean Temperature of Driest Quarter: This index provides mean tem- peratures during the driest three months of the year, which can be useful for examining how such environmental factors may affect species seasonal distributions
# BIO10 = Mean Temperature of Warmest Quarter : This index provides mean tem- peratures during the warmest three months of the year, which can be useful for examining how such environmental factors may affect species seasonal distributions
# BIO11 = Mean Temperature of Coldest Quarter : This index provides mean temperatures during the coldest three months of the year, which can be useful for examining how such environmental fac- tors may affect species seasonal distributions
# BIO12 = Annual Precipitation : Annual total precipitation approximates the total water inputs and is therefore useful when ascertaining the importance of water availability to a species distribution
# BIO13 = Precipitation of Wettest Month : The wettest month is useful if extreme precipitation conditions during the year influence a species potential range
# BIO14 = Precipitation of Driest Month : The driest month is useful if extreme precipitation conditions during the year influence a species potential range
# BIO15 = Precipitation Seasonality (Coefficient of Variation) : Since species distributionscanbe strongly influenced by variability in precipitation, this index provides a percentage of precipitation variability where larger percentages represent greater variability of precipitation. 
# BIO16 = Precipitation of Wettest Quarter : This index provides total precipitation during the wettest three months of the year, which can be useful for examining how such environmental fac- tors may affect species seasonal distributions.
# BIO17 = Precipitation of Driest Quarter : This index provides total precipitation during the driest three months of the year, which can be useful for examining how such environmental fac- tors may affect species seasonal distributions
# BIO18 = Precipitation of Warmest Quarter : This index provides total precipitation during the warmest three months of the year, which can be useful for examining how such environmental factors may affect species seasonal distributions
# BIO19 = Precipitation of Coldest Quarter : his index provides total precipitation during the coldest three months of the year, which can be useful for examining how such environmental fac- tors may affect species seasonal distributions

names(pezAlpes$env)
### PCA of all present:
pca.presnent <- prcomp(na.omit(pezAlpes$env[c("present_bio_1", "present_bio_2", "present_bio_3", "present_bio_4", 
                                                "present_bio_5", "present_bio_6", "present_bio_7", "present_bio_8", 
                                                "present_bio_9", "present_bio_10", "present_bio_11", "present_bio_12", 
                                                "present_bio_13", "present_bio_14", "present_bio_15", "present_bio_16", 
                                                "present_bio_17", "present_bio_18", "present_bio_19")]))
summary(pca.presnent)
pca.presnent.score <- scores(pca.presnent)[,1]
pca.presnent.score.df <- as.data.frame(pca.presnent.score)

### PCA of all LGM:
pca.cclgmbi <- prcomp(na.omit(pezAlpes$env[c("cclgmbi1", "cclgmbi2", "cclgmbi3", "cclgmbi4", 
                                              "cclgmbi5", "cclgmbi6", "cclgmbi7", "cclgmbi8", 
                                              "cclgmbi9", "cclgmbi10", "cclgmbi11", "cclgmbi12", 
                                              "cclgmbi13", "cclgmbi14", "cclgmbi15", "cclgmbi16", 
                                              "cclgmbi17", "cclgmbi18", "cclgmbi19")]))
summary(pca.cclgmbi)
pca.cclgmbi.score <- scores(pca.cclgmbi)[,1]
pca.cclgmbi.score.df <- as.data.frame(pca.cclgmbi.score)

### PCA of all LGM:
pca.melgmbi <- prcomp(na.omit(pezAlpes$env[c("melgmbi1", "melgmbi2", "melgmbi3", "melgmbi4", 
                                             "melgmbi5", "melgmbi6", "melgmbi7", "melgmbi8", 
                                             "melgmbi9", "melgmbi10", "melgmbi11", "melgmbi12", 
                                             "melgmbi13", "melgmbi14", "melgmbi15", "melgmbi16", 
                                             "melgmbi17", "melgmbi18", "melgmbi19")]))
summary(pca.melgmbi)
pca.melgmbi.score <- scores(pca.melgmbi)[,1]
pca.melgmbi.score.df <- as.data.frame(pca.melgmbi.score)
pca.lg <- cbind(pca.presnent.score.df, pca.cclgmbi.score.df, pca.melgmbi.score.df)


### Species netural: avaialble summit area (~ Area of S. facing slope)


### Species netural: distance between summits 
spatial.dist <- vegdist(na.omit(pezAlpes$env[c("X_WGS84", "Y_WGS84")]), method = "euclid")

### Available energy : 
# BIO1 = Annual Mean Temperature
# BIO5 = Max Temperature of Warmest Month
# BIO8 = Mean Temperature of Wettest Quarter
# BIO10 = Mean Temperature of Warmest Quarter 
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO16 = Precipitation of Wettest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter
pc.availener <- princomp(na.omit(pezAlpes$env[c("present_bio_1", "present_bio_5", "present_bio_8", "present_bio_10", "present_bio_12", "present_bio_13",
                                                "present_bio_16", "present_bio_18", "present_bio_19")]))
summary(pc.availener)
print(pc.availener$loadings, cutoff=0.001)
pc.availener.score <- scores(pc.availener)[,1]
pc.availener.score.df <- as.data.frame(pc.availener.score)

### Environmental stress: 
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO9 = Mean Temperature of Driest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO14 = Precipitation of Driest Month
# BIO17 = Precipitation of Driest Quarter
pc.stress <- princomp(na.omit(pezAlpes$env[c("present_bio_6", "present_bio_7", "present_bio_9", "present_bio_11", "present_bio_14", "present_bio_17")]))
summary(pc.stress)
print(pc.stress$loadings, cutoff=0.001)
print(pc.stress$scores, cutoff=0.001)
pc.stress.score <- pc.stress$scores[,1]
pc.stress.score.df <- as.data.frame(pc.stress.score)


### Present Environmental stability:
# BIO2 = Mean Diurnal Range
# BIO3 = Isothermality (BIO2/BIO7)
# BIO4 = Temperature Seasonality
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
pc.stabil <- princomp(na.omit(pezAlpes$env[c("present_bio_2", "present_bio_3", "present_bio_4", "present_bio_15")]))
summary(pc.stabil)
print(pc.stabil$loadings, cutoff=0.001)
print(pc.stabil$scores, cutoff=0.001)
pc.stabil.score <- scores(pc.stabil)[,1]
pc.stabil.score.df <- as.data.frame(pc.stabil.score)


### Environmental heterogeneity: 
# range_elev
# range_slope
pc.hetero <- princomp(na.omit(pezAlpes$env[c("range_elev", "range_slope")]))
summary(pc.hetero)
print(pc.hetero$loadings, cutoff=0.001)
print(pc.hetero$scores, cutoff=0.001)
pc.hetero.score <- scores(pc.hetero)[,1]
pc.hetero.score.df <- as.data.frame(pc.hetero.score)

env.pca <- cbind(pc.availener.score.df, pc.stress.score.df, pc.stabil.score.df, pc.hetero.score.df) #pc.favor.score



### Environmental stability LGM to present:
# presnet_bio - cclgmbi
presnent <- as.matrix(na.omit(pezAlpes$env[c("present_bio_1", "present_bio_2", "present_bio_3", "present_bio_4", 
                                             "present_bio_5", "present_bio_6", "present_bio_7", "present_bio_8", 
                                             "present_bio_9", "present_bio_10", "present_bio_11", "present_bio_12", 
                                             "present_bio_13", "present_bio_14", "present_bio_15", "present_bio_16", 
                                             "present_bio_17", "present_bio_18", "present_bio_19")]))
head(presnent)

cclgmbi <- as.matrix(na.omit(pezAlpes$env[c("cclgmbi1", "cclgmbi2", "cclgmbi3", "cclgmbi4", 
                                            "cclgmbi5", "cclgmbi6", "cclgmbi7", "cclgmbi8", 
                                            "cclgmbi9", "cclgmbi10", "cclgmbi11", "cclgmbi12", 
                                            "cclgmbi13", "cclgmbi14", "cclgmbi15", "cclgmbi16", 
                                            "cclgmbi17", "cclgmbi18", "cclgmbi19")]))

head(cclgmbi)

present.cclgmbi.diff <- as.data.frame(presnent - cclgmbi)

present.cclgmbi.diff.availener <- princomp(present.cclgmbi.diff[c("present_bio_5", "present_bio_10", "present_bio_18", "present_bio_19", 
                                                                  "present_bio_1", "present_bio_8", "present_bio_12", "present_bio_13", "present_bio_16")])
present.cclgmbi.diff.availener
summary(present.cclgmbi.diff.availener)
print(present.cclgmbi.diff.availener$loadings, cutoff=0.001)
present.cclgmbi.diff.availener.score <- scores(present.cclgmbi.diff.availener)[,1]
present.cclgmbi.diff.availener.score.df <- as.data.frame(present.cclgmbi.diff.availener.score)

present.cclgmbi.diff.stress <- princomp(na.omit(pezAlpes$env[c("present_bio_6", "present_bio_7", "present_bio_9", "present_bio_11", "present_bio_14", "present_bio_17")]))
summary(present.cclgmbi.diff.stress)
print(present.cclgmbi.diff.stress$loadings, cutoff=0.001)
print(present.cclgmbi.diff.stress$scores, cutoff=0.001)
present.cclgmbi.diff.stress.score <- present.cclgmbi.diff.stress$scores[,1]
present.cclgmbi.diff.stress.score.df <- as.data.frame(present.cclgmbi.diff.stress.score)

present.cclgmbi.diff.stabil <- princomp(na.omit(pezAlpes$env[c("present_bio_2", "present_bio_3", "present_bio_4", "present_bio_15")]))
summary(present.cclgmbi.diff.stabil)
print(present.cclgmbi.diff.stabil$loadings, cutoff=0.001)
print(present.cclgmbi.diff.stabil$scores, cutoff=0.001)
present.cclgmbi.diff.stabil.score <- scores(present.cclgmbi.diff.stabil)[,1]
present.cclgmbi.diff.stabil.score.df <- as.data.frame(present.cclgmbi.diff.stabil.score)

present.cclgmbi.diff.hetero <- princomp(na.omit(pezAlpes$env[c("range_elev", "range_slope")]))
summary(present.cclgmbi.diff.hetero)
print(present.cclgmbi.diff.hetero$loadings, cutoff=0.001)
print(present.cclgmbi.diff.hetero$scores, cutoff=0.001)
present.cclgmbi.diff.hetero.score <- scores(present.cclgmbi.diff.hetero)[,1]
present.cclgmbi.diff.hetero.score.df <- as.data.frame(present.cclgmbi.diff.hetero.score)

env.pca.past.present <- cbind(present.cclgmbi.diff.availener.score.df, present.cclgmbi.diff.stress.score.df, 
                              present.cclgmbi.diff.stabil.score.df, present.cclgmbi.diff.hetero.score.df) 




