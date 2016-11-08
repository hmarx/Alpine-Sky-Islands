
library(ecodist)

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


## Beta diversity matrix
phylosor.alps <- read.csv("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/18Jan2016/distance_SES/SES_PhyloSor.csv", row.names=1)
phylosor.alps.summits <- as.dist(phylosor.alps[-c(1:4), -c(1:4)]) #convert to distance matric
ordering <- labels(phylosor.alps.summits) #make sure distance matrices are in the same order for mantel, etc.


### Species netural: avaialble summit area (~ Area of S. facing slope)

### Species netural: distance between summits 


### Available energy:
# BIO1 = Annual Mean Temperature
# BIO8 = Mean Temperature of Wettest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO16 = Precipitation of Wettest Quarter
pc.availener <- princomp(na.omit(pezAlpes$env[c("present_bio_1", "present_bio_8", "present_bio_12", "present_bio_13", "present_bio_16")]))
summary(pc.availener)
print(pc.availener$loadings, cutoff=0.001)
pc.availener.score <- scores(pc.availener)[,1]

edispc.availener.score <- vegdist(pc.availener.score, method = "euclid")

edispc.availener.score.mat <- as.matrix(edispc.availener.score)[ordering, ordering]
edispc.availener.score <- as.dist(edispc.availener.score.mat)

plot(phylosor.alps.summits ~ edispc.availener.score)
mantel(phylosor.alps.summits ~ edispc.availener.score)
MRM(phylosor.alps.summits ~ edispc.availener.score)


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

edispc.pc.stress.score <- vegdist(pc.stress.score, method = "euclid")
#ordering <- sort(attr(edispc.pc.stress.score, "Labels"))
edispc.pc.stress.score.mat <- as.matrix(edispc.pc.stress.score)[ordering, ordering]
edispc.pc.stress.score <- as.dist(edispc.pc.stress.score.mat)

mantel(phylosor.alps.summits ~ edispc.pc.stress.score)
MRM(phylosor.alps.summits ~ edispc.pc.stress.score)


### Environmental favourableness:
# BIO5 = Max Temperature of Warmest Month
# BIO10 = Mean Temperature of Warmest Quarter 
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter



### Environmental stability:
# BIO2 = Mean Diurnal Range
# BIO3 = Isothermality (BIO2/BIO7)
# BIO4 = Temperature Seasonality
# BIO15 = Precipitation Seasonality (Coefficient of Variation)

# presnet_bio - cclgmbi
# presnet_bio - mrlgmbi

pc.stabil <- prcomp(na.omit(pezAlpes$env[c("present_bio_2", "present_bio_3", "present_bio_4", "present_bio_15")]))
pc.stabil.score <- scores(pc.stabil)

edispc.pc.stabil.score <- vegdist(pc.stabil.score, method = "euclid")
mantel(as.dist(phylosor.alps[-c(1:4), -c(1:4)]) ~ edispc.pc.stabil.score)


  


### Environmental heterogeneity: 
# range_elev
# range_slope

(pezAlpes$env[c("range_elev", "range_slope")])

pc.hetero <- princomp(na.omit(pezAlpes$env[c("max_elev", "max_elev")]))
summary(pc.hetero)
pc.hetero.score <- scores(pc.hetero)

edis.pc.hetero.score <- vegdist(pc.hetero.score, method = "euclid")
mantel(as.dist(phylosor.alps[-c(1:4), -c(1:4)]) ~ edis.pc.hetero.score)




