#http://www.nickeubank.com/wp-content/uploads/2015/10/RGIS3_MakingMaps_part1_mappingVectorData.html

source("R/addScaleBarMap.R")

### Found NP layer on: https://inpn.mnhn.fr/telechargement/cartes-et-information-geographique/ep/pn?lg=en
## edited in qgis to include only Ecrins -> saved layer 
layerName <- "ecrinsNP"  
data_projected <- readOGR(dsn="output/7_Maps/ecrinsNP/", layer=layerName) 
plot(data_projected) #+proj=longlat +datum=WGS84 +no_defs

###### Transform shapefile to dataframe ###### 
# add to data a new column termed "id" composed of the rownames of data
data_projected@data$id <- rownames(data_projected@data)

# reproject the data onto a "longlat" projection
subsetTransform <- spTransform(data_projected, CRS("+proj=longlat"))

# create a data.frame from our spatial object
ecrinsPoints <- fortify(subsetTransform, region = "id")
head(ecrinsPoints)

# merge the "fortified" data with the data from our spatial object
ectinsDF <- merge(ecrinsPoints, subsetTransform@data, by = "id")
head(ectinsDF)
unique(ectinsDF$CODE_R_ENP)

# Core national park polygon
coreNP <- ectinsDF[ectinsDF$CODE_R_ENP == "CPN",]


################################## Map of all France with National 
map("worldHires","France", col="gray90", fill=TRUE)

france <- map_data("france")
head(france)
unique(france$region)
states <- subset(france, region %in% c("Hautes-Alpes", "Isere"))

names(alps.env.sprich)

cols = c("Brevoort" = "red4",
            "Pelvoux" = "red2",
            "Choisy" = "tomato2",
            "Occidentale Ailefroide" = "darkorange1",
            "Plat de la Selle" = "goldenrod",
            "Burlan"= "yellow",
            "La Meije"= "chartreuse",
            "Barre des Ecrins"= "limegreen",
            "Rouies"= "darkgreen",
            "Muraillette"= "darkslategrey",
            "Olan"= "darkturquoise", 
            "Sirac"= "dodgerblue",
            "Lauvitel"= "blue",    
            "Rateau"= "purple2", 
            "Rocher de la Selle"= "hotpink")



pdf(file="output/7_Maps/FranceMapSummitsColor.pdf") 
p <- ggplot()
#p <- p + geom_polygon(data=europe, aes(x=long, y=lat, group=group), colour="grey", fill="grey")
p <- p + geom_polygon(data=france, aes(x=long, y=lat, group=group), colour="grey", fill="grey")
#p <- p + geom_polygon(data=states, aes(x=long, y=lat, group=group), colour="snow4", fill="snow4")
p <- p + geom_polygon(data=coreNP, aes(x=long, y=lat, group=group), colour="grey30", fill="snow4", lty=1, lwd=.1) 
p <- p + geom_point(data = alps.env.sprich.summits, aes(x=X_WGS84, y=Y_WGS84, colour = Row.names), size=.5)
p <- p + scale_colour_manual(name = "Row.names",values = cols)
p <- p + geom_point(aes(x=2.3508, y=48.8567), color="black")
p <- p + geom_point(aes(x=5.7222, y=45.2002), color="black")
p <- p + annotate("text", x=5.7222, y=45.5, color="black", label="Grenoble", size= 7)
p <- p + annotate("text", x=2.3508, y=49.25, color="black", label="Paris", size= 7)
#p <- p + theme(legend.position = "none", axis.title = element_blank())
p <- p + scaleBar(lon = -5, lat = 42.5, distanceLon = 100, distanceLat = 15, distanceLegend = 30, dist.unit = "km", 
                    arrow.length = 80, arrow.distance = 50, arrow.North.size = 4)
p <- p + theme(panel.grid.minor = element_line(colour = NA),
               panel.background = element_rect(fill = NA, colour = "black"), legend.position = "none")
p
dev.off()


##################   Map of Ecrins NP zoomed in
########### Import custom raster basemap (made in qgis)
### Baselayer from:
#http://felix.rohrba.ch/en/2016/awesome-basemap-layer-for-your-qgis-project/
#http://www.naturalearthdata.com/features/

layer <- stack("output/7_Maps/EcrinsBaseZoom.tif") #EPSG: 3857
projection(layer) <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"

# Reproject layer to lat long
sr <-"+proj=longlat"
projected_raster <- projectRaster(layer, crs = sr)
plot(projected_raster$EcrinsBaseZoom.1, col= grey(1:99/100))

## bounding box
domain <- c(5.9, 6.6, 44.55,  45.1)

## crop layer to bounding box
ecrins.base.crop <- crop(projected_raster, y=extent(domain))
plot(ecrins.base.crop$EcrinsBaseZoom.1, col= grey(1:99/100))

## create raster dataframe
rast.table <- data.frame(xyFromCell(ecrins.base.crop, 1:ncell(ecrins.base.crop)), getValues(ecrins.base.crop/255))
head(rast.table)

## plot
ggplot(data = rast.table, aes(x = x, y = y)) +
  geom_raster(aes(fill=EcrinsBaseZoom.1)) +
  scale_fill_gradientn(colours=c("grey61","grey100")) +
  geom_polygon(data=coreNP, aes(x=long, y=lat, group=group), colour="white", fill="grey10", alpha=0.4) + 
  geom_point(data = alps.env.sprich.summits, aes(x=X_WGS84, y=Y_WGS84, colour = Row.names, size = elevation), pch=17) + #size=6, 
  scale_colour_manual(name = "Row.names",values = cols) +
  scale_alpha_discrete(range=c(1,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('long') + ylab('lat') + 
  theme(legend.position = "none") +
  scaleBar(lon = 5.94, lat = 44.6, distanceLon = 5, distanceLat = 1, distanceLegend = 2, dist.unit = "km", 
           arrow.length = 4, arrow.distance = 3, arrow.North.size = 4) +
  scale_size_continuous(range = c(3, 8))

# et voila!
ggsave("output/7_Maps/EcrinsSummitsColor.pdf", width = 8, height = 8)





### Other sites searched for layers
#http://www.eea.europa.eu/data-and-maps/data/digital-elevation-model-of-europe#tab-gis-data
#http://www.webgis.com/srtm3.html
#http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Eurasia/
#http://viewfinderpanoramas.org/dem3.html#alps

################# Plot SES metrics on summits

# get and plot a map
ecrins <- get_map(location = c(5.9, 44.7, 6.6, 45), source="google", maptype = "satellite", color="bw", zoom=10, crop=T)
#ecrins <- get_map(location = c(5.9, 44.7, 6.6, 45), maptype = "terrain", color="bw", zoom=10, crop=T)

#pdf(file="output/7_Maps/EcrinsSummitsColored.pdf")
ggmap(ecrins) + 
  geom_point(data = alps.env.sprich.summits, aes(x=X_WGS84, y=Y_WGS84, size=ntax, colour = Row.names)) +  
  geom_polygon(data=coreNP, aes(x=long, y=lat, group=group), colour="white", fill=NA, lty=3) +
  scale_colour_manual(name = "Row.names",values = cols) +
  scale_size_continuous(range = c(3, 10)) +
  ggtitle("SES mpd")
#dev.off()

#pdf(file="output/7_Maps/EcrinsSummitsColored.pdf")
ggmap(ecrins) + 
  geom_point(data = alps.env.sprich.summits, aes(x=X_WGS84, y=Y_WGS84, colour = Row.names, size = ntax)) +  
  geom_polygon(data=coreNP, aes(x=long, y=lat, group=group), colour="white", fill=NA, lty=3) +
  scale_colour_manual(name = "Row.names",values = cols) +
  scale_size_continuous(range = c(3, 5)) + 
  theme(legend.position = "none", axis.title = element_blank()) 
#dev.off()


dispersion.phylonull <- merge(ecrins.sesmpd.phylonull, ecrins.sesmntd.phylonull, by=0)
rownames(dispersion.phylonull) <- dispersion.phylonull$Row.names
env.dispersion <- merge(pezAlpes$env, na.omit(dispersion.phylonull[2:ncol(dispersion.phylonull)]), by=0, all.x=F)

#pdf(file="output/7_Maps/EcrinsSummitsTopoMNTD.pdf")
ggmap(ecrins) + 
  geom_point(data = env.dispersion, aes(x=X_WGS84, y=Y_WGS84, color=mntd.obs.z, size=ntaxa.y)) +  
  geom_polygon(data=coreNP, aes(x=long, y=lat, group=group), colour="white", fill=NA, lty=3) +
  scale_colour_gradientn(colours=rainbow(5)) +
  scale_size_continuous(range = c(3,10)) +
  ggtitle("SES mntd")
#dev.off()

#pdf(file="output/7_Maps/EcrinsSummitsTopoMPD.pdf")
ggmap(ecrins) + 
  geom_point(data = env.dispersion, aes(x=X_WGS84, y=Y_WGS84, color=mpd.obs.z, size=ntaxa.x)) +  
  geom_polygon(data=coreNP, aes(x=long, y=lat, group=group), colour="white", fill=NA, lty=3) +
  scale_colour_gradientn(colours=rainbow(5)) +
  scale_size_continuous(range = c(3, 10)) +
  ggtitle("SES mpd")
#dev.off()


