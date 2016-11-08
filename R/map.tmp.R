


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
  geom_point(data = alps.env.sprich.summits, aes(x=X_WGS84, y=Y_WGS84, colour = Row.names), size=6, pch=17) +
  scale_colour_manual(name = "Row.names",values = cols) +
  scale_alpha_discrete(range=c(1,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('long') + ylab('lat') + 
  theme(legend.position = "none") +
  scaleBar(lon = 6, lat = 44.6, distanceLon = 5, distanceLat = 1, distanceLegend = 2, dist.unit = "km", 
            arrow.length = 4, arrow.distance = 3, arrow.North.size = 4)

# et voila!
ggsave("output/7_Maps/EcrinsSummitsColor.pdf", width = 8, height = 8)




