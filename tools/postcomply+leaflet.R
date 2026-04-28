# 1. Load necessary libraries
library(akima)
library(sf)
library(leafpop)
library(leaflet)
library(leaflet.extras)
library(maps)
library(mapdata)
library(webshot)
library(lattice)
library(htmlwidgets)
library(pandoc)
library(htmltools)
#library(ncdf4)
library(raster)



lookup <- c(N=0,NNE=22.5,NE=45,ENE=67.5,E=90,ESE=112.5,SE=135,SSE=157.5,S=180,SSW=202.5,SW=225,WSW=247.5,W=270,WNW=292.5,NW=315,NNW=337.5)


## Read COMPLY input and output files

## Get source location from input file
a <- "/work/ROMO/users/kpc/comply/dev/v20/input.level4.txt"
lat <- read.csv(a, skip = 8, nrows = 1,header=F, sep="|")
a <- "/work/ROMO/users/kpc/comply/dev/v20/input.level4.txt"
long <- read.csv(a, skip = 9, nrows = 1,header=F, sep="|")
both <- as.data.frame(t(rbind(lat$V2,long$V2)))
colnames(both) <- c("Latitude","Longitude")
center <- both

xlat <- as.vector(lat$V2)
xlong <- as.vector(long$V2)
xlat <- 35.887591
xlong <- -78.84179


## Get key receptor locations from the input file
a <- "/work/ROMO/users/kpc/comply/dev/v20/input.level4.txt"
dist <- read.csv(a, skip = 19, nrows = 3,header=F, sep="|")
dirs <- read.csv(a, skip = 22, nrows = 3,header=F, sep="|")
both <- as.data.frame(t(rbind(dirs$V2,dist$V2)))
colnames(both) <- c("dirs","dist")
both$dirs <- trimws(both$dirs)
both$angle <- lookup[both$dirs]
both$dist <- as.numeric(as.character(both$dist))
both$x <- both$dist * sin(both$angle * 3.1416 / 180) # *pi/180 to convert degrees to radiens
both$y <- both$dist * cos(both$angle * 3.1416 / 180)

type <- c("Veg","Milk","Meat")
for(i in 1:nrow(both)) {
  both$new_var[i] <- type[i]
}


# convert to lat/long
temp <- both[c("x","y")]
crs.old <- paste("+proj=lcc +lat_1=",xlat," +lat_2=",xlat," +lat_0=",xlat," +lon_0=",xlong," +datum=NAD83",sep="")
crs.new <- c("+proj=longlat +datum=NAD83")
g <- st_as_sf(temp,coords=c("x","y"),crs = crs.old)
p <- st_transform(g, crs.new)
g2 <- matrix(unlist(p),ncol=2,byrow=TRUE)
colnames(g2) <- c("Longitude","Latitude")
keyreceptors <- cbind(both,g2)


## Get the MEI receptor for each direction from the input file
a <- "/work/ROMO/users/kpc/comply/dev/v20/input.level4.txt"
temp <- read.csv(a, skip = 36, nrows = 16,header=F, sep="|")
temp2 <- unlist(strsplit(temp$V2,","))
dist <- c(temp2[4],temp2[9],temp2[14],temp2[19],temp2[24],temp2[29],temp2[34],temp2[39],temp2[44],temp2[49],temp2[54],temp2[59],temp2[64],temp2[69],temp2[74],temp2[79])
dist
dirs <- c(temp2[1],temp2[6],temp2[11],temp2[16],temp2[21],temp2[26],temp2[31],temp2[36],temp2[41],temp2[46],temp2[51],temp2[56],temp2[61],temp2[66],temp2[71],temp2[76])
dirs

both <- as.data.frame(t(rbind(dirs,dist)))
both$dirs <- trimws(both$dirs)
both$angle <- lookup[both$dirs]
both$dist <- as.numeric(as.character(both$dist))
both$x <- both$dist * sin(both$angle * 3.1416 / 180) # *pi/180 to convert degrees to radiens
both$y <- both$dist * cos(both$angle * 3.1416 / 180)
mei <- both

# convert to lat/long
temp <- mei[c("x","y")]
crs.old <- paste("+proj=lcc +lat_1=",xlat," +lat_2=",xlat," +lat_0=",xlat," +lon_0=",xlong," +datum=NAD83",sep="")
crs.new <- c("+proj=longlat +datum=NAD83")
g <- st_as_sf(temp,coords=c("x","y"),crs = crs.old)
p <- st_transform(g, crs.new)
g2 <- matrix(unlist(p),ncol=2,byrow=TRUE)
colnames(g2) <- c("Longitude","Latitude")
xmei <- cbind(mei,g2)


## Read output diagnostic file
a <- "/work/ROMO/users/kpc/comply/dev/v20/diag.level4.txt"
comply <- read.csv(file=a,skip=0,header=T)
comply$SECTOR <- trimws(comply$SECTOR)
comply$RADNAME <- trimws(comply$RADNAME)
comply$DOSE <- comply$IMMERSION + comply$INHALATION + comply$SURFACE + comply$MEAT + comply$MILK + comply$VEG
head(comply)

one <- subset(comply,STACK==1 & RADNAME=="I-131")
#lookup <- c(N=0,NNE=22.5,NE=45,ENE=67.5,E=90,ESE=112.5,SE=135,SSE=157.5,S=180,SSW=202.5,SW=225,WSW=247.5,W=270,WNW=292.5,NW=315,NNW=337.5)
one$angle <- lookup[one$SECTOR]
one$x <- one$DISTANCE * sin(one$angle * 3.1416 / 180) # *pi/180 to convert degrees to radiens
one$y <- one$DISTANCE * cos(one$angle * 3.1416 / 180)
head(one)

# Create a new data frame of the dose data by downwind receptor location
x_data <- one$x
y_data <- one$y
z_data <- one$DOSE
dat <- data.frame(x = one$x, y = one$y, z = z_data, dist = one$DISTANCE)

# convert to lat/long
temp <- dat[c("x","y")]
crs.old <- paste("+proj=lcc +lat_1=",xlat," +lat_2=",xlat," +lat_0=",xlat," +lon_0=",xlong," +datum=NAD83",sep="")
crs.new <- c("+proj=longlat +datum=NAD83")
g <- st_as_sf(temp,coords=c("x","y"),crs = crs.old)
p <- st_transform(g, crs.new)
g2 <- matrix(unlist(p),ncol=2,byrow=TRUE)
colnames(g2) <- c("Longitude","Latitude")
dose <- cbind(dat,g2)

# Define the new grid points for the interpolated receptors
x_out <- seq(-1000, 1000, by=5)
y_out <- x_out

# Perform bilinear interpolation
# The 'linear = TRUE' argument specifies bilinear interpolation
interpolated_data <- with(dat, akima::interp(x, y, z, xo = x_out, yo = y_out, linear = TRUE, extrap = TRUE))
# The result is a list with components x, y, and a matrix z
interpolated_data$z[interpolated_data$z > 10] <- 9.99

# Convert to lat/long for leaflet plot
xtemp <- as.data.frame(interpolated_data[c("x","y")])
crs.old <- paste("+proj=lcc +lat_1=",xlat," +lat_2=",xlat," +lat_0=",xlat," +lon_0=",xlong," +datum=NAD83",sep="")
crs.new <- c("+proj=longlat +datum=NAD83")
g <- st_as_sf(xtemp,coords=c("x","y"),crs = crs.old)
p <- st_transform(g, crs.new)
g2 <- matrix(unlist(p),ncol=2,byrow=TRUE)
colnames(g2) <- c("Longitude","Latitude")
g3 <- as.data.frame(g2)
#xdose <- cbind(interpolated_data$z,g2)

### Leaflet process grid to reprojected raster

rlon     <- range(g3$Longitude)
rlat     <- range(g3$Latitude)
xy       <- cbind(as.vector(g3$Longitude), as.vector(g3$Latitude))

r <- raster(ncols=401, nrows=401, xmn=rlon[1], xmx=rlon[2], ymn=rlat[1], ymx=rlat[2])
r1 <- rasterize(xy, r)
values(r1) <- as.vector(interpolated_data$z)


### PLOT
require(fields)

my.colors <- colorRampPalette(c("white","purple","deepskyblue2","green","yellow","orange","red","brown"))
ncols <- 50

##### LEAFLET

leaflet_map <- c("OpenStreetMap.Mapnik","OpenTopoMap","Esri.WorldImagery","Esri.WorldStreetMap","USGS.USTopo","USGS.USImagery","USGS.USImageryTopo")

GetURL <- function(service, host = "basemap.nationalmap.gov") {
  sprintf("https://%s/arcgis/services/%s/MapServer/WmsServer", host, service)
}

att <- paste0("<a href='https://www.usgs.gov/'>",
              "U.S. Geological Survey</a> | ",
              "<a href='https://www.usgs.gov/laws/policies_notices.html'>",
              "Policies</a>")

mapStates <- map("state",fill=T,plot=F)
my.leaf.base <- leaflet(options = leafletOptions(zoomSnap = 0.25, zoomDelta =0.25), data=mapStates) %>%
  addProviderTiles(leaflet_map[1],group="Street Map")  %>%
  addProviderTiles(leaflet_map[2],group="Topo Map") %>%
  addProviderTiles(leaflet_map[3],group="ESRI World Imagery") %>%
  addProviderTiles(leaflet_map[4],group="ESRI Street Map") %>%
  addWMSTiles(GetURL("USGSImageryOnly"),group="USGS Imagery",layers="0",attribution=att) %>%
  addWMSTiles(GetURL("USGSTopo"),group="USGS Topo",layers="0",attribution=att) %>%
  addWMSTiles(GetURL("USGSImageryTopo"),group="USGS Topo Imagery",layers="0",attribution=att)

base_Groups = c("Street Map","Topo Map","ESRI World Imagery","ESRI Street Map")


### LEAF PREP

sinfo_data <- data.frame(stat_id="NA",lat=dose$Latitude,lon=dose$Longitude,Dose=dose$z,Distance=dose$dist)
sinfo_mei_data <- data.frame(stat_id="MEI receptor",lat=xmei$Latitude,lon=xmei$Longitude,Distance=xmei$dist)
sinfo_key_data <- data.frame(stat_id=keyreceptors$new_var,lat=keyreceptors$Latitude,lon=keyreceptors$Longitude,Distance=keyreceptors$dist)

units <- 'mrem/year'
contents_data <- paste("Prediction:", round(sinfo_data$Dose, 1), units,
                       "at distance:",sinfo_data$Distance, sep =" ")

contents_mei_data <- paste("Receptor:", mei_data$stat_id,
                       "distance:", round(mei_data$Distance, 1),"m", sep=" ")
contents_key_data <- paste("Receptor:", key_data$stat_id,
                           "distance:", round(key_data$Distance, 1),"m", sep=" ")

sinfo_data$Dose[sinfo_data$Dose > 10] <- 9.99
omax <- max(sinfo_data$Dose, na.rm=T)
obins<- seq(0.0,10,by=0.5)
obins <- c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,6,7,8,9,10)
my.colors <- colorRampPalette(c(grey(.8),"mediumpurple","darkorchid4", "#002FFF","lightblue", "green", "yellow", "orange", "red", "brown"))
binpal2 <- colorBin(my.colors(length(obins)), c(min(obins),max(obins)), length(obins)-1 , pretty = FALSE)
pal1 <- colorBin(my.colors(length(obins)), c(min(obins),max(obins)), na.color = "transparent",pretty=T,bins=obins)

### LEAFLET

fill_opacity <- 0.66
transp <- 0.66
plot_rad <- 5
plot_names <- c("DOSE","INTERP","MEI","KEYREC")
legend_title <- c("DOSE","INTERP","MEI","KEYREC")
main_title <- paste("COMPLY v2.0",sep="/")
base_Groups = c("ESRI Street Map","ESRI World Imagery","Topo Map")

my.leaf.base <- leaflet(options = leafletOptions(zoomSnap = 0.25, zoomDelta =0.25), data=mapStates)
my.leaf <- my.leaf.base
my.leaf <- my.leaf %>% addProviderTiles("Esri.WorldStreetMap",group="ESRI Street Map")
my.leaf <- my.leaf %>% addProviderTiles("Esri.WorldImagery",group="ESRI World Imagery")
my.leaf <- my.leaf %>% addProviderTiles("OpenTopoMap",group="Topo Map")

my.leaf <- my.leaf %>% addCircleMarkers(sinfo_data$lon,sinfo_data$lat,color="black",fillColor=~pal1(sinfo_data$Dose),group="DOSE",radius=plot_rad,data=sinfo_data,opacity=1,fillOpacity=fill_opacity,stroke=TRUE,weight=1,popup=contents_data, label=contents_data, labelOptions = labelOptions(noHide = F, textsize = "15px"))
my.leaf <- my.leaf %>% addCircleMarkers(sinfo_mei_data$lon,sinfo_mei_data$lat,color="black",fillColor="yellow",group="MEI",radius=plot_rad*1.66,data=sinfo_mei_data,opacity=0.66,fillOpacity=fill_opacity,stroke=TRUE,weight=1,popup=contents_mei_data, label=contents_mei_data, labelOptions = labelOptions(noHide = F, textsize = "15px"))
my.leaf <- my.leaf %>% addCircleMarkers(sinfo_key_data$lon,sinfo_key_data$lat,color="black",fillColor="red",group="KEYREC",radius=plot_rad*1.66,data=sinfo_key_data,opacity=0.85,fillOpacity=fill_opacity,stroke=TRUE,weight=1,popup=contents_key_data, label=contents_key_data, labelOptions = labelOptions(noHide = F, textsize = "15px"))

my.leaf <- my.leaf %>% addRasterImage(flip(r1,"y"), colors = pal1, opacity = transp, group="INTERP")

my.leaf <- my.leaf %>% leaflet::addLegend(position="bottomright", pal=pal1,values=obins,group="DOSE",layerId="DOSE",title="DOSE (mrem/yr)",opacity=1)
my.leaf <- my.leaf %>% leaflet::addLegend(position="bottomright", pal=pal1,values=obins,group="INTERP",layerId="INTERP",title="DOSE (mrem/yr)",opacity=1)

my.leaf <- my.leaf %>% addControl(main_title,position="topleft",className="map-title")
my.leaf <- my.leaf %>% addLayersControl(baseGroups = base_Groups,overlayGroups = plot_names, options = layersControlOptions(collapsed = FALSE,position="topleft"))
my.leaf <- my.leaf %>% hideGroup(c("DOSE","MEI","KEYREC"))

saveWidget(my.leaf, file='/home/kbaker/comply/comply.html',selfcontained=T)

