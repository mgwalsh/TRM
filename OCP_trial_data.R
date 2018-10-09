# OCP/IITA 2017 Maize fertilizer response trials, Central Nigeria
# Yield response data courtesy of IITA
# M. Walsh & J. Huising, July 2018

# Required packages
# install.packages(c("downloader","rgdal","raster","quantreg","arm","leaflet","htmlwidgets")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(sp)
  require(raster)
  require(quantreg)
  require(arm)
  require(leaflet)
  require(htmlwidgets)
})

# Create a data folder in your current working directory
dir.create("OCP_data", showWarnings=F)
setwd("./OCP_data")

# Data downloads -----------------------------------------------------------
# download IITA/OCP yield data
download("https://www.dropbox.com/s/hi75cnp3ejr4srk/OCP_trials.zip?raw=1", "OCP_trials.zip", mode = "wb")
unzip("OCP_trials.zip", overwrite = T)
sites <- read.table("sites.csv", header=T, sep=",")
trial <- read.table("trials.csv", header=T, sep=",")
tresp <- merge(sites, trial, by="sid")

# download GADM-L2 shapefile (courtesy: http://www.gadm.org)
download("https://www.dropbox.com/s/y3h6l7yu00orm78/NGA_adm2.zip?raw=1", "NGA_adm2.zip", mode = "wb")
unzip("NGA_adm2.zip", overwrite = T)
shape <- shapefile("NGA_adm2.shp")

# download raster stack (note this is a big 800+ Mb download)
download("https://www.dropbox.com/s/u5fyjbujf0d7q43/NG_250m_2017.zip?raw=1", "NG_250m_2017.zip", mode = "wb")
unzip("NG_250m_2017.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# set ROI grid extent
ext <- data.frame(lat = c(8.5,8.5,11.7,11.7), lon = c(5.3,11.8,5.3,11.8)) ## set ROI extent in decimal degrees
names(ext) <- c("lat","lon")
coordinates(ext) <- ~ lon + lat
proj4string(ext) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
ext <- spTransform(ext, CRS("+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
bb <- extent(ext)
grids <- crop(grids, bb)

# Data setup --------------------------------------------------------------
# attach GADM-L2 admin unit names from shape
coordinates(tresp) <- ~lon+lat
projection(tresp) <- projection(shape)
gadm <- tresp %over% shape
tresp <- as.data.frame(tresp)
tresp <- cbind(gadm[ ,c(5,7)], tresp)
colnames(tresp) <- c("state","lga","sid","lat","lon","alt","team","trt","ccob","tcob","twgt","cyld","tyld","ayld")

# project survey coords to grid CRS
tresp.proj <- as.data.frame(project(cbind(tresp$lon, tresp$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(tresp.proj) <- c("x","y")
tresp <- cbind(tresp, tresp.proj)
coordinates(tresp) <- ~x+y
projection(tresp) <- projection(tresp)

# extract gridded variables at survey locations
trespgrid <- extract(grids, tresp)
gsdat <- as.data.frame(cbind(tresp, trespgrid)) 
# plot(alt~MDEM, gsdat) ## gps altitude/location check against MDEM 

# Classify by site indices ------------------------------------------------
si.lmer <- lmer(log(tyld)~trt+(1|sid), gsdat) ## random intercept (site-level) model
display(si.lmer)
si.ran <- ranef(si.lmer) ## extract random effects
si <- as.data.frame(rownames(si.ran$sid))
si$SI <- si.ran$sid[,1]
colnames(si) <- c("sid","SI")
si$SIC <- ifelse(si$SI > 0, "A", "B") ## classify above/below average site indices
gsdat <- merge(gsdat, si, by="sid")

# Plots
boxplot(tyld~trt, notch=T, ylab="Cob yield (kg/ha)", ylim=c(0,8000), gsdat) ## treatment differences
boxplot(tyld~SIC, notch=T, gsdat) ## yield differences between site index groups
boxplot(tcob~trt*SIC, notch=T, ylab="Number of cobs", ylim=c(0,800), gsdat) ## treatment differences
boxplot(tyld~trt*SIC, notch=T, ylab="Cob yield (kg/ha)", ylim=c(0,8000), gsdat) ## treatment differences
plot(tyld~cyld, xlab="Cob yield (kg/ha), circular plot", ylab="Cob yield (kg/ha), total plot", gsdat)

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/OCP_tdat.csv", row.names = F)

# Yield survey map widget -------------------------------------------------
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'OCP_trials.html', selfcontained = T) ## save widget
