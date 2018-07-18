# OCP/IITA 2017 Maize fertilizer response trials, Central Nigeria
# Yield response data courtesy of IITA
# M. Walsh & J. Huising, July 2018

# Required packages
# install.packages(c("downloader","rgdal","raster","quantreg","leaflet","htmlwidgets")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(sp)
  require(raster)
  require(quantreg)
  require(leaflet)
  require(htmlwidgets)
})

# Create a data folder in your current working directory
dir.create("OCP_data", showWarnings=F)
setwd("./OCP_data")

# Data downloads -----------------------------------------------------------
# download IITA/OCP yield data
download("https://www.dropbox.com/s/efr02hlxn3n1yvn/OCP_trial_data.zip?raw=1", "OCP_trial_data.zip", mode = "wb")
unzip("OCP_trial_data.zip", overwrite = T)
sites <- read.table("location.csv", header=T, sep=",")
trial <- read.table("trial.csv", header=T, sep=",")
tresp <- merge(sites, trial, by="tid")

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
ext <- data.frame(lat = c(8.5,8.5,11.7,11.7), lon = c(5.3,11.8,5.3,11.8)) ## set ROI extent in degrees
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
colnames(tresp) <- c("state","lga","tid","lat","lon","alt","acc","team","trt","pdens","cdens","yc","yo")
# boxplot(yc~trt, tresp, notch=T)
# boxplot(yo~trt, tresp, notch=T)
# plot(yo~yc, tresp)

# project survey coords to grid CRS
tresp.proj <- as.data.frame(project(cbind(tresp$lon, tresp$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(tresp.proj) <- c("x","y")
tresp <- cbind(tresp, tresp.proj)
coordinates(tresp) <- ~x+y
projection(tresp) <- projection(tresp)

# extract gridded variables at survey locations
trespgrid <- extract(grids, tresp)
gsdat <- as.data.frame(cbind(tresp, trespgrid)) 
gsdat <- gsdat[complete.cases(gsdat[,c(9:11, 14:56)]),] ## removes incomplete cases

# Classify yield propensities by conditional quantile ---------------------
qy.rq <- rq(log(yo)~factor(trt)+pdens*cdens, tau = 0.5, data = gsdat) ## try quantiles other than the median
summary(qy.rq)
gsdat$qy <- as.factor(ifelse(exp(predict(qy.rq, gsdat)) > gsdat$yield, "B", "A"))
table(gsdat$qy)
table(gsdat$district, gsdat$qy)
boxplot(yield~qy, notch=T, gsdat)

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/OCP_gsdat.csv", row.names = F)

# Yield survey map widget -------------------------------------------------
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'OCP_response.html', selfcontained = T) ## save widget
