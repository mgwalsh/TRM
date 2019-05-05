# TAMASA maize fertilizer nutrient ommission trials, Central Ethiopia, 2015/16
# Yield trial data courtesy of CIMMYT
# M. Walsh, Z. Ahamed, J. Chamberlin & P. Craufurd, May 2019

# Required packages
# install.packages(c("downloader","rgdal","raster","arm","leaflet","htmlwidgets")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(sp)
  require(raster)
  require(arm)
  require(leaflet)
  require(htmlwidgets)
})

# Create a data folder in your current working directory
dir.create("TAMASA_data", showWarnings=F)
setwd("./TAMASA_data")

# Data downloads -----------------------------------------------------------
# download TAMASA yield data
download("https://www.dropbox.com/s/0cefqlb0g584mfd/ET_no_trials.zip?raw=1", "ET_no_trials.zip", mode = "wb")
unzip("ET_no_trials.zip", overwrite = T)
sites <- read.table("sites.csv", header=T, sep=",")
trial <- read.table("trials.csv", header=T, sep=",")
tresp <- merge(sites, trial, by="tid")

# download GADM-L3 shapefile (courtesy: http://www.gadm.org)
download("https://www.dropbox.com/s/25kw13359f8l5tr/ET_adm_shp.zip?raw=1", "ET_adm_shp.zip", mode = "wb")
unzip("ET_adm_shp.zip", overwrite = T)
shape <- shapefile("ETH_adm3.shp")

# download raster stack (note this is a big 1+ Gb download)
download("https://www.dropbox.com/s/iqix6sn66w04jo0/ET_250m_2018.zip?raw=1", "ET_250m_2018.zip", mode = "wb")
unzip("ET_250m_2018.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# Data setup --------------------------------------------------------------
# attach GADM-L3 admin unit names from shape
coordinates(tresp) <- ~lon+lat
projection(tresp) <- projection(shape)
gadm <- tresp %over% shape
tresp <- as.data.frame(tresp)
tresp <- cbind(gadm[ ,c(5,7,9)], tresp)
colnames(tresp) <- c("region","zone","woreda","tid","sid","lon","lat","cyld","year","loc","trt","tyld")

# project survey coords to grid CRS
tresp.proj <- as.data.frame(project(cbind(tresp$lon, tresp$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(tresp.proj) <- c("x","y")
tresp <- cbind(tresp, tresp.proj)
coordinates(tresp) <- ~x+y
projection(tresp) <- projection(tresp)

# extract gridded variables at survey locations
trespgrid <- extract(grids, tresp)
gsdat <- as.data.frame(cbind(tresp, trespgrid)) 

# Classify by site indices ------------------------------------------------
si.lmer <- lmer(log(tyld)~log(cyld)*trt+(1|sid)+(1|year), gsdat) ## random intercept (site-level) model
summary(si.lmer)
si.ran <- ranef(si.lmer) ## extract random effects
si <- as.data.frame(rownames(si.ran$sid))
si$si <- si.ran$sid[,1]
colnames(si) <- c("sid","si")
si$sic <- ifelse(si$si > 0, "A", "B") ## classify above/below average site indices (sic = A/B)
gsdat <- merge(gsdat, si, by="sid")
si <- merge(si, sites, by="sid")

# Plots
boxplot(tyld~trt, notch=T, ylab="Yield (kg/ha)", ylim=c(0,14000), gsdat) ## treatment differences
boxplot(tyld~sic, notch=T, ylab="Yield (kg/ha)", ylim=c(0,14000), gsdat) ## yield differences between site index classes
boxplot(tyld~trt*sic, notch=T, ylab="Yield (kg/ha)", ylim=c(0,14000), gsdat) ## treatment differences
par(pty="s")
plot(tyld~cyld, xlab="Control yield (kg/ha)", ylab="Treatment yield (kg/ha)", cex.lab=1.3, gsdat)

# extract gridded variables at trial locations
si.proj <- as.data.frame(project(cbind(si$lon, si$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(si.proj) <- c("x","y")
si <- cbind(si, si.proj)
coordinates(si) <- ~x+y
projection(si) <- projection(tresp)
sigrid <- extract(grids, si)
sidat <- as.data.frame(cbind(si, sigrid)) 

# Write data frames -------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/TAMASA_gsdat.csv", row.names = F)
write.csv(sidat, "./Results/TAMASA_sidat.csv", row.names = F)

# Yield trial map widget --------------------------------------------------
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(si$lon, si$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'TAMASA_trials.html', selfcontained = T) ## save widget

