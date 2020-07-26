# OCP/IITA 2017 Maize fertilizer response trials, Central Nigeria
# Yield response data courtesy of IITA
# M. Walsh & J. Huising, July 2018

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
dir.create("OCP_data", showWarnings=F)
setwd("./OCP_data")

# Data downloads -----------------------------------------------------------
# download IITA/OCP yield data
download("https://www.dropbox.com/s/hi75cnp3ejr4srk/OCP_trials.zip?raw=1", "OCP_trials.zip", mode = "wb")
unzip("OCP_trials.zip", overwrite = T)
sites <- read.table("sites.csv", header=T, sep=",")
trial <- read.table("trials.csv", header=T, sep=",")
tresp <- merge(sites, trial, by="sid")
tresp <- tresp[complete.cases(tresp[ ,c(11:12)]),] ## removes incomplete cases

# download GADM-L2 shapefile (courtesy: http://www.gadm.org)
download("https://www.dropbox.com/s/y3h6l7yu00orm78/NGA_adm2.zip?raw=1", "NGA_adm2.zip", mode = "wb")
unzip("NGA_adm2.zip", overwrite = T)
shape <- shapefile("NGA_adm2.shp")

# download raster stack (note this is a big 800+ Mb download)
download("https://www.dropbox.com/s/u5fyjbujf0d7q43/NG_250m_2017.zip?raw=1", "NG_250m_2017.zip", mode = "wb")
unzip("NG_250m_2017.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

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
si$si <- si.ran$sid[,1]
colnames(si) <- c("sid","si")
si$sic <- ifelse(si$si > 0, "A", "B") ## classify above/below average site indices (sic = A or B)

# Classify by harvest indices ---------------------------------------------
gsdat$hi <- gsdat$tyld/gsdat$ayld ## grain yield / biomass
hi.lmer <- lmer(hi~trt+(1|sid), gsdat) ## random intercept (site-level) model
display(hi.lmer)
hi.ran <- ranef(hi.lmer) ## extract random effects
hi <- as.data.frame(rownames(hi.ran$sid))
si$hi <- hi.ran$sid[,1]
si$hic <- ifelse(si$hi > 0, "H", "L") ## classify above/below average harvest indices (hic = H or L)
gsdat <- merge(gsdat, si, by="sid")
si <- merge(si, sites, by="sid")

# Diagnostic plots --------------------------------------------------------
boxplot(tyld~trt, notch=T, ylab="Maize grain yield (kg / ha)", ylim=c(0,8000), gsdat) ## treatment differences
boxplot(tyld~sic, notch=T, ylab="Maize grain yield (kg / ha)", ylim=c(0,8000), gsdat) ## yield differences between site index classes
boxplot(tyld~hic, notch=T, ylab="Maize grain yield (kg / ha)", ylim=c(0,8000), gsdat) ## yield differences between site index classes
boxplot(tcob~trt*sic, notch=T, ylab="Number of cobs", ylim=c(0,800), gsdat) ## treatment differences
boxplot(tyld~trt*sic, notch=T, ylab="Maize grain yield (kg / ha)", ylim=c(0,8000), gsdat) ## treatment differences
boxplot(tyld~trt*hic, notch=T, ylab="Maize grain yield (kg / ha)", ylim=c(0,8000), gsdat) ## treatment differences

# Yield component plots
par(pty="s")
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(tyld~cyld, xlab="Maize grain yield (kg / ha), circular plot", ylab="Maize grain yield (kg / ha), total plot", gsdat)
plot(tyld~ayld, xlab="Maize biomass (kg DM / ha)", ylab="Maize grain yield (kg / ha)", gsdat)
plot(hi~si, xlab="Site index", ylab="Harvest index", si)

# Extract gridded variables at trial locations ----------------------------
si.proj <- as.data.frame(project(cbind(si$lon, si$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(si.proj) <- c("x","y")
si <- cbind(si, si.proj)
coordinates(si) <- ~x+y
projection(si) <- projection(tresp)
sigrid <- extract(grids, si)
sidat <- as.data.frame(cbind(si, sigrid)) 

# Write data frames -------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/OCP_gsdat.csv", row.names = F)
write.csv(sidat, "./Results/OCP_sidat.csv", row.names = F)

# Yield trial map widget --------------------------------------------------
w <- leaflet() %>% 
  setView(lng = mean(sidat$lon), lat = mean(sidat$lat), zoom = 7) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(sidat$lon, sidat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'OCP_trials.html', selfcontained = T) ## save widget
