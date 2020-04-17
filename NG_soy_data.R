# Soybean response trials and site indices, Central Nigeria
# Yield response data courtesy of IITA
# M. Walsh & J. van Heerwaarden, April 2020

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
dir.create("NG_soy_data", showWarnings=F)
setwd("./NG_soy_data")

# Data downloads -----------------------------------------------------------
# download Nigeria soybean yield trial data
download("https://osf.io/w6zhr?raw=1", "NG_soy_rct.zip", mode = "wb")
unzip("NG_soy_rct.zip", overwrite = T)
rct <- read.table("NG_soy_rct.csv", header=T, sep=",")

# download GADM-L2 shapefile (courtesy: http://www.gadm.org)
download("https://www.dropbox.com/s/y3h6l7yu00orm78/NGA_adm2.zip?raw=1", "NGA_adm2.zip", mode = "wb")
unzip("NGA_adm2.zip", overwrite = T)
shape <- shapefile("NGA_adm2.shp")

# download raster stack (note this is a big 900+ Mb download)
download("https://osf.io/6g8aq?raw=1", "NG_250m_2017.zip", mode = "wb")
unzip("NG_250m_2017.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# set ROI grid extent
ext <- data.frame(lat = c(5.5,5.5,12.6,12.6), lon = c(4.5,9.8,4.5,9.8)) ## set ROI extent in degrees
names(ext) <- c("lat","lon")
coordinates(ext) <- ~ lon + lat
proj4string(ext) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
ext <- spTransform(ext, CRS("+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
bb <- extent(ext)
grids <- crop(grids, bb)

# Data setup --------------------------------------------------------------
# attach Nigeria GADM-L2 admin unit names from shape
coordinates(rct) <- ~lon+lat
projection(rct) <- projection(shape)
gadm <- rct %over% shape
rct <- as.data.frame(rct)
rct <- cbind(gadm[ ,c(5,7)], rct)
colnames(rct) <- c("state","lga","year","source","study","id","lon","lat","yc","yt")

# project survey coords to grid CRS
rct.proj <- as.data.frame(project(cbind(rct$lon, rct$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(rct.proj) <- c("x","y")
rct <- cbind(rct, rct.proj)
coordinates(rct) <- ~x+y
projection(rct) <- projection(rct)

# extract gridded variables at survey locations
rctgrid <- extract(grids, rct)
gsdat <- as.data.frame(cbind(rct, rctgrid)) 

# Define unique grid ID's (GID)
# Specify pixel scale (res.pixel, in m)
res.pixel <- 250

# Grid ID (GID) definition
xgid <- ceiling(abs(gsdat$x)/res.pixel)
ygid <- ceiling(abs(gsdat$y)/res.pixel)
gidx <- ifelse(gsdat$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(gsdat$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
gsdat <- cbind(GID, gsdat)

# Classify yield response propensities by conditional mean ----------------
yt.lme <- lmer(log(yt)~log(yc)+(1|year)+(1|GID), data = gsdat)
summary(yt.lme) 
gsdat$yte <- exp(fitted(yt.lme))
par(pty="s")
plot(yt~exp(fitted(yt.lme)), xlab="Expected yield (kg/ha)", ylab="Measured yield (kg/ha)", xlim = c(-5, 4505), ylim = c(-5, 4505), cex.lab=1.1, gsdat)
abline(c(0,1), col="red", lwd=2)
gsdat$ysi <- as.factor(ifelse(exp(fitted(yt.lme, gsdat)) > gsdat$yt, "B", "A"))
table(gsdat$ysi)
boxplot(yt~ysi, notch=T, gsdat)

