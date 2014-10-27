# Case studies of Maize yield responses to fertilizer applications
# Kenya Vihiga/Siaya district NPK fertilizer response trial data courtesy of N2AFRICA
# M. Walsh, October 2014

# Set local working directory e.g.
dat_dir <- "/Users/markuswalsh/Documents/Projects/N2AFRICA"
setwd(dat_dir)

# Required packages
# install.packages(c("downloader","proj4","rgeos","dismo")), dependencies=TRUE)
require(downloader)
require(proj4)
require(rgeos)
require(dismo)

# Data downloads ----------------------------------------------------------

download("https://www.dropbox.com/s/u8o9935k3r4hdyx/VISI_trials.csv?dl=0", "VISI_trials.csv", mode="wb")
geos <- read.table("VISI_trials.csv", header=T, sep=",")

# Kenya Gtif download (~22 Mb)
download("https://www.dropbox.com/s/wg7k2ff1snge8h9/KE_grids.zip?dl=0", "KE_grids.zip", mode="wb")
unzip("KE_grids.zip", overwrite=T)

# Generate LAEA CRS & grid ID's -------------------------------------------

# Project to Africa LAEA from LonLat
geos.laea <- as.data.frame(project(cbind(geos$Lon, geos$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.laea) <- c("x","y")
geos <- cbind(geos, geos.laea)

# Generate AfSIS grid cell ID's (GID)
res.pixel <- 1000
xgid <- ceiling(abs(geos$x)/res.pixel)
ygid <- ceiling(abs(geos$y)/res.pixel)
gidx <- ifelse(geos$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(geos$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
geos.gid <- cbind(geos, GID)

# Distribution model (DM) setup -------------------------------------------

pres <- aggregate(geos.gid[,9:10], by=list(GID=geos.gid$GID), mean)
glist <- list.files(pattern='tif', full.names=T)
grids <- stack(glist)

# Generate a "x" (specify) km Region of Interest (ROI) buffer around existing GID's
x <- 25000
coordinates(pres) <- ~x+y
proj4string(pres) <- CRS("+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=5 +lon_0=20 +units=m +no_defs")
buffer <- circles(pres, d=x, lonlat=F)
roi <- gUnaryUnion(buffer@polygons)

# Randomly sample the ROI background with B (specify) * no. of trial GID's samples
B <- 10
ext <- extent(roi)
set.seed(1235813)
back <- randomPoints(grids, n=B*length(pres), ext=ext, extf = 1)
gback <- extract(grids, back)

# Plot ROI, trial GID's and background sample locations
plot(ext, xlab="Easting (m)", ylab="Northing (m)")
plot(roi, add=T)
points(back, pch=3, col="grey", cex=0.5)
points(pres, pch=21, col="red", bg="red")

# Dataframe
gpres <- extract(grids, pres)
pb <- c(rep(1, nrow(pres)), rep(0, nrow(back)))
presback <- data.frame(cbind(pb, rbind(gpres, gback)))

# Dismo profile models (just for illustration) ----------------------------

# Bioclim (for comparison to Mahalanobis distance only)
bc <- bioclim(grids, pres)
ebc <- evaluate(pres, back, bc, grids)
plot(ebc, "ROC")

# Mahalanobis distances from trial locations
mh <- mahal(grids, pres)
emh <- evaluate(pres, back, mh, grids)
plot(emh, "ROC")

# Plots
mhp <- predict(grids, mh, ext=ext)
mhp[mhp < -1000] <- -1000
plot(mhp)
mht <- threshold(emh, "spec_sens")
plot(mhp>mht)

# Regression models -------------------------------------------------------

# stepwise main effects GLM
require(MASS)
pres.glm <- glm(pb ~ ., family = binomial(link="logit"), data=presback)
step <- stepAIC(pres.glm)
summary(step)
pglm <- predict(grids, step, type="response", ext=ext)
plot(ext, xlab="Easting (m)", ylab="Northing (m)")
plot(pglm, add=T)
plot(roi, add=T)
points(back, pch=3, col="black", cex=0.5)
points(pres, pch=21, col="red", bg="red")

# Random forest (no tuning default)
require(randomForest)
set.seed(1235)
pres.rf <- randomForest(pb ~ ., importance=T, proximity=T, data=presback)
print(pres.rf)
prf <- predict(grids, pres.rf, ext=ext)
plot(ext, xlab="Easting (m)", ylab="Northing (m)")
plot(prf, add=T)
plot(roi, add=T)
points(back, pch=3, col="black", cex=0.5)
points(pres, pch=21, col="red", bg="red")

# Export Gtifs
writeRaster(pglm, filename="VISI_glm", format="Gtif", overwrite=T)
writeRaster(prf, filename="VISI_glm", format="Gtif", overwrite=T)
