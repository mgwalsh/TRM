# Case studies of Maize yield responses to fertilizer applications
# Kenya Vihiga/Siaya district NPK fertilizer response trial data courtesy of N2AFRICA
# M. Walsh, October 2014

# Set local working directory e.g.
dat_dir <- "/Users/markuswalsh/Documents/Projects/N2AFRICA"
setwd(dat_dir)

# Required packages
# install.packages(c("downloader","proj4","dismo")), dependencies=TRUE)
require(downloader)
require(proj4)
require(dismo)

# Data downloads ----------------------------------------------------------

download("https://www.dropbox.com/s/u8o9935k3r4hdyx/VISI_trials.csv?dl=0", "VISI_trials.csv", mode="wb")
geot <- read.table("VISI_trials.csv", header=T, sep=",")

# Kenya Gtif download (~22 Mb)
download("https://www.dropbox.com/s/wg7k2ff1snge8h9/KE_grids.zip?dl=0", "KE_grids.zip", mode="wb")
unzip("KE_grids.zip", overwrite=T)

# generate LAEA CRS & grid ID's -------------------------------------------

# Project to Africa LAEA from LonLat
geot.laea <- as.data.frame(project(cbind(geot$Lon, geot$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geot.laea) <- c("x","y")
geot <- cbind(geot, geot.laea)

# Generate AfSIS grid cell ID's (GID)
res.pixel <- 1000
xgid <- ceiling(abs(geot$x)/res.pixel)
ygid <- ceiling(abs(geot$y)/res.pixel)
gidx <- ifelse(geot$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(geot$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
geot.gid <- cbind(geot, GID)
coordinates(geot.gid) = ~x+y
proj4string(geot.gid) = CRS("+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=5 +lon_0=20 +no_defs")

# Grid overlay ------------------------------------------------------------

grid.list <- c("BSANs.tif","BSASs.tif","BSAVs.tif","CTIs.tif","ELEVs.tif","EVIAs.tif","LSTDs.tif","LSTNs.tif","REF1s.tif","REF2s.tif","REF3s.tif","REF7s.tif","RELIs.tif","TMAPs.tif","TMFIs.tif")
for (i in 1:length(grid.list)){
  print(paste("extracting", grid.list[i]))
  grids <- raster(grid.list[i])
  geot.gid@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract(
    x = grids, 
    y = geot.gid,
    method = "simple")
}
tgrid <- as.data.frame(geot.gid)
write.csv(tgrid, "VISI_dat.csv")

# Environmental similarities to trial loacations (LID'S) -------------------

LID <- aggregate(tgrid[,1:2], by=list(LID=tgrid$LID), mean)
glist <- list.files(pattern='tif', full.names=TRUE)
grids <- stack(glist)
mahal <- mahal(grids, LID[,2:3])
preds <- predict(grids, mahal)



