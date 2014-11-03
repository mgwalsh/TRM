# Spatial models of site-level control yields (Yc) and treatment response ratio indices (SRI)
# Malawi LREP response trial data (courtesy of LREP & Todd Benson)
# LREP data documentation at: https://www.dropbox.com/s/4qbxnz4mdl92pdv/Malawi%20area-specific%20fertilizer%20recs%20report.pdf?dl=0
# Data pre-processing with: https://github.com/mgwalsh/TRM/blob/master/MW_LREP_SI.R
# M. Walsh, October 2014

# Set local working directory e.g.
dat_dir <- "~/Documents/Data/Malawi/Fert_resp_models"
setwd(dat_dir)

# Required packages
# install.packages(c("downloader","proj4","raster")), dependencies=TRUE)
require(downloader)
require(proj4)
require(raster)

# Data download -----------------------------------------------------------

download("https://www.dropbox.com/s/o9588q2wci8mtiv/MW_Site_Indices.csv?dl=0", "MW_Site_Indices.csv", mode="wb")
mwsite <- read.table("MW_Site_Indices.csv", header=T, sep=",")

# Malawi grids download (~7.2 Mb)
download("https://www.dropbox.com/s/54di5f37yp30bz4/MW_grids.zip?dl=0", "MW_grids.zip", mode="wb")
unzip("MW_grids.zip", overwrite=T)
glist <- list.files(pattern="tif", full.names=T)
mwgrid <- stack(glist)

# Overlay gridded covariates & generate dataframes ------------------------

coordinates(mwsite) <- ~Easting+Northing
projection(mwsite) <- projection(mwgrid)
exgrid <- extract(mwgrid, mwsite)
Yc <- mwsite$Yc
SRI <- mwsite$SRI
ycdat <- data.frame(cbind(Yc, exgrid))
ycdat <- na.omit(ycdat)
srdat <- data.frame(cbind(SRI, exgrid))
srdat <- na.omit(srdat)

# Regression models -------------------------------------------------------

# Stepwise main effects GLM's
require(MASS)
## Control yield predictions (Yc)
Yc.glm <- glm(Yc ~ ., family=gaussian(link="log"), data=ycdat)
Yc.step <- stepAIC(Yc.glm)
ycglm <- predict(mwgrid, Yc.step, type="response")
plot(ycglm)
## Site response index predictions (SRI)
SRI.glm <- glm(SRI ~ ., family=gaussian, data=srdat)
SRI.step <- stepAIC(SRI.glm)
sriglm <- predict(mwgrid, SRI.step, type="response")
plot(sriglm)

# Regression trees
require(rpart)
## Control yield predictions (Yc)
Yc.rt <- rpart(Yc ~ ., data=ycdat)
ycrt <- predict(mwgrid, Yc.rt)
plot(ycrt)
## Site response index predictions (SRI)
SRI.rt <- rpart(SRI ~ ., data=srdat)
srirt <- predict(mwgrid, SRI.rt)
plot(srirt)

# Random forests (no tuning default)
require(randomForest)
## Control yield predictions (Yc)
Yc.rf <- randomForest(Yc ~ ., importance=T, proximity=T, data=ycdat)
ycrf <- predict(mwgrid, Yc.rf)
plot(ycrf)
## Site response index predictions (SRI)
SRI.rf <- randomForest(SRI ~ ., importance=T, proximity=T, data=srdat)
srirf <- predict(mwgrid, SRI.rf)
plot(srirf)

# Unweighted mean prediction ensemble (glm, rt & rf models)
## Control yield predictions (Yc)
myc <- mean(ycglm, ycrt, ycrf)
plot(myc)
## Site response index predictions (SRI)
msri <- mean(sriglm, srirt, srirf)
plot(msri)

# Write regression predictions --------------------------------------------

ycpred <- stack(ycglm, ycrt, ycrf)
names(ycpred) <- c("ycglm", "ycrt", "ycrf")
if (require(rgdal)){
  yc <- writeRaster(ycpred, filename="ycpred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
}

sripred <- stack(sriglm, srirt, srirf)
names(ycpred) <- c("sriglm", "srirt", "srirf")
if (require(rgdal)){
  yc <- writeRaster(sripred, filename="sripred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
}

