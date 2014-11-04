# Ensemble regression predictions of site-level control yields (Yc) 
# and treatment response ratio indices (SRI)
# Malawi LREP response trial data (courtesy of LREP & Todd Benson)
# LREP data documentation at: https://www.dropbox.com/s/4qbxnz4mdl92pdv/Malawi%20area-specific%20fertilizer%20recs%20report.pdf?dl=0
# Data pre-processing with: https://github.com/mgwalsh/TRM/blob/master/MW_LREP_SI.R
# M. Walsh, J. Chen & A. Verlinden, November 2014

# Set local working directory e.g.
dat_dir <- "~/Documents/Data/Malawi/Fert_resp_models"
setwd(dat_dir)

# Required packages
# install.packages(c("downloader","proj4","raster")), dependencies=TRUE)
require(downloader)
require(proj4)
require(raster)

# Data downloads ----------------------------------------------------------
download("https://www.dropbox.com/s/o9588q2wci8mtiv/MW_Site_Indices.csv?dl=0", "MW_Site_Indices.csv", mode="wb")
mwsite <- read.table("MW_Site_Indices.csv", header=T, sep=",")

# Malawi grids download (~7.5 Mb)
download("https://www.dropbox.com/s/54di5f37yp30bz4/MW_grids.zip?dl=0", "MW_grids.zip", mode="wb")
unzip("MW_grids.zip", overwrite=T)
glist <- list.files(pattern="tif", full.names=T)
mwgrid <- stack(glist)

# Split the site data into train and test sets ----------------------------
# Hold-out ~1/4 to 1/3 of the site data
set.seed(1385321)
index <- 1:nrow(mwsite)
testn <- sample(index, trunc(length(index)/4))
test <- na.omit(mwsite[testn,-6])
train <- na.omit(mwsite[-testn,-6])

# Overlay grids & generate dataframes -------------------------------------
coordinates(train) <- ~Easting+Northing
projection(train) <- projection(mwgrid)
exgrid <- extract(mwgrid, train)
Yc <- train$Yc
SRI <- train$SRI
ycdat <- data.frame(cbind(Yc, exgrid))
ycdat <- na.omit(ycdat)
srdat <- data.frame(cbind(SRI, exgrid))
srdat <- na.omit(srdat)

# Regression models -------------------------------------------------------
# Stepwise main effects GLM's
require(MASS)
# Control yield predictions (Yc)
Yc.glm <- glm(Yc ~ ., family=gaussian(link="log"), ycdat)
Yc.step <- stepAIC(Yc.glm)
ycglm <- predict(mwgrid, Yc.step, type="response")

# Site response index predictions (SRI)
SRI.glm <- glm(SRI ~ ., family=gaussian, data=srdat)
SRI.step <- stepAIC(SRI.glm)
sriglm <- predict(mwgrid, SRI.step, type="response")

# Regression trees
require(rpart)
# Control yield predictions (Yc)
Yc.rt <- rpart(Yc ~ ., data=ycdat)
ycrt <- predict(mwgrid, Yc.rt)

# Site response index predictions (SRI)
SRI.rt <- rpart(SRI ~ ., data=srdat)
srirt <- predict(mwgrid, SRI.rt)

# Random forests (no tuning default)
require(randomForest)
# Control yield predictions (Yc)
Yc.rf <- randomForest(Yc ~ ., importance=T, proximity=T, data=ycdat)
ycrf <- predict(mwgrid, Yc.rf)

# Site response index predictions (SRI)
SRI.rf <- randomForest(SRI ~ ., importance=T, proximity=T, data=srdat)
srirf <- predict(mwgrid, SRI.rf)

# Regression ensembles ---------------------------------------------------
# Dataframe setup
ycpred <- stack(ycglm, ycrt, ycrf)
names(ycpred) <- c("ycglm", "ycrt", "ycrf")
sripred <- stack(sriglm, srirt, srirf)
names(sripred) <- c("sriglm", "srirt", "srirf")
exyc <- extract(ycpred, train)
exyc <- data.frame(cbind(Yc, exyc))
exyc <- na.omit(exyc)
exsri <- extract(sripred, train)
exsri <- data.frame(cbind(SRI, exsri))
exsri <- na.omit(exsri)

# Control yield predictions (Yc) 
YCwgt.glm <- glm(Yc~log(ycglm)+log(ycrt)+log(ycrf), family=gaussian(link="log"), data=exyc)
summary(YCwgt.glm)
plot(Yc~fitted(YCwgt.glm), exyc)
ycwgt <- predict(ycpred, YCwgt.glm, type="response")
quantile(ycwgt, prob=c(0.025,0.25,0.5,0.75,0.975))
plot(ycwgt)

# Site response index predictions (SRI)
SRIwgt.glm <- glm(SRI ~ ., family=gaussian, data=exsri)
summary(SRIwgt.glm)
plot(SRI~fitted(SRIwgt.glm), exsri)
sriwgt <- predict(sripred, SRIwgt.glm, type="response")
quantile(sriwgt, prob=c(0.025,0.25,0.5,0.75,0.975))
plot(sriwgt)

# Test set validation -----------------------------------------------------
coordinates(test) <- ~Easting+Northing
projection(test) <- projection(mwgrid)
ycpred <- extract(ycwgt, test)
sripred <- extract(sriwgt, test)
pretest <- cbind(as.data.frame(test), ycpred, sripred)
plot(Yc~ycpred, pretest)
plot(SRI~sripred, pretest)

# Write predictions -------------------------------------------------------
dir.create("Results", recursive=F)
writeRaster(ycpred, filename="./Results/MW_ycpred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(sripred, filename="./Results/MW_sripred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
enspred <- stack(ycwgt, sriwgt)
names(enspred) <- c("ycwgt", "sriwgt")
writeRaster(enspred, filename="./Results/MW_enspred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)