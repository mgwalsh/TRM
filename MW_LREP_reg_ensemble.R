# Ensemble regression predictions of site-level control yields (Yc) 
# and treatment response ratio indices (SRI) in Malawi.
# Malawi LREP response trial data (courtesy of LREP & Todd Benson)
# LREP data documentation at: https://www.dropbox.com/s/4qbxnz4mdl92pdv/Malawi%20area-specific%20fertilizer%20recs%20report.pdf?dl=0
# Data pre-processing with: https://github.com/mgwalsh/TRM/blob/master/MW_LREP_SI.R
# M.Walsh, J.Chen & A.Verlinden, November 2014

# Required packages
# install.packages(c("downloader","raster","MASS","rpart","randomForest")), dependencies=TRUE)
require(downloader)
require(raster)
require(MASS)
require(rpart)
require(randomForest)

# Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("Data", showWarnings=F)
dat_dir <- "./Data"

# Site index data download to "./Data"
download("https://www.dropbox.com/s/o9588q2wci8mtiv/MW_Site_Indices.csv?dl=0", "./Data/MW_Site_Indices.csv", mode="wb")
mwsite <- read.table(paste(dat_dir, "/MW_Site_Indices.csv", sep=""), header=T, sep=",")

# Malawi grids download to "./Data" (~7.5 Mb)
download("https://www.dropbox.com/s/54di5f37yp30bz4/MW_grids.zip?dl=0", "./Data/MW_grids.zip", mode="wb")
unzip("./Data/MW_grids.zip", exdir="./Data", overwrite=T)
glist <- list.files(path="./Data", pattern="tif", full.names=T)
mwgrid <- stack(glist)

# Randomly split the site index data into train and test sets -------------
# Hold-out ~1/4 to 1/3 of the sites, set this with n
n <- 4
set.seed(1385321)
index <- 1:nrow(mwsite)
testn <- sample(index, trunc(length(index)/n))
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
# Control yield predictions (Yc)
Yc.glm <- glm(Yc ~ ., family=gaussian(link="log"), ycdat)
Yc.step <- stepAIC(Yc.glm)
ycglm <- predict(mwgrid, Yc.step, type="response")

# Site response index predictions (SRI)
SRI.glm <- glm(SRI ~ ., family=gaussian, data=srdat)
SRI.step <- stepAIC(SRI.glm)
sriglm <- predict(mwgrid, SRI.step, type="response")

# Regression trees
# Control yield predictions (Yc)
Yc.rt <- rpart(Yc ~ ., data=ycdat)
ycrt <- predict(mwgrid, Yc.rt)

# Site response index predictions (SRI)
SRI.rt <- rpart(SRI ~ ., data=srdat)
srirt <- predict(mwgrid, SRI.rt)

# Random forests (no tuning default)
# Control yield predictions (Yc)
Yc.rf <- randomForest(Yc ~ ., importance=T, proximity=T, data=ycdat)
ycrf <- predict(mwgrid, Yc.rf)

# Site response index predictions (SRI)
SRI.rf <- randomForest(SRI ~ ., importance=T, proximity=T, data=srdat)
srirf <- predict(mwgrid, SRI.rf)

# Test set ensemble predictions ------------------------------------------
# Dataframe setup
coordinates(test) <- ~Easting+Northing
projection(test) <- projection(mwgrid)
Yc_test <- test$Yc
SRI_test <- test$SRI
ycpred <- stack(ycglm, ycrt, ycrf)
names(ycpred) <- c("ycglm", "ycrt", "ycrf")
sripred <- stack(sriglm, srirt, srirf)
names(sripred) <- c("sriglm", "srirt", "srirf")

# Overlay training set predictions w. test data
exyc <- extract(ycpred, test)
exyc <- data.frame(cbind(Yc_test, exyc))
exyc <- na.omit(exyc)
exsri <- extract(sripred, test)
exsri <- data.frame(cbind(SRI_test, exsri))
exsri <- na.omit(exsri)

# Control yield predictions (Yc) 
YCwgt.glm <- glm(Yc_test~log(ycglm)+log(ycrt)+log(ycrf), family=gaussian(link="log"), data=exyc)
YCwgt.step <- stepAIC(YCwgt.glm)
summary(YCwgt.step)
plot(Yc_test~fitted(YCwgt.step), exyc)
ycwgt <- predict(ycpred, YCwgt.step, type="response")
quantile(ycwgt, prob=c(0.025,0.25,0.5,0.75,0.975))
plot(ycwgt)

# Site response index predictions (SRI)
SRIwgt.glm <- glm(SRI_test ~ ., family=gaussian, data=exsri)
SRIwgt.step <- stepAIC(SRIwgt.glm)
summary(SRIwgt.step)
plot(SRI_test~fitted(SRIwgt.step), exsri)
sriwgt <- predict(sripred, SRIwgt.step, type="response")
quantile(sriwgt, prob=c(0.025,0.25,0.5,0.75,0.975))
plot(sriwgt)

# Write spatial predictions -----------------------------------------------
# Create a "Results" folder in your current working directory
dir.create("Results", showWarnings=F)

# Export Gtif's to "./Results"
writeRaster(ycpred, filename="./Results/MW_ycpred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(sripred, filename="./Results/MW_sripred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
enspred <- stack(ycwgt, sriwgt)
names(enspred) <- c("ycwgt", "sriwgt")
writeRaster(enspred, filename="./Results/MW_enspred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)