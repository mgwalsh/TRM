# Ensemble  predictions of site-level control yields (Yc) and
# treatment response ratio indices (SRI) in Malawi
# Malawi LREP response trial data (courtesy of LREP & Todd Benson)
# LREP data documentation at: https://www.dropbox.com/s/4qbxnz4mdl92pdv/Malawi%20area-specific%20fertilizer%20recs%20report.pdf?dl=0
# Data pre-processing with: https://github.com/mgwalsh/TRM/blob/master/MW_LREP_SI.R
# J.Chen & M.Walsh, December 2014

# Required packages
# install.packages(c("downloader","raster","rgdal","caret","MASS","randomForest","gbm","nnet","elasticnet")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(caret)
require(MASS)
require(randomForest)
require(gbm)
require(nnet)
require(elasticnet)

# Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("LREP_Data", showWarnings=F)
dat_dir <- "./LREP_Data"

# Site index data download to "./LREP_Data"
download("https://www.dropbox.com/s/o9588q2wci8mtiv/MW_Site_Indices.csv?dl=0", "./LREP_Data/MW_Site_Indices.csv", mode="wb")
mwsite <- read.table(paste(dat_dir, "/MW_Site_Indices.csv", sep=""), header=T, sep=",")

# Malawi grids download to "./LREP_Data" (~7.6 Mb)
download("https://www.dropbox.com/s/54di5f37yp30bz4/MW_grids.zip?dl=0", "./LREP_Data/MW_grids.zip", mode="wb")
unzip("./LREP_Data/MW_grids.zip", exdir="./LREP_Data", overwrite=T)
glist <- list.files(path="./LREP_Data", pattern="tif", full.names=T)
mwgrid <- stack(glist)

# Data setup --------------------------------------------------------------
# Extract gridded variables at LREP locations
coordinates(mwsite) <- ~Easting+Northing
projection(mwsite) <- projection(mwgrid)
sitegrid <- extract(mwgrid, mwsite)

# Assemble dataframes
# Average control yield estimates (Yc, kg/ha)
Yc <- mwsite$Yc
ycdat <- cbind.data.frame(Yc, sitegrid)
ycdat <- na.omit(ycdat)

# Average site response index estimates (SRI, dimensionless)
SRI <- mwsite$SRI
sidat <- cbind.data.frame(SRI, sitegrid)
sidat <- na.omit(sidat)

# set train/test set randomization seed
seed <- 1385321
set.seed(seed)

# Split data into train and test sets ------------------------------------
# Control yield train/test split
ycIndex <- createDataPartition(ycdat$Yc, p = 0.7, list = FALSE, times = 1)
ycTrain <- ycdat[ ycIndex,]
ycTest  <- ycdat[-ycIndex,]

# Site response index train/test split
siIndex <- createDataPartition(sidat$SRI, p = 0.7, list = FALSE, times = 1)
siTrain <- sidat[ siIndex,]
siTest  <- sidat[-siIndex,]

# Stepwise main effects GLM's <MASS> --------------------------------------
# 5-fold CV
step <- trainControl(method = "cv", number = 5)

# Average control yields (Yc, kg/ha)
Yc.glm <- train(log(Yc) ~ ., data = ycTrain,
                method = "glmStepAIC",
                trControl = step)
summary(Yc.glm)
ycglm.pred <- predict(mwgrid, Yc.glm) ## spatial predictions

# Average site response indices (SRI, dimensionless)
SRI.glm <- train(SRI ~ ., data = siTrain,
                 method = "glmStepAIC",
                 trControl = step)
summary(SRI.glm)
siglm.pred <- predict(mwgrid, SRI.glm) ## spatial predictions
