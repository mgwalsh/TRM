#+ Example site index models from AfSIS omission trial data
#+ AfSIS ommission trial sample data, courtesy of CIAT
#+ M. Walsh, May 2015

#+ Required packages
# install.packages(c("downloader","caret","fastICA","glmnet")), dependencies=TRUE)
require(downloader)
require(caret)
require(fastICA)
require(glmnet)

#+ Data download ---------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("OT_data", showWarnings=F)
dat_dir <- "./OT_data"

# AfSIS omission trial data download to ... "./OT_data"
download("https://www.dropbox.com/s/rs68t4077iafo17/AfSIS_OT.zip?dl=0", "./OT_data/AfSIS_OT.zip", mode="wb")
unzip("./OT_data/AfSIS_OT.zip", exdir="./OT_data", overwrite=T)
field <- read.table(paste(dat_dir, "/Fields.csv", sep=""), header=T, sep=",")
spect <- read.table(paste(dat_dir, "/Spectra.csv", sep=""), header=T, sep=",")

#+ Spectral data transforms <caret> --------------------------------------
# Principal components analysis (PCA) transform
pcatr <- preProcess(spect[,3:1762], method=c("BoxCox", "center", "scale", "pca"), pcaComp=3)
sppca <- predict(pcatr, spect[,3:1762])

# not run: Independent components analysis (ICA) transform
# icatr <- preProcess(spect[,3:1762], method=c("BoxCox", "center", "scale", "ica"), n.comp=3)
# spica <- predict(icatr, spect[,3:1762])

#+ Merge spectral transforms with field measurements ---------------------
speci <- cbind(spect[,1:2], sppca)
avesp <- aggregate(speci[,3:5], by=list(speci$FieldID), FUN="mean")
colnames(avesp) <- c("FieldID","PCA1","PCA2","PCA3")
afotd <- merge(field, avesp, by="FieldID")

#+ Write data file --------------------------------------------------------
write.csv(afotd, "./OT_data/OT_data.csv", row.names=F)

#+ Regularized regression models ------------------------------------------
set.seed(1235813)

tc <- trainControl(
  method = "repeatedCV",
  number = 10,
  repeats = 5,
  returnResamp = "all"
 )

# Control yields (Yc)
# Yc with spectral covariates
Yc.spca <- train(log(Yc) ~ PCA1 + PCA2 + PCA3, data = afotd,
                family = "gaussian", 
                method = "glmnet",
                tuneGrid = expand.grid(.alpha=seq(0.1,0.5, by=0.1),.lambda=seq(0,0.3,by=0.01)),
                trControl = tc)

# Yc with wet chemistry covariates
Yc.wetc <- train(log(Yc) ~ pH + Sand + N + P + K, data = afotd,
                standardize = TRUE,
                family = "gaussian", 
                method = "glmnet",
                tuneGrid = expand.grid(.alpha=seq(0.1,0.5, by=0.1),.lambda=seq(0,0.3,by=0.01)),
                trControl = tc)

# Treatment response ratios (TRR)
# TRR with spectral covariates
TRR.spca <- train(TRR ~ PCA1 + PCA2 + PCA3, data = afotd,
                  family = "gaussian", 
                  method = "glmnet",
                  tuneGrid = expand.grid(.alpha=seq(0.5,1, by=0.1),.lambda=seq(0,0.3,by=0.01)),
                  trControl = tc)

# TRR with wet chemistry covariates
TRR.wetc <- train(TRR ~ pH + Sand + N + P + K, data = afotd,
                  standardize = TRUE,
                  family = "gaussian", 
                  method = "glmnet",
                  tuneGrid = expand.grid(.alpha=seq(0.1,1, by=0.1),.lambda=seq(0,0.3,by=0.01)),
                  trControl = tc)

