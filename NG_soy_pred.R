# Stacked spatial predictions of Nigeria soybean trial site indices
# M. Walsh, April 2020

# Required packages
# install.packages(c("devtools","caret","MASS","randomForest","gbm","nnet","glmnet","plyr","doParallel","dismo")), dependencies=T)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(MASS)
  require(randomForest)
  require(gbm)
  require(nnet)
  require(glmnet)
  require(plyr)
  require(doParallel)
  require(dismo)
})

# Data setup --------------------------------------------------------------
# SourceURL <- "https://github.com/mgwalsh/TRM/blob/master/OAF_yield_data%202016:17.R"
# source_url(SourceURL)
rm(list=setdiff(ls(), c("sidat","grids"))) ## scrub extraneous objects in memory

# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(sidat$sic, p = 4/5, list = F, times = 1)
gs_cal <- sidat[ gsIndex,]
gs_val <- sidat[-gsIndex,]

# Survey calibration labels
gl_cal <- gs_cal$sic

# raster calibration features
gf_cal <- gs_cal[,11:53]

# GLM <MASS> --------------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
gl <- train(gf_cal, gl_cal, 
            method = "glmStepAIC",
            family = "binomial",
            preProc = c("center","scale"), 
            trControl = tc,
            metric ="ROC")

# model outputs & predictions
summary(gl)
print(gl) ## ROC's accross cross-validation
gl.pred <- predict(grids, gl, type = "prob") ## spatial predictions
stopCluster(mc)

# Random forest <randomForest> --------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)
tg <- expand.grid(mtry = seq(1,5, by=1)) ## model tuning steps

# model training
rf <- train(gf_cal, gl_cal,
            preProc = c("center","scale"),
            method = "rf",
            ntree = 501,
            metric = "ROC",
            tuneGrid = tg,
            trControl = tc)

# model outputs & predictions
print(rf) ## ROC's accross tuning parameters
plot(varImp(rf))
rf.pred <- predict(grids, rf, type = "prob") ## spatial predictions
stopCluster(mc)

# Generalized boosting <gbm> ----------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, summaryFunction = twoClassSummary,
                   allowParallel = T)

## for initial <gbm> tuning guidelines see @ https://stats.stackexchange.com/questions/25748/what-are-some-useful-guidelines-for-gbm-parameters
tg <- expand.grid(interaction.depth = seq(6,14, by=2), shrinkage = 0.01, n.trees = 501,
                  n.minobsinnode = 25) ## model tuning steps

# model training
gb <- train(gf_cal, gl_cal, 
            method = "gbm", 
            preProc = c("center", "scale"),
            trControl = tc,
            tuneGrid = tg,
            metric = "ROC")

# model outputs & predictions
print(gb) ## ROC's accross tuning parameters
plot(varImp(gb))
gb.pred <- predict(grids, gb, type = "prob") ## spatial predictions
stopCluster(mc)

# Neural network <nnet> ---------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)
tg <- expand.grid(size = seq(6,14, by=2), decay = 0.01) ## model tuning steps

# model training
nn <- train(gf_cal, gl_cal, 
            method = "nnet",
            preProc = c("center","scale"), 
            tuneGrid = tg,
            trControl = tc,
            metric ="ROC")

# model outputs & predictions
print(nn) ## ROC's accross tuning parameters
plot(varImp(nn))
nn.pred <- predict(grids, nn, type = "prob") ## spatial predictions
stopCluster(mc)

# Model stacking setup ----------------------------------------------------
preds <- stack(gl.pred, rf.pred, gb.pred, nn.pred)
names(preds) <- c("gl","rf", "gb","nn")
plot(preds, axes = F)

# extract model predictions
coordinates(gs_val) <- ~x+y
projection(gs_val) <- projection(preds)
gspred <- extract(preds, gs_val)
gspred <- as.data.frame(cbind(gs_val, gspred))

# stacking model validation labels and features
gl_val <- gspred$sic
gf_val <- gspred[,54:57] ## subset validation features

# Model stacking ----------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, 
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
st <- train(gf_val, gl_val,
            method = "glm",
            family = "binomial",
            metric = "ROC",
            trControl = tc)

# model outputs & predictions
print(st)
summary(st)
plot(varImp(st))
st.pred <- predict(preds, st, type = "prob") ## spatial predictions of soy site indices
plot(st.pred, axes = F)
stopCluster(mc)

# Write prediction grids --------------------------------------------------
sipreds <- stack(gl.pred, rf.pred, gb.pred, nn.pred, st.pred)
names(sipreds) <- c("gl","rf","gb","nn","st")
writeRaster(sipreds, filename="./results/NG_soy_si_preds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Receiver-Operator characteristics ---------------------------------------
preds <- stack(gl.pred, rf.pred, gb.pred, nn.pred, st.pred)
names(preds) <- c("gl","rf", "gb","nn","st")

# extract model predictions
coordinates(sidat) <- ~x+y
projection(sidat) <- projection(preds)
sipre <- extract(preds, sidat)
sidat <- as.data.frame(cbind(sidat, sipre))

# ROC
siA <- subset(sidat, sidat$sic=='A', select=st)
siB <- subset(sidat, sidat$sic=='B', select=st)
si_eval <- evaluate(p=siA[,1], a=siB[,1]) ## calculate ROC's on full dataset
plot(si_eval, 'ROC') ## plot the ROC curve
t <- threshold(si_eval) ## calculate thresholds based on ROC
sidat$sit <- ifelse(sidat$st > t[,1], "A", "B") ## predicted classification threshold using kappa
confusionMatrix(data = sidat$sit, reference = sidat$sic, positive = "A")

# ECDF plot of soybean yield responses by modelled site index categories
soyA <- subset(sidat, sit=='A', select=c(yc,yt)) 
soyB <- subset(sidat, sit=='B', select=c(yc,yt)) 
par(pty="s")
plot(ecdf(soyA$yt-soyA$yc), verticals=T, lty=1, lwd=1, col="dark green", do.points=F, main="",
     xlab="Soybean yield response (kg/ha)", ylab="Cum. proportion of observations", xlim=c(-500,2500), cex.lab=1.2)
plot(ecdf(soyB$yt-soyB$yc), add=T, verticals=T, lty=1, lwd=1, col="red", do.points=F)
abline(0.5,0, lty=1, col="grey")

# Write output data frame ---------------------------------------------------
write.csv(sidat, "./results/NG_soy_si_pred.csv", row.names = F)

