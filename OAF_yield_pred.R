# Stacked spatial predictions of 2016/2017 OAF maize yield site indices & yield potentials
# M. Walsh, July 2020

# Required packages
# install.packages(c("devtools","caret","mgcv","MASS","randomForest","gbm","nnet","plyr","doParallel","dismo")), dependencies=T)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(mgcv)
  require(MASS)
  require(randomForest)
  require(gbm)
  require(nnet)
  require(plyr)
  require(doParallel)
  require(dismo)
})

# Data setup --------------------------------------------------------------
rm(list=setdiff(ls(), c("gsdat","grids","glist"))) ## scrub extraneous objects in memory

# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(gsdat$qy, p = 4/5, list = F, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# Site index calibration labels
labs <- c("qy")
lcal <- as.vector(t(gs_cal[labs]))

# raster calibration features
fcal <- gs_cal[,13:31,35:58]

# Spatial trend model <mgcv> -----------------------------------------------
# select spatial coordinates
gf_cpv <- gs_cal[,32:34]

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, 
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
gm <- train(gf_cpv, lcal, 
            method = "gam",
            preProc = c("center","scale"), 
            family = "binomial",
            metric = "ROC",
            trControl = tc)

# model outputs & predictions
summary(gm)
gm.pred <- predict(grids, gm, type = "prob") ## spatial predictions
stopCluster(mc)
# plot(gm.pred, axes=F) ## plot of the spatial prediction
fname <- paste("./Results/", labs, "_gm.rds", sep = "")
saveRDS(gm, fname)

# Central place theory model <glm> -----------------------------------------
# select central place covariates
gf_cpv <- gs_cal[,20:31]

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
gl1 <- train(gf_cpv, lcal, 
             method = "glmStepAIC",
             family = "binomial",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="ROC")

# model outputs & predictions
summary(gl1)
print(gl1) ## ROC's accross cross-validation
gl1.pred <- predict(grids, gl1, type = "prob") ## spatial predictions
stopCluster(mc)
# plot(gl1.pred, axes=F)
fname <- paste("./Results/", labs, "_gl1.rds", sep = "")
saveRDS(gl1, fname)

# GLM with all covariates -------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
gl2 <- train(fcal, lcal, 
             method = "glmStepAIC",
             family = "binomial",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="ROC")

# model outputs & predictions
summary(gl2)
print(gl2) ## ROC's accross cross-validation
gl2.pred <- predict(grids, gl2, type = "prob") ## spatial predictions
stopCluster(mc)
# plot(gl2.pred, axes=F)
fname <- paste("./Results/", labs, "_gl2.rds", sep = "")
saveRDS(gl2, fname)

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
rf <- train(fcal, lcal,
            preProc = c("center","scale"),
            method = "rf",
            ntree = 501,
            metric = "ROC",
            tuneGrid = tg,
            trControl = tc)

# model outputs & predictions
print(rf) ## ROC's accross tuning parameters
rf.pred <- predict(grids, rf, type = "prob") ## spatial predictions
stopCluster(mc)
# plot(rf.pred, axes=F)
fname <- paste("./Results/", labs, "_rf.rds", sep = "")
saveRDS(rf, fname)

# Generalized boosting <gbm> ----------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, summaryFunction = twoClassSummary,
                   allowParallel = T)

## for initial <gbm> tuning guidelines see @ https://stats.stackexchange.com/questions/25748/what-are-some-useful-guidelines-for-gbm-parameters
tg <- expand.grid(interaction.depth = seq(2,5, by=1), shrinkage = 0.01, n.trees = seq(101,501, by=50),
                  n.minobsinnode = 50) ## model tuning steps

# model training
gb <- train(fcal, lcal, 
            method = "gbm", 
            preProc = c("center", "scale"),
            trControl = tc,
            tuneGrid = tg,
            metric = "ROC")

# model outputs & predictions
print(gb) ## ROC's accross tuning parameters
gb.pred <- predict(grids, gb, type = "prob") ## spatial predictions
stopCluster(mc)
# plot(gb.pred, axes=F)
fname <- paste("./Results/", labs, "_gb.rds", sep = "")
saveRDS(gb, fname)

# Neural network <nnet> ---------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)
tg <- expand.grid(size = seq(2,10, by=2), decay = c(0.001, 0.01, 0.1)) ## model tuning steps

# model training
nn <- train(fcal, lcal, 
            method = "nnet",
            preProc = c("center","scale"), 
            tuneGrid = tg,
            trControl = tc,
            metric ="ROC")

# model outputs & predictions
print(nn) ## ROC's accross tuning parameters
nn.pred <- predict(grids, nn, type = "prob") ## spatial predictions
stopCluster(mc)
# plot(nn.pred, axes=F)
fname <- paste("./Results/", labs, "_nn.rds", sep = "")
saveRDS(nn, fname)

# Model stacking setup ----------------------------------------------------
preds <- stack(gm.pred, gl1.pred, gl2.pred, rf.pred, gb.pred, nn.pred)
names(preds) <- c("gm","gl1","gl2","rf","gb","nn")
# plot(preds, axes = F)

# extract model predictions
coordinates(gs_val) <- ~x+y
projection(gs_val) <- projection(preds)
gspred <- extract(preds, gs_val)
gspred <- as.data.frame(cbind(gs_val, gspred))

# stacking model validation labels and features
gs_val <- as.data.frame(gs_val)
lval <- as.vector(t(gs_val[labs]))
fval <- gspred[,60:65] ## subset validation features

# Model stacking ----------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, 
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
si <- train(fval, lval,
            method = "glmStepAIC",
            family = "binomial",
            metric = "ROC",
            trControl = tc)

# model outputs & predictions
summary(si)
print(si)
si.pred <- predict(preds, si, type = "prob") ## spatial predictions
plot(si.pred, axes = F)
stopCluster(mc)
fname <- paste("./Results/", labs, "_si.rds", sep = "")
saveRDS(si, fname)

# Receiver-operator characteristics ---------------------------------------
cp_pre <- predict(si, fval, type="prob")
cp_val <- cbind(lval, cp_pre)
cpa <- subset(cp_val, cp_val=="A", select=c(A))
cpb <- subset(cp_val, cp_val=="B", select=c(A))
cp_eval <- evaluate(p=cpa[,1], a=cpb[,1]) ## calculate ROC's on test set
plot(cp_eval, 'ROC') ## plot ROC curve

# Generate feature mask ---------------------------------------------------
t <- threshold(cp_eval) ## calculate thresholds based on ROC
r <- matrix(c(0, t[,1], 0, t[,1], 1, 1), ncol=3, byrow = T) ## set threshold value <kappa>
mask <- reclassify(si.pred, r) ## reclassify stacked predictions

# Write prediction grids --------------------------------------------------
gspreds <- stack(preds, si.pred, mask)
names(gspreds) <- c("gm","gl1","gl2","rf","gb","nn","si","mk")
fname <- paste("./Results/","OAF_", labs, "_preds_2020.tif", sep = "")
writeRaster(gspreds, filename=fname, datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Site index prediction check ---------------------------------------------
coordinates(gsdat) <- ~x+y
projection(gsdat) <- projection(grids)
gspre <- extract(gspreds, gsdat)
gsout <- as.data.frame(cbind(gsdat, gspre))
gsout$mzone <- as.factor(ifelse(gsout$mk == 1, "A", "B"))
confusionMatrix(gsout$mzone, gsout$qy) ## overall prediction accuracy stats
boxplot(yield~mzone, notch=T, xlab="SI zone", ylab="Measured yield (t/ha)",
        cex.lab=1.3, gsout) ## yield differences between predicted site index zones

# Maize yield potentials (t/ha) ------------------------------------------
yld.lme <- lmer(log(yield)~factor(trt)*log(si/(1-si))+log(can+1)*log(dap+1)+(1|year)+(1|location), data = gsout)
summary(yld.lme) ## mixed model yield estimate results
gsout$yldf <- exp(fitted(yld.lme, gsout))

# Quantile regression (uncertainty) plot
par(pty="s")
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(yield~yldf, xlab="Production function (t/ha)", ylab="Measured yield (t/ha)", cex.lab=1.3, 
     xlim=c(0,8), ylim=c(0,8), gsout)
stQ <- rq(yield~yldf, tau=c(0.025,0.5,0.975), data=gsout)
print(stQ)
curve(stQ$coefficients[2]*x+stQ$coefficients[1], add=T, from=0, to=8, col="blue", lwd=2)
curve(stQ$coefficients[4]*x+stQ$coefficients[3], add=T, from=0, to=8, col="red", lwd=2)
curve(stQ$coefficients[6]*x+stQ$coefficients[5], add=T, from=0, to=8, col="blue", lwd=2)
abline(c(0,1), col="grey", lwd=1)

# Write output data frame -------------------------------------------------
fname <- paste("./Results/","OAF_", labs, "_out.csv", sep = "")
write.csv(gsout, fname, row.names = F)

library(leaflet)
library(sp)

sill <- projectRaster(si.pred, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

si <- leaflet() %>% 
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 8) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addRasterImage(sill, colors = pal, opacity = 0.5) %>%
  addLegend(colors=pal, values = values(sill), title = "SI")
si ## plot widget 