# Stacked spatial predictions of 2016/2017 OAF maize yield gap potentials
# M. Walsh, July 2018

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
rm(list=setdiff(ls(), c("gsdat","grids"))) ## scrub extraneous objects in memory

# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(gsdat$qy, p = 4/5, list = F, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# Survey calibration labels
cp_cal <- gs_cal$qy

# raster calibration features
gf_cal <- gs_cal[,13:44]

# Central place theory model <glmnet> --------------------------------------
# select central place covariates
gf_cpv <- gs_cal[,18:26]

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
rr <- train(gf_cpv, cp_cal, 
            method = "glmnet",
            family = "binomial",
            preProc = c("center","scale"), 
            trControl = tc,
            metric ="ROC")

# model outputs & predictions
print(rr)
plot(varImp(rr))

stopCluster(mc)

# GLMNET with all covariates ----------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
rr1 <- train(gf_cal, cp_cal, 
             method = "glmnet",
             family = "binomial",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="ROC")

# model outputs & predictions
summary(rr1)
print(rr1) ## ROC's accross cross-validation
rr1.pred <- predict(grids, rr1, type = "prob") ## spatial predictions

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
rf <- train(gf_cal, cp_cal,
            preProc = c("center","scale"),
            method = "rf",
            ntree = 501,
            metric = "ROC",
            tuneGrid = tg,
            trControl = tc)

# model outputs & predictions
print(rf) ## ROC's accross tuning parameters
plot(varImp(rf)) ## relative variable importance
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
gb <- train(gf_cal, cp_cal, 
            method = "gbm", 
            preProc = c("center", "scale"),
            trControl = tc,
            tuneGrid = tg,
            metric = "ROC")

# model outputs & predictions
print(gb) ## ROC's accross tuning parameters
plot(varImp(gb)) ## relative variable importance
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
nn <- train(gf_cal, cp_cal, 
            method = "nnet",
            preProc = c("center","scale"), 
            tuneGrid = tg,
            trControl = tc,
            metric ="ROC")

# model outputs & predictions
print(nn) ## ROC's accross tuning parameters
plot(varImp(nn)) ## relative variable importance
nn.pred <- predict(grids, nn, type = "prob") ## spatial predictions

stopCluster(mc)

# Model stacking setup ----------------------------------------------------
preds <- stack(rr1.pred, rf.pred, gb.pred, nn.pred)
names(preds) <- c("rr1","rf", "gb","nn")
plot(preds, axes = F)

# extract model predictions
coordinates(gs_val) <- ~x+y
projection(gs_val) <- projection(preds)
gspred <- extract(preds, gs_val)
gspred <- as.data.frame(cbind(gs_val, gspred))

# stacking model validation labels and features
cp_val <- gspred$qy
gf_val <- gspred[,46:49] ## subset validation features

# Model stacking ----------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, 
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
st <- train(gf_val, cp_val,
            method = "glm",
            family = "binomial",
            metric = "ROC",
            trControl = tc)

# model outputs & predictions
print(st)
plot(varImp(st))
st.pred <- predict(preds, st, type = "prob") ## spatial predictions of maize yield propensities
plot(st.pred, axes = F)

stopCluster(mc)

# Receiver-operator characteristics ---------------------------------------
cp_pre <- predict(st, gf_val, type="prob")
cp_val <- cbind(cp_val, cp_pre)
cpp <- subset(cp_val, cp_val=="B", select=c(B))
cpa <- subset(cp_val, cp_val=="A", select=c(B))
cp_eval <- evaluate(p=cpp[,1], a=cpa[,1]) ## calculate ROC's on test set
plot(cp_eval, 'ROC') ## plot ROC curve

# Generate feature mask ---------------------------------------------------
t <- threshold(cp_eval) ## calculate thresholds based on ROC
r <- matrix(c(0, t[,1], 0, t[,1], 1, 1), ncol=3, byrow = T) ## set threshold value <kappa>
mask <- reclassify(st.pred, r) ## reclassify stacked predictions

# Write prediction grids --------------------------------------------------
gspreds <- stack(preds, st.pred, mask)
names(gspreds) <- c("rr1","rf","gb","nn","st","mk")
writeRaster(gspreds, filename="./Results/KE_preds_2017.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Write output data frame -------------------------------------------------
coordinates(gsdat) <- ~x+y
projection(gsdat) <- projection(grids)
gspre <- extract(gspreds, gsdat)
gsout <- as.data.frame(cbind(gsdat, gspre))

# prediction summaries
gsout$mzone <- ifelse(gsout$mk == 1, "A", "B")
boxplot(yield~mzone, notch=T, gsout)
table(gsout$mzone, gsout$qy)
write.csv(gsout, "./Results/OAF_preds_2017.csv", row.names = F)

# ECDF plot of predicted management zone maize yields
mzA <- subset(gsout, mzone=='A', select=yield) 
mzB <- subset(gsout, mzone=='B', select=yield) 
plot(ecdf(mzA$yield), verticals=T, lty=1, lwd=1, col="dark green", do.points=F, main="",
     xlab="Expected maize yield (Mg/ha)", ylab="Cum. proportion of observations")
plot(ecdf(mzB$yield), add=T, verticals=T, lty=1, lwd=1, col="red", do.points=F)
abline(0.5,0, lty=2, col="grey")

# Prediction map widget ---------------------------------------------------
pred <- st.pred ## management zone ensemble probability
pal <- colorBin("Greens", domain = 0:1) ## set color palette
w <- leaflet() %>% 
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 8) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addRasterImage(pred, colors = pal, opacity = 0.7, maxBytes=6000000) %>%
  addLegend(pal = pal, values = values(pred), title = "SI = p(A) ") %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'KE_high_SI_prob.html', selfcontained = T) ## save html ... change feature names here

