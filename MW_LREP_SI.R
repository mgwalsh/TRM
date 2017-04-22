#' Site-level Maize yield responses to fertilizer applications
#' Malawi LREP response trial data (courtesy of LREP)
#' LREP data documentation at: https://www.dropbox.com/s/4qbxnz4mdl92pdv/Malawi%20area-specific%20fertilizer%20recs%20report.pdf?dl=0
#' M. Walsh, September 2014

# Required packages
# install.packages(c("downloader","rgdal","arm")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(arm)
})

# Data download -----------------------------------------------------------
# Create a "LREP_data" folder in your current working directory
dir.create("LREP_data", showWarnings=F)
setwd("./LREP_data")

# LREP fertilizer response data download to "./Data"
download("https://www.dropbox.com/s/i4dby04fl9j042a/MW_fert_trials.zip?raw=1", "MW_fert_trials.zip", mode="wb")
unzip("MW_fert_trials.zip", overwrite=T)
mwsite <- read.table("Location.csv", header=T, sep=",")
mtrial <- read.table("Trial.csv", header=T, sep=",")

# Georeference and specify site ID's --------------------------------------
# Project to Africa LAEA from UTM36S
mw <- SpatialPoints(cbind(mwsite$Easting, mwsite$Northing), proj4string = CRS("+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs"))
tr <- as.data.frame(spTransform(mw, '+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs'))
colnames(tr) <- c("x","y")
mwsite <- cbind(mwsite, tr)

# Define unique grid cell / site ID's (GID)
# Specify pixel scale (res.pixel, in m)
res.pixel <- 1000

# Grid ID (GID) definition
xgid <- ceiling(abs(mwsite$x)/res.pixel)
ygid <- ceiling(abs(mwsite$y)/res.pixel)
gidx <- ifelse(mwsite$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(mwsite$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
mwsite.gid <- cbind(mwsite, GID)

# Merge location w trial data
mwresp <- merge(mwsite.gid, mtrial, by="LID")

# Exploratory plots -------------------------------------------------------
# ECDF plots of control and treatment yields
trt1 <- subset(mwresp, NPS==1 & Urea==1, select=c(Yt,Yc)) 
trt2 <- subset(mwresp, NPS==2 & Urea==2, select=c(Yt,Yc)) 
trt3 <- subset(mwresp, NPS==2 & Urea==3, select=c(Yt,Yc))
plot(ecdf(mwresp$Yc), main="", xlab="Maize yield (kg/ha)", ylab="Cum. proportion of observations", xlim=c(-50, 8050), verticals=T, lty=1, lwd=2, col="red", do.points=F)
abline(0.5,0, lty=2, col="grey")
plot(ecdf(trt1$Yt), add=T, verticals=T, lty=1, lwd=1, col="grey", do.points=F)
plot(ecdf(trt2$Yt), add=T, verticals=T, lty=1, lwd=1, col="grey", do.points=F)
plot(ecdf(trt3$Yt), add=T, verticals=T, lty=1, lwd=1, col="grey", do.points=F)

# ECDF plots of treatment response ratios
# plot(ecdf(log(trt1$Yt/trt1$Yc)), main="", verticals=T, lty=1, lwd=1, xlim=c(-0.5,3), xlab="Treatment response ratio = log(Yt/Yc)", ylab="Cum. proportion of observations", do.points=F)
# abline(0.5,0, lty=2)
# plot(ecdf(log(trt2$Yt/trt2$Yc)), add=T, main="", verticals=T, lty=1, lwd=1, do.points=F)
# plot(ecdf(log(trt3$Yt/trt3$Yc)), add=T, main="", verticals=T, lty=1, lwd=1, do.points=F)

# Response ratio plot
# plot(log(Yt/Yc)~log(Yc), ylab="Treatment response ratio = log(Yt/Yc)", xlab="Unfertilized control yield = log(Yc)", mwresp)
# abline(0,0, lwd=2, col="red")
# abline(log(2),0, col="grey")

# REML models -------------------------------------------------------------
# Treatment response ratio models
mlm1 <- lmer(log(Yt/Yc)~log(Yc)+NPS+Urea+(1|GID)+(1|Year/GID), data=mwresp)
display(mlm1)
mlm2 <- lmer(log(Yt/Yc)~log(Yc)+NPS+Urea+log(Yc)*NPS+log(Yc)*Urea+(1|GID)+(1|Year/GID), data=mwresp)
display(mlm2)
anova(mlm1, mlm2)

# Diagnostic plots of mlm2 model fit & residuals
plot(log(Yt/Yc)~fitted(mlm2), xlim=c(-2,8), ylim=c(-2,8), xlab="Modeled log(Yt/Yc)", ylab="Observed log(Yt/Yc)", mwresp)
abline(0,1, col="red")
# plot(residuals(mlm2)~fitted(mlm2), xlim=c(-2,8), ylim=c(-2,2), xlab="Modeled RR", ylab="Model residuals", mwresp)

# Not run: Conditional odds model of doubling yield
# mlm3 <- glmer(I(log(Yt/Yc)>log(2))~log(Yc)+NPS+Urea+log(Yc)*NPS+log(Yc)*Urea+(1|GID)+(1|Year/GID), family=binomial(link="logit"), data=mwresp)
# display(mlm3)

# Treatment yield model
mlm4 <- lmer(log(Yt)~log(Yc)+NPS+Urea+log(Yc)*NPS+log(Yc)*Urea+(1|GID)+(1|Year/GID), data=mwresp)
display(mlm4)
plot(log(Yt)~fitted(mlm4), xlim=c(4,10), ylim=c(4,10), xlab="Modeled log(Yt)", ylab="Observed log(Yt)", mwresp)
abline(0,1, col="red")

# Extract control yield and site response ratio indices at GID's ----------
mlm2.ran <- ranef(mlm2)
gidsrr <- as.data.frame(rownames(mlm2.ran$GID))
colnames(gidsrr) <- c("GID")
x <- aggregate(mwsite.gid$x, by=list(mwsite.gid$GID), FUN="mean")
colnames(x) <- c("GID", "Easting")
y <- aggregate(mwsite.gid$y, by=list(mwsite.gid$GID), FUN="mean")
colnames(y) <- c("GID", "Northing")
Yc <- aggregate(mwresp$Yc, by=list(mwresp$GID), FUN="mean")
colnames(Yc) <- c("GID", "Yc")
gidsrr <- merge(gidsrr, x, by="GID")
gidsrr <- merge(gidsrr, y, by="GID")
gidsrr <- merge(gidsrr, Yc, by="GID")
gidsrr$SRI <- mlm2.ran$GID[,1]

# Write site index file ---------------------------------------------------
write.csv(gidsrr, "MW_Site_Indices.csv", row.names=F)
