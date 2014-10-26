# Site-level Maize yield responses to fertilizer applications
# Malawi LREP response trial data (courtesy of LREP & Todd Benson)
# LREP data documentation at: https://www.dropbox.com/s/4qbxnz4mdl92pdv/Malawi%20area-specific%20fertilizer%20recs%20report.pdf?dl=0
# M. Walsh, September 2014

# Set local working directory e.g.
# dat_dir <- "/Users/markuswalsh/Documents/LDSF/Malawi/Fert_resp_models/data"
# setwd(dat_dir)

# Required packages
# install.packages(c("downloader","proj4","raster","arm", "gstat")), dependencies=TRUE)
require(downloader)
require(proj4)
require(raster)
require(arm)
require(gstat)

# Response trial data ------------------------------------------------------

download("https://www.dropbox.com/s/243n844p9kep3e6/MW_fert_trials.zip?dl=0", "MW_fert_trials.zip", mode="wb")
unzip("MW_fert_trials.zip")
mtrial <- read.table("Trial.csv", header=T, sep=",")
mwsite <- read.table("Location.csv", header=T, sep=",")

# Define "coordinate reference system" (CRS)
# Project to Africa LAEA from UTM36S
mw <- cbind(mwsite$Easting, mwsite$Northing)
tr <- ptransform(mw, '+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs', '+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs')
colnames(tr) <- c("x","y","z")
mwsite <- cbind(mwsite, tr)

# Specify grid cell ID's (GID)
# Define pixel scale (res.pixel, in m)
res.pixel <- 1000

# Grid ID (GID) definition
xgid <- ceiling(abs(mwsite$x)/res.pixel)
ygid <- ceiling(abs(mwsite$y)/res.pixel)
gidx <- ifelse(mwsite$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(mwsite$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
mwsite.gid <- cbind(mwsite, GID)

# Grid cell center point locations
mwsite.gid$X <- (ceiling(mwsite.gid$x/res.pixel)*1000)-0.5*res.pixel
mwsite.gid$Y <- (ceiling(mwsite.gid$y/res.pixel)*1000)-0.5*res.pixel

# Merge location w trial data
mwresp <- merge(mwsite.gid, mtrial, by="LID")

# ECDF plots of control and treatment yields
trt1 <- subset(mwresp, NPS==1 & Urea==1, select=c(Yt,Yc)) 
trt2 <- subset(mwresp, NPS==2 & Urea==2, select=c(Yt,Yc)) 
trt3 <- subset(mwresp, NPS==2 & Urea==3, select=c(Yt,Yc))
plot(ecdf(mwresp$Yc), main="", xlab="Maize yield (kg/ha)", ylab="Cum. proportion of observations", xlim=c(-50, 8050), verticals=TRUE, lty=1, lwd=2, col="red", do.points=FALSE)
abline(0.5,0, lty=2, col="grey")
plot(ecdf(trt1$Yt), add=T, verticals=TRUE, lty=1, lwd=1, col="grey", do.points=FALSE)
plot(ecdf(trt2$Yt), add=T, verticals=TRUE, lty=1, lwd=1, col="grey", do.points=FALSE)
plot(ecdf(trt3$Yt), add=T, verticals=TRUE, lty=1, lwd=1, col="grey", do.points=FALSE)

# ECDF plots of treatment response ratios
plot(ecdf(log(trt1$Yt/trt1$Yc)), main="", verticals=TRUE, lty=1, lwd=1, xlim=c(-0.5,3), xlab="Treatment response ratio = log(Yt/Yc)", ylab="Cum. proportion of observations", do.points=FALSE)
abline(0.5,0, lty=2)
plot(ecdf(log(trt2$Yt/trt2$Yc)), add=T, main="", verticals=TRUE, lty=1, lwd=1, do.points=FALSE)
plot(ecdf(log(trt3$Yt/trt3$Yc)), add=T, main="", verticals=TRUE, lty=1, lwd=1, do.points=FALSE)

# Response ratio plot
plot(log(Yt/Yc)~log(Yc), ylab="Treatment response ratio = log(Yijk/Y0jk)", xlab="Unfertilized control yield = log(Y0jk)", mwresp)
abline(0,0, lwd=2, col="red")
abline(log(2),0, col="grey")

# REML models -------------------------------------------------------------

mlm1 <- lmer(log(Yt/Yc)~log(Yc)+NPS+Urea+(1|GID)+(1|Year/GID), data=mwresp)
display(mlm1)
mlm2 <- lmer(log(Yt/Yc)~log(Yc)+NPS+Urea+log(Yc)*NPS+log(Yc)*Urea+(1|GID)+(1|Year/GID), data=mwresp)
summary(mlm2)
anova(mlm1, mlm2)

# Aside: Conditional odds model of doubling yield
# mlm3 <- glmer(I(log(Yt/Yc)>log(2))~log(Yc)+NPS+Urea+log(Yc)*NPS+log(Yc)*Urea+(1|GID)+(1|Year/GID), family=binomial(link="logit"), data=mwresp)
# display(mlm3)

# Diagnostic plots of mlm2 model fit & residuals
plot(log(Yt/Yc)~fitted(mlm2), xlim=c(-2,8), ylim=c(-2,8), xlab="Modeled log(Yt/Yc)", ylab="Observed log(Yt/Yc)", mwresp)
abline(0,1, col="red")
# plot(residuals(mlm2)~fitted(mlm2), xlim=c(-2,8), ylim=c(-2,2), xlab="Modeled RR", ylab="Model residuals", mwresp)

# Extract mean control yields and site response ratios at GID's  ----------

mlm2.ran <- ranef(mlm2)
gidsrr <- as.data.frame(rownames(mlm2.ran$GID))
colnames(gidsrr) <- c("GID")
x <- aggregate(mwsite.gid$X, by=list(mwsite.gid$GID), FUN="mean")
colnames(x) <- c("GID", "Easting")
y <- aggregate(mwsite.gid$Y, by=list(mwsite.gid$GID), FUN="mean")
colnames(y) <- c("GID", "Northing")
Yc <- aggregate(mwresp$Yc, by=list(mwresp$GID), FUN="mean")
colnames(Yc) <- c("GID", "Yc")
gidsrr <- merge(gidsrr, x, by="GID")
gidsrr <- merge(gidsrr, y, by="GID")
gidsrr <- merge(gidsrr, Yc, by="GID")
gidsrr$SRR <- mlm2.ran$GID[,1]

# ECDF plot of GID-level Site Response Ratios
plot(ecdf(gidsrr$SRR), main="", verticals=TRUE, col="red", xlab="Site Response Ratio", ylab="Cum. proportion of observations", lty=1, lwd=1, do.points=FALSE)

# Overlay gridded covariates ----------------------------------------------

# Malawi grids download (~7.2 Mb)
download("https://www.dropbox.com/s/54di5f37yp30bz4/MW_grids.zip?dl=0", "MW_grids.zip", mode="wb")
unzip("MW_grids.zip", overwrite=T)

# Grid overlay
coordinates(gidsrr) <- ~Easting+Northing
proj4string(gidsrr) <- CRS("+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=5 +lon_0=20 +no_defs")

grid.list <- c("BSANs.tif","BSASs.tif","BSAVs.tif","CTIs.tif","ELEVs.tif","EVIs.tif","LSTDs.tif","LSTNs.tif","REF1s.tif","REF2s.tif","REF3s.tif","REF7s.tif","RELIs.tif","TMAPs.tif","TMFIs.tif")
for (i in 1:length(grid.list)){
  print(paste("extracting", grid.list[i]))
  grid.cov <- raster(grid.list[i]) 
  gidsrr@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract(
    x = grid.cov, 
    y = gidsrr,
    method = "simple")
}
MW_SI <- as.data.frame(gidsrr)
mwgrids <- stack(grid.list)
# plot(mwgrids)

# Write csv
write.csv(MW_SI, "MW_SI.csv")

# Control yield & SRR variograms ------------------------------------------

MW_SI$x <- MW_SI$x/1000
MW_SI$y <- MW_SI$y/1000
coordinates(MW_SI) <- ~x+y

# Mean GID-level control yield variogram 
yc.var <- variogram(I(Yc/1000)~1, MW_SI, cutoff=50)
yc.fit <- fit.variogram(yc.var, model = vgm(1, "Sph", 50, 1))
plot(yc.var, yc.fit)

# Mean GID-level response ratio variogram
srr.var <- variogram(SRR~1, MW_SI, cutoff=50)
srr.fit <- fit.variogram(srr.var, model = vgm(1, "Sph", 50, 1))
plot(srr.var, srr.fit)
