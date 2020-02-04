# TAMASA maize fertilizer nutrient ommission trials, Central Ethiopia, 2015/16
# Yield trial data courtesy of CIMMYT
# M. Walsh, J. Chamberlin, January 2019

# Required packages
# install.packages(c("downloader","rgdal","arm","MASS","quantreg")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(arm)
  require(MASS)
  require(quantreg)
})

# Create a data folder in your current working directory
dir.create("TAMASA_NOT_data", showWarnings=F)
setwd("./TAMASA_NOT_data")
dir.create("Results", showWarnings = F)
rm(list = ls()) ## scrub extraneous objects in memory

# Data downloads -----------------------------------------------------------
# download TAMASA yield data
download("https://www.dropbox.com/s/0cefqlb0g584mfd/ET_no_trials.zip?raw=1", "ET_no_trials.zip", mode = "wb")
unzip("ET_no_trials.zip", overwrite = T)
sites <- read.table("sites.csv", header=T, sep=",")
trial <- read.table("trials.csv", header=T, sep=",")
tresp <- merge(sites, trial, by="tid")

# download GADM-L3 shapefile (courtesy: http://www.gadm.org)
download("https://www.dropbox.com/s/25kw13359f8l5tr/ET_adm_shp.zip?raw=1", "ET_adm_shp.zip", mode = "wb")
unzip("ET_adm_shp.zip", overwrite = T)
shape <- shapefile("ETH_adm3.shp")

# Data setup --------------------------------------------------------------
# attach GADM-L3 admin unit names from shape
coordinates(tresp) <- ~lon+lat
projection(tresp) <- projection(shape)
gadm <- tresp %over% shape
tresp <- as.data.frame(tresp)
tresp <- cbind(gadm[ ,c(5,7,9)], tresp)
colnames(tresp) <- c("region","zone","woreda","tid","sid","lon","lat","cyld","year","loc","trt","tyld")
tresp$resp <- tresp$tyld-tresp$cyld
rlevel <- 2 ## set intended response-level relative to current untreated yield (2 = double control yield)
tresp$trti <- ifelse(tresp$tyld > rlevel*tresp$cyld, 1, 0)

# Regressions -------------------------------------------------------------
# Quantile regression, control vs NOT treatment yields
par(pty="s")
plot(tyld~cyld, xlab="Control yield (kg/ha)", ylab="Treatment yield (kg/ha)", xlim=c(-5,10005), cex.lab=1.3, tresp)
yq <- rq(log(tyld)~log(cyld), tau=c(0.25,0.5,0.75), data=tresp)
print(yq)
curve(exp(yq$coefficients[1])*x^yq$coefficients[2], add=T, from=0, to=10000, col="blue", lwd=2)
curve(exp(yq$coefficients[3])*x^yq$coefficients[4], add=T, from=0, to=10000, col="red", lwd=2)
curve(exp(yq$coefficients[5])*x^yq$coefficients[6], add=T, from=0, to=10000, col="blue", lwd=2)

# GLM yield response probability
rp <- glm(trti~trt+log(cyld), family=binomial(link="logit"), data=tresp)
display(rp)

# GLMER yield response probability ... site (response) index
rm <- glmer(trti~trt+log(cyld)+(1|sid), family=binomial(link="logit"), data=tresp)
display(rm)
ranef(rm)



