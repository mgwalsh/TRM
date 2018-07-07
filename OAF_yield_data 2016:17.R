# OAF 2016/2017 Maize yield gap, Western Kenya data setup
# Yield data courtesy of One Acre Fund. See @ ??
# M. Walsh, July 2018

# Required packages
# install.packages(c("downloader","rgdal","raster","quantreg","leaflet","htmlwidgets")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
  require(quantreg)
  require(leaflet)
  require(htmlwidgets)
})

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("OAF_data", showWarnings=F)
setwd("./OAF_data")

# download OAF yield data
# download("", "", mode = "wb") ## insert OAF data link here
# unzip("oafdata.csv.zip", overwrite = T)
yield <- read.table("oafyga.csv", header = T, sep = ",")
yield <- yield[!duplicated(yield), ] ## removes duplicates

# download GADM-L3 shapefile (courtesy: http://www.gadm.org)
download("https://www.dropbox.com/s/otspr9b9jtuyneh/KEN_adm3.zip?raw=1", "KEN_adm3.zip", mode = "wb")
unzip("KEN_adm3.zip", overwrite = T)
shape <- shapefile("KEN_adm3.zip")

# download raster stack
download("https://www.dropbox.com/s/mz1t0zyq8uoqrhq/KE_250m_2017.zip?raw=1", "KE_250m_2017.zip", mode = "wb")
unzip("KE_250m_2017.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# Data setup --------------------------------------------------------------
# attach GADM-L3 admin unit names from shape
coordinates(yield) <- ~lon+lat
projection(yield) <- projection(shape)
gadm <- yield %over% shape
yield <- as.data.frame(yield)
yield <- cbind(gadm[ ,c(5,7,9)], yield)
colnames(yield) <- c("district","division","location","id","lat","lon","year","trt","can","dap","fsize","yield")

# project survey coords to grid CRS
yield.proj <- as.data.frame(project(cbind(yield$lon, yield$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(yield.proj) <- c("x","y")
yield <- cbind(yield, yield.proj)
coordinates(yield) <- ~x+y
projection(yield) <- projection(yield)

# extract gridded variables at survey locations
yieldgrid <- extract(grids, yield)
gsdat <- as.data.frame(cbind(yield, yieldgrid)) 
gsdat <- gsdat[complete.cases(gsdat[,c(1:3,13:44)]),] ## removes incomplete cases
gsdat <- gsdat[which(gsdat$can < 100 & gsdat$dap < 100), ] ## removes outlier fertilizer treatments
gsdat <- gsdat[which(gsdat$fsize > 0), ] ## removes field size = 0

# Classify yield propensities by conditional quantile ---------------------
# this is the conditional yield gap based on the current data at median values
qy.rq <- rq(log(yield)~factor(year)+factor(trt)+I(dap/fsize)*I(can/fsize), tau = 0.5, data = gsdat) ## try quantiles other than the median
summary(qy.rq)
gsdat$qy <- as.factor(ifelse(exp(predict(qy.rq, gsdat)) > gsdat$yield, "B", "A"))
table(gsdat$qy)
boxplot(yield~qy, notch=T, gsdat)

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/OAF_gsdat.csv", row.names = F)

# Yield survey map widget -------------------------------------------------
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'OAF_yield_survey.html', selfcontained = T) ## save widget
