# OAF 2016/2017 Maize yields, Western Kenya data setup
# Yield data courtesy of One Acre Fund.
# M. Walsh, July 2020

# Required packages
# install.packages(c("downloader","rgdal","raster","quantreg","arm","leaflet","htmlwidgets")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(sp)
  require(raster)
  require(quantreg)
  require(arm)
  require(leaflet)
  require(htmlwidgets)
})

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("OAF_data", showWarnings=F)
setwd("./OAF_data")

# download OAF yield data
download("https://www.dropbox.com/s/lhniws8prfvjw9r/OAF_yield_data.csv.zi?raw=1", "OAF_yield_data.csv.zip", mode = "wb")
unzip("OAF_yield_data.csv.zip", overwrite = T)
yield <- read.table("OAF_yield_data.csv", header = T, sep = ",")
yield$trt <- ifelse(yield$trt == 1, "oaf", "control")
# yield <- yield[!duplicated(yield), ] ## removes duplicates if needed

# download GADM-L3 shapefile (@ http://www.gadm.org)
download("https://www.dropbox.com/s/otspr9b9jtuyneh/KEN_adm3.zip?raw=1", "KEN_adm3.zip", mode = "wb")
unzip("KEN_adm3.zip", overwrite = T)
shape <- shapefile("KEN_adm3.zip")

# download raster stack
download("https://osf.io/4jvnu?raw=1", "KE_250m_2020.zip", mode = "wb")
unzip("KE_250m_2020.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# set ROI grid extent
ext <- data.frame(lat = c(-1.2,-1.2,1.2,1.2), lon = c(33.9,35.5,33.9,35.5)) ## set ROI extent in degrees
names(ext) <- c("lat","lon")
coordinates(ext) <- ~ lon + lat
proj4string(ext) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
ext <- spTransform(ext, CRS("+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
bb <- extent(ext)
grids <- crop(grids, bb)

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
quant <- quantile(gsdat$yield, probs=c(0.025,0.975))
gsdat <- gsdat[which(gsdat$yield > quant[1]), ] ## removes outlier low yields (0.025 quantile)
gsdat <- gsdat[which(gsdat$yield < quant[2]), ] ## removes outlier high yields (0.975 quantile)

# Fit production function -------------------------------------------------
# production function using quantile regression
qy.rq <- rq(log(yield)~year+factor(trt)+log(dap+1)*log(can+1), tau = 0.5, data = gsdat)
summary(qy.rq)
gsdat$qy <- as.factor(ifelse(exp(predict(qy.rq, gsdat)) > gsdat$yield, "B", "A")) ## classify site index zones
table(gsdat$qy)
boxplot(yield~qy, notch=T, gsdat)

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/OAF_gsdat.csv", row.names = F)

# Yield survey map widget -------------------------------------------------
w <- leaflet() %>%
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 8) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
saveWidget(w, 'OAF_yield_survey.html', selfcontained = T) ## save widget
w ## plot widget 

