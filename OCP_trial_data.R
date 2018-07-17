# OCP/IITA 2017 Maize fertilizer trials, Central Nigeria
# Yield response data courtesy of IITA
# J. Huising and M. Walsh, July 2018

# Required packages
# install.packages(c("downloader","rgdal","raster","quantreg","leaflet","htmlwidgets")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(sp)
  require(raster)
  require(quantreg)
  require(leaflet)
  require(htmlwidgets)
})

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("OCP_data", showWarnings=F)
setwd("./OCP_data")

# download IITA/OCP yield data
download("https://www.dropbox.com/s/efr02hlxn3n1yvn/OCP_trial_data.zip?raw=1", "OCP_trial_data.zip", mode = "wb")
unzip("OCP_trial_data.zip", overwrite = T)
sites <- read.table("location.csv", header=T, sep=",")
trial <- read.table("trial.csv", header=T, sep=",")
tresp <- merge(sites, trial, by="tid")
# boxplot(Yc~trt, tresp, notch=T)
# boxplot(Yo~trt, tresp, notch=T)
# plot(Yo~Yc, tresp)

