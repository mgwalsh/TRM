#+ Example maize yield site index spatial variance components from omission trial data.
#+ AfSIS ommission trial sample data, courtesy of CIAT
#+ M. Walsh, May 2015

#+ Required packages
# install.packages(c("downloader","proj4","arm")), dependencies=TRUE)
require(downloader)
require(proj4)
require(arm)

#+ Data download ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("OT_data", showWarnings=F)
dat_dir <- "./OT_data"

# AfSIS omission trial data download to ... "./OT_data"
download("https://www.dropbox.com/s/rs68t4077iafo17/AfSIS_OT.zip?dl=0", "./OT_data/AfSIS_OT.zip", mode="wb")
unzip("./OT_data/AfSIS_OT.zip", exdir="./OT_data", overwrite=T)
loc <- read.table(paste(dat_dir, "/Location.csv", sep=""), header=T, sep=",")
tri <- read.table(paste(dat_dir, "/Trial.csv", sep=""), header=T, sep=",")

