#' Quantile regressions of Maize yield responses to fertilizer applications
#' Malawi LREP response trial data (courtesy of LREP)
#' LREP data documentation at: https://www.dropbox.com/s/4qbxnz4mdl92pdv/Malawi%20area-specific%20fertilizer%20recs%20report.pdf?dl=0
#' M. Walsh, December 2015

# Required packages
# install.packages(c("downloader","quantreg")), dependencies=TRUE)
require(downloader)
require(quantreg)

# Data setup --------------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("MW_data", showWarnings=F)
setwd("./MW_data")

# download LREP fertilizer response data
download("https://www.dropbox.com/s/i4dby04fl9j042a/MW_fert_trials.zip?dl=0", "MW_fert_trials.zip", mode="wb")
unzip("MW_fert_trials.zip", overwrite=T)
sites <- read.table("Location.csv", header=T, sep=",")
trial <- read.table("Trial.csv", header=T, sep=",")
mresp <- merge(sites, trial, by="LID")
mresp <- mresp[order(mresp$Yt),] ## order dataframe based on treated yield (Yt)
mresp$Year <- mresp$Year-1996

# Exploratory plots -------------------------------------------------------
# Treatment/Control plot
plot(Yt ~ Yc, data = mresp, cex= .5, col = "grey", 
     xlim = c(-200, 8200), ylim = c(-200, 8200),
     xlab = "Unfertilized yield (kg/ha)", ylab = "Fertilized yield (kg/ha)")
abline(c(0,1), col = "red")
tau <- c(.025,.5,.975)
for(i in 1:length(tau)) {
  abline(rq(Yt~Yc, tau=tau[i], data = mresp), col = "blue", lty = 2)
}

# ECDF plot
trt1 <- subset(mresp, NPS==1 & Urea==1, select=c(Yt,Yc)) 
trt2 <- subset(mresp, NPS==2 & Urea==2, select=c(Yt,Yc)) 
trt3 <- subset(mresp, NPS==2 & Urea==3, select=c(Yt,Yc))
plot(ecdf(mresp$Yc), main="", xlab="Maize yield (kg/ha)", ylab="Cum. proportion of observations", xlim=c(-50, 8050), verticals=T, lty=1, lwd=2, col="red", do.points=F)
abline(0.5,0, lty=2, col="grey")
plot(ecdf(trt1$Yt), add=T, verticals=T, lty=1, lwd=1, col="grey", do.points=F)
plot(ecdf(trt2$Yt), add=T, verticals=T, lty=1, lwd=1, col="grey", do.points=F)
plot(ecdf(trt3$Yt), add=T, verticals=T, lty=1, lwd=1, col="grey", do.points=F)

# Quantile regression -----------------------------------------------------
attach(mresp)
Yt.rq <- rq(Yt~Yc+NPS+Urea, tau = seq(0.05, 0.95, by = 0.05), data = mresp)
detach(mresp)

# Result plots
plot(summary(Yt.rq), main = c("Intercept","Unfertilized yield","NPS","Urea"))
