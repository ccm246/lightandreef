##- Miller, Fobert, & Rice 2023 Code -##
##- Light Pollution and Coral Reefs --##
##- C R Miller 2023 ------------------##


#### libraries ####

library(ggplot2)
library(tidyverse)
library('ggspatial')
library(cowplot)
library(ggrepel)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(googleway)
library(sf)
library(spData)
library(sjPlot)
library(mgcv)
library(MuMIn)
library(plot.gam)
library(gstat)
library(nlme)
library(sp)
library(parallel)
library(progress)

#### load in and set up data ####

reef <- read.csv('nc_reef_data_dhw.csv')

#create dark and light datasets
darkknown <- subset(reef, diff < 0 & measuredlight == 1) 
darkNA <- reef[reef$measuredlight == 0,]
darkall <- rbind(darkknown, darkNA)
light <- subset(reef, lightpres == 1)

head(darkall)
head(light)

#### find percentage of reefs experiencing light pollution at their own depth ####
darkknowncalc <- rbind(light, darkknown)
nrow(light)/nrow(reef) # % of US EEZ corals experience light pollution including uncalculated light areas
nrow(light)/nrow(darkknowncalc) # % of US EEZ corals experience light pollution only including calculated light areas
nrow(light[light$light>30,])/nrow(light) # % mesophotic reefs experiencing light pollution

#### find percentage of reefs experiencing noise pollution ####

head(reef)
nrow(reef[reef$vesseltotal > 0,])/nrow(reef)

#### find percentage of reefs experiencing warming ####

head(reef)
nrow(reef[reef$hsmax > 1,])/nrow(reef)

#### find percentage of reefs experiencing dhw warming ####

head(reef)
nrow(reef[reef$dhw > 1,])/nrow(reef)

#### find percentage of reefs experiencing noise, light, short term heat stress

head(reef)
nrow(reef[reef$hsmax > 1 & reef$lightpres ==1 & reef$vesseltotal > 0,])/nrow(reef)


#### clean up environment ####

darkall <- NULL
darkNA <- NULL
darkknown <- NULL
darkknowncalc <- NULL

reef$fid <- NULL
reef$sstmin <- NULL
reef$coord <- NULL

#### figure out correlation structure of data ####

#first set up initial model and take a look
model <- lm(hsmax ~ vesseltotal * diff * sstmax, reef)
summary(model)

#initial bubble plot
E <- rstandard(model)
mydata <- data.frame(E, reef$lon, reef$lat)
coordinates(mydata) <- c('reef.lon', 'reef.lat')

Variol <- variogram(E~1, data = mydata)
plot(Variol)
Variol <- NULL

#structure options
spA <- corSpher(form = ~ lon + lat, nugget = T)
spC <- corRatio(form = ~ lon + lat, nugget = T)
spD <- corGaus(form = ~ lon + lat, nugget = T)
spE <- corExp(form = ~ lon + lat, nugget = T)

#see which is best structure

#subsample
lightsample <- light[sample(nrow(light), 1000), ]

#formula
f1 <- formula(hsmax ~ vesseltotal * lightpres * sstmax)
#lmeControl(msTol = 1e-7)

#options
B1A <- gls(f1, correlation = spA, data = lightsample)
B1C <- gls(f1, correlation = spC, data = lightsample)
B1D <- gls(f1, correlation = spD, data = lightsample)
B1E <- gls(f1, correlation = spE, data = lightsample)
B1.gls <- gls(f1, data = lightsample)
AIC(B1.gls, B1C, B1D, B1E)

#id column
reef$fid <- 1:nrow(reef)

#cleanup
B1A <- NULL
B1C <- NULL
B1D <- NULL
B1E <- NULL
B1.gls <- NULL

light$fid <- NULL
light$sstmin <- NULL
light$coord <- NULL

#### run models ####

cl <- makeCluster(3)

#spatial correlation set up
spE <- corExp(form = ~ lon + lat, nugget = F)

#first create random dataset
randomcoral <- reef[sample(nrow(reef), 1000),]
randomcoral <- randomcoral[randomcoral$sstmax > 10,] #REMOVE OUTLIER
randomcoral$dhw_sqrt <- sqrt(randomcoral$dhw)

#run simple model with light only hs
model_light <- gls(hsmax ~ lightpres, correlation = spE, method='ML', data=randomcoral)
summary(model_light)
plot(model_light)

#run simple model with noise only hs
model_noise <- gls(hsmax ~ vesseltotal, correlation = spE, method='ML', data=randomcoral)
summary(model_noise)
plot(model_noise)

#run simple model with sstmax only hs
model_sst <- gls(hsmax ~ sstmax, correlation = spE, method='ML', data=randomcoral)
summary(model_sst)
plot(model_sst)

#run simple model with light only dhw
model_light <- gls(dhw_sqrt ~ lightpres, correlation = spE, method='ML', data=randomcoral)
summary(model_light)
plot(model_light)

#run simple model with noise only dhw
model_noise <- gls(dhw_sqrt ~ vesseltotal, correlation = spE, method='ML', data=randomcoral)
summary(model_noise)
plot(model_noise)

#run simple model with sstmax only dhw
model_sst <- gls(dhw_sqrt ~ sstmax, correlation = spE, method='ML', data=randomcoral)
summary(model_sst)
plot(model_sst)

#spatial correlation set up
spE <- corExp(form = ~ lon + lat, nugget = F)
spA <- corSpher(form = ~ lon + lat, nugget = F)

#hotspot model
hsmaxmodel <- gls(hsmax.x ~ scale(vesseltotal) * lightpres + scale(bath), correlation = spA, randomcoral)
summary(hsmaxmodel)

#degree heating model with sea surface temperature
dhwsst <- gls(dhw_sqrt ~ scale(vesseltotal) * lightpres * scale(sstmax) + scale(bath.x), correlation = spA, randomcoral)
summary(dhwsst)

#degree heating model 
dhwmod <- gls(dhw_sqrt ~ scale(vesseltotal) * lightpres + scale(bath), correlation = spA, randomcoral)
summary(dhwmod)

#plot

tiff('figure2.tiff', units="in", width=5, height=4, res=600, compression = 'lzw')

p <- ggplot(randomcoral, aes(x = sstmax, y = hsmax, colour = lightpres)) + 
  geom_point() + labs(x='Maximum Sea Surface Temperature (\u00B0C)', y='Hot Spot Likelihood ((\u00B0C) of Heat Stress)') +
  geom_smooth(method = 'lm', fullrange=T, fill='gray') + theme_classic() + 
  scale_color_manual(values=c(alpha('darkslategray',.3), 'goldenrod1')) +
  theme(text=element_text(size=10)) 
p + labs(color='Light Pollution \n Presence') 

dev.off()



