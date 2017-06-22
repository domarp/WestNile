##=========================================================================================
## Script for Data Preparation
##=========================================================================================

rm(list = ls())

library(dplyr)
library(lubridate)
library(data.table)
library(reshape2)
library(geosphere)

##=========================================================================================
## Set working directory
##=========================================================================================
setwd("C:\\Users\\prbalakrishnan\\Desktop\\New folder\\")
source("eda_utilities.r")

##=========================================================================================
## read Files
##=========================================================================================
train.data = read.csv("train.csv")
spray.data = read.csv("spray.csv")
weather.data = read.csv("weather.csv", na.strings = c(" ","M","-"))
test.data = read.csv("test.csv")

##=========================================================================================
## Fix train data to remove duplicate rows (Date, Trap, Species)
## NumMosquitos = sum(NumMosquitos)
## WnVPresent = max(WnvPresent)
## WnVTotal = sum(WnvTotal)
##=========================================================================================
train.data = train.data %>% group_by(Date, Trap, Latitude, Longitude, Species) %>% 
  summarize(Address = last(Address), Block = last(Block), Street = last(Street), 
            AddressNumberAndStreet = last(AddressNumberAndStreet), 
            AddressAccuracy = mean(AddressAccuracy), NumMosquitos = sum(NumMosquitos), 
            WnvPresent = max(WnvPresent))#,WnvTotal = sum(WnvPresent))

##=========================================================================================
## Append train & test into one dataset
##=========================================================================================

train.data$Id <- -999
test.data$WnvPresent <- NA
test.data$NumMosquitos <- -999
test.data <- test.data[,names(train.data)]
# test.data$WnvTotal <- NA
 
# test.data = test.data %>% group_by(Date, Trap, Latitude, Longitude, Species) %>% 
#   summarize(Address = last(Address), Block = last(Block), Street = last(Street), 
#             AddressNumberAndStreet = last(AddressNumberAndStreet), 
#             AddressAccuracy = mean(AddressAccuracy), NumMosquitos = sum(NumMosquitos), 
#             WnvPresent = max(WnvPresent))#,WnvTotal = sum(WnvPresent))

tags.train.data <- rbind(data.frame(train.data), data.frame(test.data))

##=========================================================================================
## Prep Data for merging & merge
##=========================================================================================

# Create Zipcode variable & Remove redudant columns in train data
ind <- regexpr("IL",tags.train.data$Address)
tags.train.data$Zipcode <- as.numeric(substr(tags.train.data$Address, ind + 3, ind +7))
tags.train.data$Zipcode <- as.character(tags.train.data$Zipcode)
drop.columns.trn <- c("Address","AddressNumberAndStreet","Street")
tags.train.data <- tags.train.data[,!(names(tags.train.data) %in% drop.columns.trn)]

# Remove redudant columns in weather data
drop.columns.weather <- c("CodeSum","Water1","SnowFall", "Depart", "Sunrise", "Sunset", "Depth")
weather.data <- weather.data[,!(names(weather.data) %in% drop.columns.weather)]

# Clean weather data preciptotal for 'T'
weather.data$PrecipTotal = as.character(weather.data$PrecipTotal)
weather.data$PrecipTotal[which(trimws(weather.data$PrecipTotal) %in% c('T',' T'))] = 0.0
weather.data$PrecipTotal = as.numeric(weather.data$PrecipTotal)

# Split weather data by weather stations
weather.data.1 <- subset(weather.data, Station == 1)
names(weather.data.1) <- sapply(names(weather.data.1), function(x) paste0(x, ".1"))
names(weather.data.1)[2] = "Date"
weather.data.2 <- subset(weather.data, Station == 2)
names(weather.data.2) <- sapply(names(weather.data.2), function(x) paste0(x, ".2"))
names(weather.data.2)[2] = "Date"

# Merge
model.data <- merge(tags.train.data, weather.data.1, on = c("Date"))
model.data <- merge(model.data, weather.data.2, on = c("Date"))

# Quick EDA check
df.stats.num <- statsNum(model.data)
df.stats.nnum <- statsNonNum(model.data)

##=========================================================================================
## Derive Variables - 1
##=========================================================================================

# Date Variable
model.data$Month <- month(model.data$Date)
model.data$Year <- year(model.data$Date)
model.data$day <- yday(model.data$Date)
model.data$WoY <- as.integer(yday(model.data$Date)/7)
model.data$DoW <- wday(model.data$Date)

# Trap variable
model.data$isSat <- as.integer(nchar(as.character(model.data$Trap)) != 4)
model.data$orig.Trap <- model.data$Trap
model.data$orig.Trap[which(model.data$isSat == 1)] <- 
  substr(model.data$orig.Trap[which(model.data$isSat == 1)], 1,4)

# Refine Species
model.data$Species.1 <- as.character(model.data$Species)
valid.Species <- c("CULEX PIPIENS","CULEX PIPIENS/RESTUANS","CULEX RESTUANS")
model.data$Species.1[which(!(model.data$Species.1 %in% valid.Species))] = "MISC"
model.data$Species <- factor(model.data$Species.1)
model.data$Species.1 <- NULL

# NumMosquitos
df.num.mosquito <- subset(model.data, Id == -999) %>% group_by(Year, orig.Trap, Species,Month) %>% 
  summarize(mean.mos.count = mean(NumMosquitos))
df.num.mosquito$Year <- df.num.mosquito$Year + 1

model.data <- merge(model.data, df.num.mosquito, by = c("Year","orig.Trap","Species","Month"), all.x = TRUE)
model.data$NumMosquitos[which(model.data$NumMosquitos == -999)] = 
  model.data$mean.mos.count[which(model.data$NumMosquitos == -999)]
model.data$NumMosquitos[which(model.data$NumMosquitos == -999)] = 
  mean(model.data$NumMosquitos[which(model.data$Id == -999)])
model.data$mean.mos.count <- NULL
model.data$NumMosquitos[is.na(model.data$NumMosquitos)] <- mean(model.data$NumMosquitos[model.data$Id == -999])

# Block
model.data$Block <- as.factor(model.data$Block)

# Station Distances
station.1.coords = c(41.995, -87.933)
station.2.coords = c(41.786, -87.752)
model.data$Station.1 = distHaversine(as.matrix(cbind(model.data$Latitude, model.data$Longitude)), station.1.coords)
model.data$Station.2 = distHaversine(as.matrix(cbind(model.data$Latitude, model.data$Longitude)), station.2.coords)

# Refine weather variables with stations distance
weather.vars <- c("Tmax","Tmin","Tavg","DewPoint","WetBulb","Heat","Cool","PrecipTotal",
                  "StnPressure","SeaLevel","ResultSpeed","ResultDir","AvgSpeed")

wt.1 = model.data$Station.1 / (model.data$Station.1 + model.data$Station.2)
wt.2 = 1.0 - wt.1
for(var in weather.vars){
  var.1 = paste0(var,".1",sep="")
  var.2 = paste0(var,".2",sep="")
  model.data[[var]] = (wt.1 * model.data[[var.1]]) + 
    (wt.2 * model.data[[var.2]])
  model.data = model.data[,-which(names(model.data) %in% c(var.1,var.2))]
}


  
# Quick EDA check
df.stats.num <- statsNum(model.data)
df.stats.nnum <- statsNonNum(model.data)

##=========================================================================================
## Write model.data & score.data
##=========================================================================================
# score.data will be the test data set for submission
score.data = model.data[is.na(model.data$WnvPresent),]
saveRDS(score.data,"score.data.rds") 

# model.data
model.data = model.data[!is.na(model.data$WnvPresent),]
model.data$WnvPresent = factor(model.data$WnvPresent, levels = c("0","1"))
saveRDS(model.data,"model.data.rds") 
