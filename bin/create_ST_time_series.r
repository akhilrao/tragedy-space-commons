#####################################################################
# Script to process Space-Track data JSONs into a time series of satellites, debris, launches, losses, and launch failures
#####################################################################

rm(list=ls())

### load packages
library(data.table) # Used for rbindlist()
#library(httr) # Used to pull from Space-Track API.
library(jsonlite) # Used to work with JSON data.
library(tidyr) # Used to reshape data from wide to long.
library(plyr)
library(dplyr) # Used to recode categorical variables from lookup tables.
library(stringr)

options(stringsAsFactors=FALSE)

setwd("~/Documents/git-repos/dissertation/data/Space-Track\ data")
historical_leo_satcat <- read_json("historical_leo_satcat.json",simplifyVector=TRUE)
setwd("~/Documents/git-repos/dissertation/data/JSR\ data")
JSR <- read.csv("launchlog.csv")
setwd("~/Documents/git-repos/dissertation/data/destructions\ data")
asat_destr <- read.csv("asat_destr.csv")
satsat_destr <- read.csv("satsat_destr.csv")
debsat_destr <- read.csv("debsat_destr.csv")
setwd("~/Documents/git-repos/tragedy-space-commons/data/")
fragrisk <- read.csv("fragrisk.csv")

############################################
### function definitions

up_now_count <- function(objects,year_vec) {
	count <- as.data.frame(matrix(0,nrow=length(year_vec),ncol=3))
	colnames(count) <- c("sat_count","launches","sat_decays")
	for(y in 1:length(year_vec)) {
		if(y==1) {
			up_this_year <- which(objects$LAUNCH_YEAR==year_vec[y])
			down_this_year <- which(objects$DECAY_YEAR==year_vec[y])
			netdiff <- length(up_this_year) - length(down_this_year)
			count[y,1] <- netdiff
		}
		if(y>1) {
			up_this_year <- which(objects$LAUNCH_YEAR==year_vec[y])
			down_this_year <- which(objects$DECAY_YEAR==year_vec[y])
			netdiff <- length(up_this_year) - length(down_this_year)
			#print(netdiff)
			count[y,1] <- count[y-1,1] + netdiff
			count[y-1,2] <- length(up_this_year)
			count[y-1,3] <- length(down_this_year)
		}
	}
	return(count)
}

launch_count <- function(objects,year_vec,JSR) {
	count <- rep(0,length=length(year_vec))
	for(y in 1:length(year_vec)) {
		if(JSR==1) {up_this_year <- intersect(which(objects$Suc=="S"),which(objects$LAUNCH_YEAR==year_vec[y])) }
		if(JSR==0) {up_this_year <- which(objects$LAUNCH_YEAR==year_vec[y])}
		count[y] <- length(up_this_year)
	}
	return(count)
}

decay_count <- function(objects,year_vec) {
	count <- rep(0,length=length(year_vec))
	for(y in 1:length(year_vec)) {
		down_this_year <- which(objects$DECAY_YEAR==year_vec[y])
		count[y] <- length(down_this_year)
	}
	return(count)
}

launchfail_count <- function(dfrm,year_vec) {
	count <- rep(0,length=length(year_vec))
	for(y in 1:length(year_vec)) {
		thisyear <- which(dfrm$LAUNCH_YEAR==year_vec[y])
		fails <- which(dfrm$Suc=="F")
		num_failed <- intersect(thisyear,fails)
		count[y] <- length(num_failed)
	}
	return(count)
}

############################################
###  create LEO object time series

year_vec <- seq(from=as.numeric(min(historical_leo_satcat$LAUNCH_YEAR[complete.cases(historical_leo_satcat$LAUNCH_YEAR)])), to=as.numeric(max(historical_leo_satcat$LAUNCH_YEAR[complete.cases(historical_leo_satcat$LAUNCH_YEAR)])), by=1)

## ST data cleaning
historical_leo_satcat$LAUNCH_YEAR <- as.numeric(historical_leo_satcat$LAUNCH_YEAR)
historical_leo_satcat$DECAY_YEAR <- as.numeric(format(as.Date(historical_leo_satcat$DECAY,"%Y-%m-%d"),"%Y"))

## JSR data cleaning
JSR$LAUNCH_YEAR <- substr(JSR$Launch.Date..UTC., 1, 4)
# impute left-out fields with whatever was going on just before (missing values are multiple-manifested launches)
for(i in 1:length(JSR$LAUNCH_YEAR) ) {
	if(JSR$LAUNCH_YEAR[i]=="") {JSR$LAUNCH_YEAR[i] <- JSR$LAUNCH_YEAR[i-1]}
}
for(i in 1:length(JSR$Suc) ) {
	if(JSR$Suc[i]=="") {JSR$Suc[i] <- JSR$Suc[i-1]}
}

## Processing
satellites <- historical_leo_satcat[which(historical_leo_satcat$OBJECT_TYPE=="PAYLOAD"),]
debris <- historical_leo_satcat[which(historical_leo_satcat$OBJECT_TYPE!="PAYLOAD"),]

sat_stock <- up_now_count(satellites,year_vec)

deb_count <- up_now_count(debris,year_vec)

#launches <- launch_count(satellites,year_vec,0)
#sat_decays <- decay_count(satellites,year_vec)
launch_fails <- launchfail_count(JSR,year_vec)

stock_series_1 <- data.frame(year=year_vec,
							payloads_in_orbit=sat_stock$sat_count,
							launch_successes=sat_stock$launches, 
							payloads_decayed=sat_stock$sat_decays,
							launch_failures=launch_fails,
							jspoc_debris=deb_count$sat_count)

stock_series <- merge(stock_series_1,asat_destr,by=c("year"),all.x=TRUE,all.y=TRUE)
stock_series <- merge(stock_series,satsat_destr,by=c("year"),all.x=TRUE,all.y=TRUE)
stock_series <- merge(stock_series,debsat_destr,by=c("year"),all.x=TRUE,all.y=TRUE)
stock_series <- merge(stock_series,fragrisk,by=c("year"),all.x=TRUE,all.y=TRUE)

stock_series[is.na(stock_series)] <- 0

write.csv(stock_series,file="ST_stock_series.csv",row.names=FALSE)
