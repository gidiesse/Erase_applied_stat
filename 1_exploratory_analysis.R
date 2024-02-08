# Exploratory Analysis: this is the preliminary investigation into the Metadata. 
# We were particularly interested in figuring out where our NAs lay and if there
# was any way to interpolate them. We noticed that the features date and M
# acro_Period were linked and also that the features Macroarea and lat_num,
# lat_long were linked so we thought that for individuals where only one of the 
# features in the pair was missing we could use the other to interpolate it.


### SET WD AND LOAD USED PACKAGES ###
setwd("~/Desktop/ERASE_project2023_desgrp")
library("Dict")

### LOAD DATASET ###
metadata <- read.csv("~/Desktop/ERASE_project2023_desgrp/data/erase_metadata.csv", header=TRUE)
head(metadata)
dim(metadata)

### EDA ON NAs ###
# total number of NAs in dataset
sum(is.na(metadata))

# how many in each column?
count_na <- function(x) sum(is.na(x))
missing_num <- sapply(metadata, count_na)
missing <- rbind(missing_num)
missing

# we remove the three columns that are completely empty, no need to analyse them 
# now since Dr.Raveane will probably fill them in later (ask him!)
sub_md = subset(metadata, select = -c(5,6,14) )

# use Dict package to get a data structure which has the index of NAs for each column
# these indexes will be used later to identify the entries to be interpolated
na_index = dict(keyvalues=colnames(sub_md))

for (col in colnames(sub_md)) {
  na_index[col] = which(is.na(sub_md[[col]]), arr.ind = TRUE)
  }

# this dict contains the indexes of the rows with NAs for each column
na_index

# we notice that Macro_Period and date are two linked features and that if one is 
# missing, it can be interpolated using the other 

# PROBLEM: if both are missing, we can't interpolate in this way.
# We want to see how many rows are missing both columns

# variable that contains indexes of rows with missing values on both columns
na_date_period <- intersect(na_index["Macro_Period"],na_index["date"])

# all the rows that have NAs in Macro_Period also have NAs in date which makes sense
# because if the Macro_Period is unknown, it is impossible to have the exact date 
na_date_period_len <- length(intersect(na_index["Macro_Period"],na_index["date"]))
na_period_len <- length(na_index["Macro_Period"])
na_period_len == na_date_period_len

# for 158 observations, we have the Macro_Period but no date, can we interpolate 
# the date?
# ie for an observation missing the date feature but not the Macro_Period feature
# we take the mode of the dates associated with that Macro-Period and use that to 
# fill the missing date
# eg: MiddleAges <-> [mode_date_MA]
# Obs0$Macro_Period == "MiddleAges" --> Obs0$date = mode_date_MA

unique(sub_md$Macro_Period)
no_interp_candidate_len = length(na_index["date"])-length(na_index["Macro_Period"])

# extract indexes of the candidate rows to be interpolated (no date, yes Macro_Period!)
check_presence_dp <- function(x) x %in% na_date_period
na_interp_indexes <- subset(na_index["date"], !check_presence_dp(na_index["date"]))

# we can apply a similar logic to the one above to the columns lat, long and Macroarea

# since the intersection between the columns lat and long has the same number of 
# NAs as lat, we conclude that lat and long always have NAs in the same rows
length(intersect(na_index["lat"], na_index["long"])) == length(na_index["lat"])

# for these 440 observations, we are missing both Macroarea and lat, long values
# so the interpolation will have to be done in another way 
na_area_coord <- intersect(na_index["lat"], na_index["Macroarea"])
length(intersect(na_index["lat"], na_index["Macroarea"]))

# extract indexes of the candidate rows to be interpolated (no Macroarea, yes coord)
# lat and long!)
check_presence_ac <- function(x) x %in% na_area_coord

na_interp_area <- subset(na_index["Macroarea"], !check_presence_ac(na_index["Macroarea"]))
length(na_interp_area)

na_interp_coord <- subset(na_index["lat"], !check_presence_ac(na_index["lat"]))
length(na_interp_coord)

# we want to find the names of all the possible Macroareas
M_area = unique(sub_md$Macroarea)

# this is unfinished since there is a problem with the formatting of the lat and long, 
# we will finish when this has been fixed
est_eur_index <- c(which(sub_md$Macroarea == "Eastern_Europe"))

# idea for interpolation of missing Macroarea from coordinates: 
# find the max and min of lat and long amongst the data for each Macroarea and 
# build a grid mapping out the macroareas based on this information.
# eg: take Eastern Europe <-> {[min_lat, max_lat], [min_long, max_long]} 
# obs0 = [lat0, long0] belongs to Eastern Europe bounds -> interpolate obs0
# Macroarea as Eastern Europe.

# idea for interpolation of missing coordinates from Macroarea: 
# (1): for each Macroarea, we will calculate the mode of the long and the lat and use
# these values to fill in the NAs 
# eg: Eastern Europe <-> [mode_lat, mode_long] 
# if obs0$Macroarea == "Eastern Europe" ==> obs0$lat = mode_lat, obs0$long = mode_long

# (2): for each Macroarea, decide on a "capital" and use the coordinates of that place 
# to fill in the NAs
# eg: Eastern Europe <-> [capital_lat, capital_long]
# if obs0$Macroarea == "Eastern Europe" ==> obs0$lat = capital_lat, 
# obs0$long = capital_long