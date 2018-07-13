# Analysis of Lees Ferry drift data, forecasting for bug flows -----------------
# Class project for Near-term Ecological Forecasting Course
# Boston, MA, 2018
# Document created by J. D. Muehlbauer, jmuehlbauer@usgs.gov, 13 July 2018


# Load packages ----------------------------------------------------------------

# Check for foodbase package, install if not found
if('foodbase' %in% rownames(installed.packages()) == FALSE) {
	if('devtools' %in% rownames(installed.packages()) == FALSE) {
		install.packages('devtools', repos = 'https://cran.cnr.berkeley.edu')
	}
	require(devtools)
	install_github(repo='jmuehlbauer-usgs/R-packages',subdir='foodbase')
}

# Check for dataRetrieval package, install if not found
if('dataRetrieval' %in% rownames(installed.packages()) == FALSE) {
	install.packages('dataRetrieval', repos = 'https://cran.cnr.berkeley.edu')
}

# Load requisite packages
require(foodbase)
require(dataRetrieval)


# Get data ---------------------------------------------------------------------

# Drift data (from Lees Ferry, river mile 0)
d0 <- readDB(type = 'Sample', gear = 'Drift')
	# Function queries drift data pushed from our lab's database to GitHub.
	# Returns drift sample info for all samples in database.
d1 <- d0[d0$RiverMile < 0.5 & d0$RiverMile > -0.5 & d0$Region == 'GrandCanyon',]
	# Reduce samples dataset to just sample collected in Lees Ferry
	# (River mile 0 +/- half a mile).
d2 <- sampspec(samp = d1, gear = 'Drift', stats = TRUE)
	# Function returns all sample info and specimen counts data for the samples. 
	# Salient list elements are the tables d1$Samples and d1$Specimens.
	
# Environmental data from Lees Ferry gage
e0 <- readNWISuv(siteNumbers = '09380000', parameterCd = 
	c('00010', '00045', '00060', '00065', '00095', '00300', '00301'), 
	startDate = min(d2$Samples$Date), endDate = max(d2$Samples$Date))
	# Function returns data on water temperature, precipitation, discharge, stage, 
	# specific conductance, dissolved oxygen, and DO percent saturation
	# for the range of drift sampling dates.
	# Takes a couple minutes to download it all. Save locally thereafter.

# Write environmental data locally to save time for future runs 
dbdir <- paste0(find.package('foodbase'),'/Data')
	# Location of the foodbase package "Data" folder.
write.csv(e0, paste0(dbdir, '/GageData.csv'))

# Read in these data (start from here after first code run)
e0 <- read.csv(paste0(dbdir, '/GageData.csv'))


	