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

# Check for data.table package, install if not found
if('data.table' %in% rownames(installed.packages()) == FALSE) {
	install.packages('data.table', repos = 'https://cran.cnr.berkeley.edu')
}

# Load requisite packages
require(foodbase)
require(dataRetrieval)
require(data.table)

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
	c('00010', '00060', '00065', '00095'), 
	startDate = min(d2$Samples$Date), endDate = max(d2$Samples$Date), 
	tz = 'America/Phoenix')
	# Function returns data on water temperature, discharge, stage, 
	# and specific conductance for the range of drift sampling dates.
	# Takes a couple minutes to download it all. Save locally thereafter.

# Get only environmental data columns of interest, rename
e1 <- e0[, c('dateTime', 'X_00010_00000', 'X_00060_00000', 'X_00065_00000',
	'X_00095_00000')]
	colnames(e1) <- c('DateTime', 'Temperature', 'Discharge', 'Stage', 'SpC')
	
# Write environmental data locally to save time for future runs 
dbdir <- paste0(find.package('foodbase'),'/Data')
	# Location of the foodbase package "Data" folder.
write.csv(e1, paste0(dbdir, '/GageData.csv'), row.names = FALSE)

# Read in these data (start from here after first code run)
e1 <- read.csv(paste0(dbdir, '/GageData.csv'))


# Clean up data ----------------------------------------------------------------

# Convert DateTimes to POSIX
e2 <- e1
	e2$DateTime <- strptime(e1$DateTime, format = '%Y-%m-%d %H:%M:%S', tz = 'MST')

# Merge sample and environmental data
d3 <- d2
e2T <- e2[!is.na(e2$Temperature),]
e2S <- e2[!is.na(e2$SpC),]
drifttime0 <- strptime(paste(d3$Samples$Date, d3$Samples$TimeBegin), 
	format = '%Y-%m-%d %H:%M:%S', tz = 'MST')
d3$Samples <- d2$Samples[order(drifttime0),]
	d3$Samples$SpC <- d3$Samples$Stage <- d3$Samples$Discharge <- 
	d3$Samples$Temperature <- NA
drifttime <- suppressWarnings(data.table(Drift = drifttime0))
	drifttime[, merge:= Drift]
envtime <- suppressWarnings(data.table(Env = e2$DateTime))
	envtime[, merge:= Env]
envtimeT <- suppressWarnings(data.table(Env = e2T$DateTime))
	envtimeT[, merge:= Env]
envtimeS <- suppressWarnings(data.table(Env = e2S$DateTime))
	envtimeS[, merge:= Env]
setkeyv(drifttime, c('merge'))
setkeyv(envtime, c('merge'))
setkeyv(envtimeT, c('merge'))
setkeyv(envtimeS, c('merge'))
closest <- envtime[drifttime, roll = 'nearest']
closestT <- envtimeT[drifttime, roll = 'nearest']
closestS <- envtimeS[drifttime, roll = 'nearest']
d3$Samples[, c('Discharge', 'Stage')] <- 
	e2[match(closest$Env, as.POSIXct(e2$DateTime)), c('Discharge', 'Stage')]
d3$Samples[, 'Temperature'] <- 
	e2T[match(closestT$Env, as.POSIXct(e2T$DateTime)), 'Temperature']
d3$Samples[, 'SpC'] <- 
	e2S[match(closestS$Env, as.POSIXct(e2S$DateTime)), 'SpC']