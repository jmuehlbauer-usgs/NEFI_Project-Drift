library(tidyverse)

# read in CSV files 
samples <- read.csv("Samples.csv"); head(samples); dim(samples)
specimens <- read.csv("Specimens.csv"); head(specimens)
stats <- read.csv("Statistics.csv"); head(stats)

samples$CHIL.count <- stats[stats$SpeciesID == "CHIL", "CountTotal"]
samples$Date <- as.Date(samples$Date, "%m/%d/%Y")
samples.test <- samples[samples$Date > as.Date("2016-01-01"),]; dim(samples.test)
samples.train <- samples[samples$Date < as.Date("2016-01-01"),]; dim(samples.train)

# This code gives the mean value for every
# day which has multiple samples. 
observation_matrix = samples.train %>%
  group_by(Date) %>%
  summarize(CHIL.count = mean(CHIL.count)) %>%
  ungroup()
##########################################

# Now fill in unsampled dates with NA
full_date_range = data.frame(Date=seq(min(observation_matrix$Date),max(observation_matrix$Date), by='days'))

# add some different discharge statistics (Q)
Q <- read.csv('Discharge.csv')
Q$Date <- as.Date(Q$Date)
obs_Q <- Q[Q$Date %in% full_date_range$Date,]
samples[, c('Qmean', 'Qmax')] <- Q[match(samples$Date, Q$Date),c('Qmean', 'Qmax')]
samples.train[, c('Qmean', 'Qmax')] <- Q[match(samples.train$Date, Q$Date),c('Qmean', 'Qmax')]
samples.test[, c('Qmean', 'Qmax')] <- Q[match(samples.test$Date, Q$Date),c('Qmean', 'Qmax')]
