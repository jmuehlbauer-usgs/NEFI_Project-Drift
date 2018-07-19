# script for analyzing grand canyon midge data
library(tidyverse)
library(rjags)
library(ecoforecastR)

# read in CSV files 
samples <- read.csv("Samples.csv"); head(samples); dim(samples)
specimens <- read.csv("Specimens.csv"); head(specimens)
stats <- read.csv("Statistics.csv"); head(stats)

samples$CHIL.count <- stats[stats$SpeciesID == "CHIL", "CountTotal"]
samples$Date <- as.Date(samples$Date, "%m/%d/%Y")
samples.test <- samples[samples$Date > as.Date("2016-01-01"),]; dim(samples.test)
samples.train <- samples[samples$Date < as.Date("2016-01-01"),]; dim(samples.train)

###########################################

# date x sample observations matrix
# where there can be >= 1 observation a day
# This code is very CPU/memory intensive so lets save it for later
# observation_matrix = samples.train %>%
#   group_by(Date) %>%
#   mutate(daily_obs = 1:n()) %>%
#   ungroup() %>%
#   select(Date, CHIL.count,daily_obs) %>%
#   spread(daily_obs, CHIL.count)

# This code gives the mean value for every
# day which has multiple samples. 
observation_matrix = samples.train %>%
  group_by(Date) %>%
  summarize(CHIL.count = mean(CHIL.count)) %>%
  ungroup()
##########################################

# Now fill in unsampled dates with NA
full_date_range = data.frame(Date=seq(min(observation_matrix$Date),as.Date('2016-07-01'), by='days'))

observation_matrix = observation_matrix %>%
  right_join(full_date_range) %>%
  arrange(Date) %>%
  select(-Date)


###############
# Setup data for jags
observation_matrix$CHIL.count <- ceiling(observation_matrix$CHIL.count)
data = list(y=observation_matrix$CHIL.count,
            num_days=nrow(observation_matrix))
###############
# This
RandomWalk = "
model{
  
  #### Data Model
  for(d in 1:num_days){
      mu[d] <- exp(x[d])
      y[d] ~ dpois(mu[d])
  }
  
  #### Process Model
  for(d in 2:num_days){
      x[d]~dnorm(x[d-1], tau_add)
  }
  
  #### Priors
  x[1]  ~ dnorm(1,1)
 # tau_obs ~ dgamma(0.01,0.01)
  tau_add ~ dgamma(0.01,0.01)
}
"

################

# nchain = 3
# init <- list()
# for(i in 1:nchain){
#   y.samp = sample(y,length(y),replace=TRUE)
#   init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),tau_obs=5/var(log(y.samp)))
# }

#############

j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                        # inits = init,
                         n.chains = 3)

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c('mu','tau_add'),
                            n.iter = 5000)
#plot(jags.out)



time.rng = c(1,nrow(observation_matrix)) ## adjust to zoom in and out
out <- as.matrix(jags.out)
x.cols <- grep("^mu",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975))

plot(full_date_range$Date,ci[2,],type='n',ylim=c(-100, 2000), xlim = c(as.Date('2015-07-01'),as.Date('2016-07-01')),
     ylab = 'Midge Count',xlab='Date (July 2015 - July 2016')

ecoforecastR::ciEnvelope(full_date_range$Date,ci[1,],ci[3,],col="lightBlue")
points(full_date_range$Date,observation_matrix$CHIL.count,pch="+",cex=2.5)
points(samples.test$Date, samples.test$CHIL.count, cex=2.5)
