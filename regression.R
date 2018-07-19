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
full_date_range  = data.frame(Date=seq(min(observation_matrix$Date),max(observation_matrix$Date), by='days'))

observation_matrix = observation_matrix %>%
  right_join(full_date_range) %>%
  arrange(Date)

######################
# add some different discharge statistics (Q)
Q <- read.csv('Discharge.csv')
Q$Date <- as.Date(Q$Date)
obs_Q <- Q[Q$Date %in% full_date_range$Date,]
samples[, c('Qmean', 'Qmax')] <- Q[match(samples$Date, Q$Date),c('Qmean', 'Qmax')]
samples.train[, c('Qmean', 'Qmax')] <- Q[match(samples.train$Date, Q$Date),c('Qmean', 'Qmax')]
samples.test[, c('Qmean', 'Qmax')] <- Q[match(samples.test$Date, Q$Date),c('Qmean', 'Qmax')]

#################
observation_matrix = observation_matrix %>%
  left_join(obs_Q, by='Date')

# Fill in some missing discharge values with the mean
observation_matrix$Qmean[c(365,376,377,391)] = mean(observation_matrix$Qmean, na.rm=TRUE)


###############
# Setup data for jags
observation_matrix$CHIL.count <- ceiling(observation_matrix$CHIL.count)
data = list(y=observation_matrix$CHIL.count,
            num_days=nrow(observation_matrix),
            discharge = observation_matrix$Qmean)
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
  x[d] <- beta[1] + beta[2]*x[d-1] + beta[3]*discharge[d-1]  ## process model
}

  
  #### Priors
  beta ~ dmnorm(b0,Vb)  	## multivariate Normal prior on vector of regression params
  x[1]  ~ dnorm(1,1)
  tau_add ~ dgamma(0.01,0.01)

}
"

## specify priors
data$b0 <- as.vector(c(0,0,0))      ## regression beta means
data$Vb <- solve(diag(10000,3))   ## regression beta precisions
################

# nchain = 3
# init <- list()
# for(i in 1:nchain){
#   y.samp = sample(y,length(y),replace=TRUE)
#   init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),tau_obs=5/var(log(y.samp)))
# }

#############
load.module('glm')
j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                        # inits = init,
                         n.chains = 3)

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c('mu','beta'),
                            n.iter = 5000)
#plot(jags.out)



time.rng = c(1,nrow(observation_matrix)) ## adjust to zoom in and out
out <- as.matrix(jags.out)
x.cols <- grep("mu",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975))

plot(full_date_range$Date,ci[2,],type='n', ylim=c(-10,60))

ecoforecastR::ciEnvelope(full_date_range$Date,ci[1,],ci[3,],col="lightBlue")
points(full_date_range$Date,log(observation_matrix$CHIL.count),pch="+",cex=0.5)

