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
  w[d] <- x[d-1] + beta*discharge[d-1]  ## process model
  x[d]~dnorm(w[d], tau_add)
}

  
  #### Priors
  beta ~ dmnorm(b0,Vb)  	## multivariate Normal prior on vector of regression params
  x[1]  ~ dnorm(1,1)
  tau_add ~ dgamma(0.01,0.01)

}
"

## specify priors
data$b0 <- as.vector(c(0))      ## regression beta means
data$Vb <- solve(diag(10000,1))   ## regression beta precisions
################

# inits <- list()
# inits[[1]] <- list(beta = c(5,0,0))
# inits[[2]] <- list(beta = c(5,0,0))
# inits[[3]] <- list(beta = c(5,0,0))


#############
load.module('glm')
j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                        # inits = inits,
                         n.chains = 3,
                        n.adapt = 1000)

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c('mu','beta','tau_add'),
                            n.iter = 5000)

# out <- list(params = NULL, predict = NULL)
# mfit <- as.matrix(jags.out, chains = TRUE)
# pred.cols <- grep("mu[", colnames(mfit), fixed = TRUE)
# chain.col <- which(colnames(mfit) == "CHAIN")
# out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
# out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
# plot(out$params)
# plot(out$predict)


time.rng = c(1,nrow(observation_matrix)) ## adjust to zoom in and out
out <- as.matrix(jags.out)
x.cols <- grep("mu",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(log(out[,x.cols]),2,quantile,c(0.025,0.5,0.975))

# Plot the full time range. This is *with* dishcarge
plot(full_date_range$Date,ci[2,],type='n',ylim=c(-5, 10), xlim = c(as.Date('2008-01-01'),as.Date('2015-12-31')),
     xlab = 'Date (Jan. 2008 - Dec. 2015)', ylab='log(Midge Count)') 
ecoforecastR::ciEnvelope(full_date_range$Date,ci[1,],ci[3,],col="lightBlue") 
points(full_date_range$Date,log(observation_matrix$CHIL.count),pch="+",cex=1.5) 
################################
# Plot zoomed into the forecast period

plot(full_date_range$Date,ci[2,],type='n',ylim=c(-5, 10), xlim = c(as.Date('2015-07-01'),as.Date('2016-07-01')),
     ylab = 'log(Midge Count)',xlab='Date (July 2015 - July 2016)')

ecoforecastR::ciEnvelope(full_date_range$Date,ci[1,],ci[3,],col="lightBlue")
points(full_date_range$Date,log(observation_matrix$CHIL.count),pch="+",cex=2.5)
points(samples.test$Date, log(samples.test$CHIL.count), cex=2.5)




######################################

forecast_discharge = function(IC, beta, tau_add, nT, nMC, disharge){
  state_storage = matrix(NA, nrow=nMC, ncol = nT)
  state_storage[,1]=IC
  for(t in 2:nT){
    temp_storage = state_storage[,t-1] + beta*discharge[t-1]
    state_storage[,t] = temp_storage + rnorm(nMC,mean=0,sd=1/sqrt(tau_add))
    #state_storage[,t] = temp_storage
  }
  return(state_storage)
}


nMC = 4000
nT=150

#########################################
# Forecast forward using only initial condition variance
mc_rows = sample.int(nrow(out), size=nMC)
# Last time step of observed data is 2015-12-20, id 2989 
IC = log(out[mc_rows,'mu[2989]'])

beta = mean(out[mc_rows,'beta'])
tau_add = mean(out[mc_rows,'tau_add'])
discharge = observation_matrix$Qmean[2989:(2989+nT)]

IC_only_forecast = forecast_discharge(IC=IC, beta=beta,tau_add=tau_add, nT=nT, nMC=nMC, disharge = disharge)

IC_only_ci <- apply(IC_only_forecast,2,quantile,c(0.025,0.5,0.975),na.rm=T)

#######################################
# forecast forward using initial condition variance AND variance and the beta parameter for dishcarge
beta = out[mc_rows,'beta']

IC_beta_forecast = forecast_discharge(IC=IC, beta=beta, tau_add=tau_add, nT=nT, nMC=nMC, disharge = disharge)
IC_beta_ci <- apply(IC_beta_forecast,2,quantile,c(0.025,0.5,0.975),na.rm=T)


################################
# Plot zoomed into the forecast period

plot(1:nT,IC_only_ci[2,],type='n',ylim=c(-5, 10), xlim = c(1,nT),
     ylab = 'log(Midge Count)',xlab='Timesteps of forecast (Days)')

ecoforecastR::ciEnvelope(1:nT,IC_beta_ci[1,],IC_beta_ci[3,],col="Red")
ecoforecastR::ciEnvelope(1:nT,IC_only_ci[1,],IC_only_ci[3,],col="grey")






