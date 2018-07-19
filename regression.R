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
                        n.adapt = 3000)

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c('mu','beta'),
                            n.iter = 10000)

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

plot(full_date_range$Date,ci[2,],type='n', ylim=c(0,10))

ecoforecastR::ciEnvelope(full_date_range$Date,ci[1,],ci[3,],col="lightBlue")
points(full_date_range$Date,log(observation_matrix$CHIL.count),pch="+",cex=0.5)

































##'  Kalman Filter
##' @param  M   = model matrix
##' @param  mu0 = initial condition mean vector
##' @param  P0  = initial condition covariance matrix
##' @param  Q   = process error covariance matrix
##' @param  R   = observation error covariance matrix
##' @param  Y   = observation matrix (with missing values as NAs), time as col's
##'
##' @return list
##'  mu.f, mu.a  = state mean vector for (a)nalysis and (f)orecast steps
##'  P.f, P.a    = state covariance matrix for a and f
KalmanFilter <- function(M,mu0,P0,Q,R,Y){
  
  ## storage
  nstates = nrow(Y)  
  nt = ncol(Y)
  mu.f  = matrix(NA,nstates,nt+1)  ## forecast mean for time t
  mu.a  = matrix(NA,nstates,nt)  ## analysis mean for time t
  P.f  = array(NA,c(nstates,nstates,nt+1))  ## forecast variance for time t
  P.a  = array(NA,c(nstates,nstates,nt))  ## analysis variance for time t
  
  ## initialization
  mu.f[,1] = mu0
  P.f[,,1] = P0
  I = diag(1,nstates)
  
  ## run updates sequentially for each observation.
  for(t in 1:nt){
    
    ## Analysis step: combine previous forecast with observed data
    KA <- KalmanAnalysis(mu.f[,t],P.f[,,t],Y[,t],R,I)
    mu.a[,t] <- KA$mu.a
    P.a[,,t] <- KA$P.a
    
    ## Forecast step: predict to next step from current
    KF <- KalmanForecast(mu.a[,t],P.a[,,t],M,Q)
    mu.f[,t+1] <- KF$mu.f
    P.f[,,t+1] <- KF$P.f
  }
  
  return(list(mu.f=mu.f,mu.a=mu.a,P.f=P.f,P.a=P.a))
}

##' Kalman Filter: Analysis step
##' @param  mu.f = Forecast mean (vector)
##' @param  P.f  = Forecast covariance (matrix)
##' @param  Y    = observations, with missing values as NAs) (vector)
##' @param  R    = observation error covariance (matrix)
##' @param  H    = observation matrix (maps observations to states)
KalmanAnalysis <- function(mu.f,P.f,Y,R,H){
  obs = !is.na(Y) ## which Y's were observed?
  if(any(obs)){
    H <- H[obs,]                                              ## observation matrix
    K <- P.f %*% t(H) %*% solve(H%*%P.f%*%t(H) + R[obs,obs])  ## Kalman gain
    mu.a <- mu.f + K%*%(Y[obs] - H %*% mu.f)                  ## update mean
    P.a <- (1-K %*% H)*P.f                                    ## update covariance
  } else {
    ##if there's no data, the posterior is the prior
    mu.a = mu.f
    P.a = P.f
  }
  return(list(mu.a=mu.a,P.a=P.a))
}

##' Kalman Filter: Forecast Step
##' @param mu.a = analysis posterior mean (vector)
##' @param P.a  = analysis posterior covariance (matrix)
##' @param M    = model (matrix)
##' @param  Q   = process error covariance (matrix)
KalmanForecast <- function(mu.a,P.a,M,Q){
  mu.f = M%*%mu.a
  P.f  = Q + M*P.a*t(M)
  return(list(mu.f=mu.f,P.f=P.f))
}


