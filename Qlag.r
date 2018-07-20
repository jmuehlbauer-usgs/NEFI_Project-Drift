###############
# Compare midge counts to previous day's discharge (1 day lag)
# Paste in lieu of similar code in other files.

###############
# Setup data for jags
observation_matrix$CHIL.count <- ceiling(observation_matrix$CHIL.count)
data = list(y=observation_matrix$CHIL.count[-1],
            num_days=nrow(observation_matrix)-1,
            discharge = observation_matrix$Qmean[-dim(observation_matrix)[1]])
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
                            variable.names = c('mu','beta'),
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

plot(full_date_range$Date[-1],ci[2,],type='n', ylim=c(0,10), ylab = 'Log Midge Count', xlab = 'Date', main = 'Regression - Discharge Lag, mean')

ecoforecastR::ciEnvelope(full_date_range$Date[-1],ci[1,],ci[3,],col="lightBlue")
points(full_date_range$Date,log(observation_matrix$CHIL.count),pch="+",cex=0.5)
