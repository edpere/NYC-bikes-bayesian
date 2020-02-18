###
### Before this script run utils.R
### All the RData containg the datasets can be build with data_preproc.R script
###


library(spBayes)
library(MBA)
library(fields)
library(geoR)
library(sp)
library(maptools)
library(rgdal)
library(MASS)
library(RColorBrewer)
library(gstat)
library(sf)


## Load training set ##
# Make sure that the current working director is correct

# Change TRUE / FALSE to run both / one training set
train_both = TRUE

load("Train_Sud_145stations.RData")
if(train_both) {
  tmp = Station.data
  load("Train_Nord_106stations.RData")
  Station.data = rbind(tmp, Station.data)
}

## NOTE ON TUNING PARAMETERS
# for smaller dataset (train_both = FALSE) use 3 for phi, 0.15 for sigma.sq, 0.2 for tau.sq
# for large dataset (train_both = TRUE) use 5 for phi, 0.05 for sigma.sq, 0.05 for tau.sq

rm(train_both, tmp)
rm(lat_max, lat_min, long_max, long_min) # for the moment I remove them (what if we use both?)


coords = as.matrix(Station.data[,c("Longitude","Latitude")])

DEMAND = Station.data[,c("N.Trips")]
beta.POPULATION = Station.data[,c("Block.population")]
beta.LANES = Station.data[,c("Lane.count")]
beta.SUBWAY = Station.data[,c("Dist.metro")]
beta.PROXIMITY = Station.data[,c("Proximity.score")]
beta.LANDMARKS = Station.data[,c("Landmarks")]

##### spLM SETUP #####
n.samples = 20000 #try also 5000,10000,20000

## Priors specification
W = st_distance(st_as_sf(Station.data, coords = c("Longitude", "Latitude"), crs = 4326))
attr(W,"units") = NULL
attr(W,"class") = NULL
diag(W) = min(W[W!=0])
# the closer is the minimum distance in the dataset, the lower the acceptance rate
min_dist = min(W)
max_dist = max(W)
### Priors for phi ###
phi.prior.a = -log(0.05)/max_dist
phi.prior.b = -log(0.05)/min_dist



# Linear model frequentist fit
freq_model = lm(DEMAND~beta.POPULATION+beta.LANES+beta.SUBWAY)
summary(freq_model)

### Priors for beta,sigma2,tau2 ###
beta.ini <- as.numeric(freq_model$coeff)
(summary(freq_model)$sigma)^2 # estimated variance of residuals
# Rate : sets the proportion of spatial and random process in the priors for sigma2 and tau2
rate = 0.8
phi.ini = 0.035 # arbitrary value between prior.a and prior.b
sigma2.ini = 1/rgamma(1,2,(summary(freq_model)$sigma)^2*rate) # beta = rate% of res std error
tau2.ini = 1/rgamma(1,2,(summary(freq_model)$sigma)^2*(1-rate)) # beta = (1-rate)% of res std error


rm(W, min_dist, max_dist)

###############################
###### 1. STANDARD MODEL ######
###############################
sp.model.v1 <- spLM(DEMAND~beta.POPULATION+beta.LANES+beta.SUBWAY,
                    data=Station.data, coords=coords,
                    starting=list("phi"=phi.ini,"sigma.sq"=sigma2.ini,
                                  "tau.sq"=tau2.ini),
                    tuning=list("phi"=5, "sigma.sq"=0.05, #try small tuning params
                                "tau.sq"=0.05),
                    priors=list("beta.flat", "phi.Unif"=c(phi.prior.a, phi.prior.b),
                                "sigma.sq.IG"=c(2, (summary(freq_model)$sigma)^2*rate),
                                "tau.sq.IG"=c(2, (summary(freq_model)$sigma)^2*(1-rate))),
                    cov.model="exponential", # EXPONENTIAL COVARIANCE function
                    n.samples=n.samples)

# Using smaller tuning parameters the acceptance rate can easily increase
# Facendo tendere la distanza minima per phi prior a 0, l'acceptance rate diminuisce
# meno di quanto non faccia aumentando i tuning parameters

#summary(mcmc(sp.model$p.theta.samples))
# Posterior samples of beta coefficients and spatial effects (w)
burn.in = floor(0.5*n.samples) # try 0.25,0.5,0.75
sp.model.v1 = spRecover(sp.model.v1, start=burn.in)
beta.samples = sp.model.v1$p.beta.recover.samples
w.samples = sp.model.v1$p.w.recover.samples
#summary(beta.samples)
#summary(w.samples)

# Traceplots and posterior marginal distribution of beta parameters
x11()
par(mai=rep(0.4,4))
plot(beta.samples[,1:4])
# Traceplots and posterior marginal distribution of covariance parameters
x11()
par(mai=rep(0.4,4))
plot(sp.model.v1$p.theta.samples[,1:3])

sp.model.v1.mc = mcmc(sp.model.v1$p.theta.samples)
# Acceptance rate
1-rejectionRate(sp.model.v1.mc)
# Autocorrelation plot
x11()
acfplot(sp.model.v1.mc, lag.max=100)
# Cumulative mean plot
x11()
cumuplot(sp.model.v1.mc)
# Effective sample size:
effectiveSize(sp.model.v1.mc)


# Goodness of fit
lpml.v1 = LPML_fun(sp.model.v1)
waic.v1 = WAIC(sp.model.v1)


################################
####### 2.with PROXIMITY #######
################################
n.samples = 20000 #try also 5000,10000,20000
freq_model = lm(DEMAND~beta.POPULATION+beta.LANES+beta.SUBWAY+beta.PROXIMITY)
(summary(freq_model)$sigma)^2 # estimated variance of residuals
sigma2.ini = 1/rgamma(1,2,(summary(freq_model)$sigma)^2*rate) # beta = rate% of res std error
tau2.ini = 1/rgamma(1,2,(summary(freq_model)$sigma)^2*(1-rate)) # beta = (1-rate)% of res std error


sp.model.v2 <- spLM(DEMAND~beta.POPULATION+beta.LANES+beta.SUBWAY+beta.PROXIMITY,
                    data=Station.data, coords=coords,
                    starting=list("phi"=phi.ini,"sigma.sq"=sigma2.ini,
                                  "tau.sq"=tau2.ini),
                    tuning=list("phi"=5, "sigma.sq"=0.05,
                                "tau.sq"=0.05),
                    priors=list("beta.flat", "phi.Unif"=c(phi.prior.a, phi.prior.b),
                                "sigma.sq.IG"=c(2, (summary(freq_model)$sigma)^2*rate),
                                "tau.sq.IG"=c(2, (summary(freq_model)$sigma)^2*(1-rate))),
                    cov.model="exponential", # EXPONENTIAL COVARIANCE function
                    n.samples=n.samples)

#summary(mcmc(sp.model.v2$p.theta.samples))
# Posterior samples of beta coefficients and spatial effects (w)
burn.in = floor(0.5*n.samples) #0.25,0.5,0.75
sp.model.v2 = spRecover(sp.model.v2, start=burn.in)
beta.samples.v2 = sp.model.v2$p.beta.recover.samples
w.samples.v2 = sp.model.v2$p.w.recover.samples

# Traceplots and posterior marginal distribution of beta parameters
x11()
par(mai=rep(0.4,4))
plot(beta.samples.v2[,1:5])
# Traceplots and posterior marginal distribution of covariance parameters
x11()
par(mai=rep(0.4,4))
plot(sp.model.v2$p.theta.samples[,1:3])

sp.model.v2.mc = mcmc(sp.model.v2$p.theta.samples)
# Acceptance rate
1-rejectionRate(sp.model.v2.mc)
# Autocorrelation plot
x11()
acfplot(sp.model.v2.mc,lag.max = 100)
# Cumulative mean plot
x11()
cumuplot(sp.model.v2.mc)
# Effective sample size:
effectiveSize(sp.model.v2.mc)
# Goodness of fit
lpml.v2 = LPML_fun(sp.model.v2)
waic.v2 = WAIC(sp.model.v2)


##############################################
####### 3.with PROXIMITY and LANDMARKS #######
##############################################
n.samples = 20000 #try also 5000,10000,20000
freq_model = lm(DEMAND~beta.POPULATION+beta.LANES+beta.SUBWAY+beta.PROXIMITY+beta.LANDMARKS)
(summary(freq_model)$sigma)^2 # estimated variance of residuals
sigma2.ini = 1/rgamma(1,2,(summary(freq_model)$sigma)^2*rate) # beta = rate% of res std error
tau2.ini = 1/rgamma(1,2,(summary(freq_model)$sigma)^2*(1-rate)) # beta = (1-rate)% of res std error

sp.model.v3 <- spLM(DEMAND~beta.POPULATION+beta.LANES+beta.SUBWAY+beta.PROXIMITY+beta.LANDMARKS,
                    data=Station.data, coords=coords,
                    starting=list("phi"=phi.ini,"sigma.sq"=sigma2.ini,
                                  "tau.sq"=tau2.ini),
                    tuning=list("phi"=5, "sigma.sq"=0.05, 
                                "tau.sq"=0.05),
                    priors=list("beta.flat", "phi.Unif"=c(phi.prior.a, phi.prior.b),
                                "sigma.sq.IG"=c(2, (summary(freq_model)$sigma)^2*rate),
                                "tau.sq.IG"=c(2, (summary(freq_model)$sigma)^2*(1-rate))),
                    cov.model="exponential", # EXPONENTIAL COVARIANCE function
                    n.samples=n.samples)

#summary(mcmc(sp.model.v3$p.theta.samples))
# Posterior samples of beta coefficients and spatial effects (w)
burn.in = floor(0.5*n.samples) #0.25,0.5,0.75
sp.model.v3 = spRecover(sp.model.v3, start=burn.in)
beta.samples.v3 = sp.model.v3$p.beta.recover.samples
w.samples.v3 = sp.model.v3$p.w.recover.samples
#summary(beta.samples.v3)
#summary(w.samples.v3)

# Traceplots and posterior marginal distribution of beta parameters
x11()
par(mai=rep(0.4,4))
plot(beta.samples.v3[,1:6])
# Traceplots and posterior marginal distribution of covariance parameters
x11()
par(mai=rep(0.4,4))
plot(sp.model.v3$p.theta.samples[,1:3])

sp.model.v3.mc = mcmc(sp.model.v3$p.theta.samples)
# Acceptance rate
1-rejectionRate(sp.model.v3.mc)
# Autocorrelation plot
x11()
acfplot(sp.model.v3.mc, lag.max = 100)
# Cumulative mean plot
x11()
cumuplot(sp.model.v3.mc)
# Effective sample size:
effectiveSize(sp.model.v3.mc)

# Goodness of fit
lpml.v3 = LPML_fun(sp.model.v3)
waic.v3 = WAIC(sp.model.v3)


#################################
####### 4. with LANDMARKS #######
#################################
n.samples = 20000 #try also 5000,10000,20000
freq_model = lm(DEMAND~beta.POPULATION+beta.LANES+beta.SUBWAY+beta.LANDMARKS)
(summary(freq_model)$sigma)^2 # estimated variance of residuals
sigma2.ini = 1/rgamma(1,2,(summary(freq_model)$sigma)^2*rate) # beta = rate% of res std error
tau2.ini = 1/rgamma(1,2,(summary(freq_model)$sigma)^2*(1-rate)) # beta = (1-rate)% of res std error

sp.model.v4 <- spLM(DEMAND~beta.POPULATION+beta.LANES+beta.SUBWAY+beta.LANDMARKS,
                    data=Station.data, coords=coords,
                    starting=list("phi"=phi.ini,"sigma.sq"=sigma2.ini,
                                  "tau.sq"=tau2.ini),
                    tuning=list("phi"=5, "sigma.sq"=0.05, #try small tuning params
                                "tau.sq"=0.05),
                    priors=list("beta.flat", "phi.Unif"=c(phi.prior.a, phi.prior.b),
                                "sigma.sq.IG"=c(2, (summary(freq_model)$sigma)^2*rate),
                                "tau.sq.IG"=c(2, (summary(freq_model)$sigma)^2*(1-rate))),
                    cov.model="exponential", # EXPONENTIAL COVARIANCE function
                    n.samples=n.samples)

#summary(mcmc(sp.model.v4$p.theta.samples))
# Posterior samples of beta coefficients and spatial effects (w)
burn.in = floor(0.5*n.samples) #0.25,0.5,0.75
sp.model.v4 = spRecover(sp.model.v4, start=burn.in)
beta.samples.v4 = sp.model.v4$p.beta.recover.samples
w.samples.v4 = sp.model.v4$p.w.recover.samples
#summary(beta.samples.v4)
#summary(w.samples.v4)

# Traceplots and posterior marginal distribution of beta parameters
x11()
par(mai=rep(0.4,4))
plot(beta.samples.v4[,1:5])
# Traceplots and posterior marginal distribution of covariance parameters
x11()
par(mai=rep(0.4,4))
plot(sp.model.v4$p.theta.samples[,1:3])

sp.model.v4.mc = mcmc(sp.model.v4$p.theta.samples)
# Acceptance rate
1-rejectionRate(sp.model.v4.mc)
# Autocorrelation plot
x11()
acfplot(sp.model.v4.mc, lag.max = 100)
# Cumulative mean plot
#x11()
#cumuplot(sp.model.v4.mc)
# Effective sample size:
effectiveSize(sp.model.v4.mc)

# Goodness of fit
lpml.v4 = LPML_fun(sp.model.v4)
waic.v4 = WAIC(sp.model.v4)


##### SUMMARIZE GOODNESS OF FIT CRITERIA #####
gof = matrix(nrow=3, ncol=4)
rownames(gof) <- c("LPML", "WAIC", "MSE")
colnames(gof) <- c("model1", "model2", "model3", "model4")
gof[1,1] = lpml.v1
gof[1,2] = lpml.v2
gof[1,3] = lpml.v3
gof[1,4] = lpml.v4
gof[2,1] = waic.v1
gof[2,2] = waic.v2
gof[2,3] = waic.v3
gof[2,4] = waic.v4
gof

rm(lpml.v1, lpml.v2, lpml.v3, lpml.v4, waic.v1, waic.v2, waic.v3, waic.v4)
rm(LPML_fun, WAIC)
rm(n.samples, phi.ini, phi.prior.a, phi.prior.b, rate, sigma2.ini, tau2.ini, beta.ini)
rm(sp.model.v1.mc, sp.model.v2.mc, sp.model.v3.mc, sp.model.v4.mc)
rm(beta.LANDMARKS, beta.LANES, beta.POPULATION, beta.PROXIMITY, beta.SUBWAY)
rm(beta.samples, beta.samples.v2, beta.samples.v3, beta.samples.v4)
rm(burn.in, freq_model, w.samples, w.samples.v2, w.samples.v3, w.samples.v4)


######################
##### PREDICTION #####
######################

# Plot the prediction surface
# sp.model: model to use for prediction
# coords: coords of the prediction points
# covars: design matrix (i.e. covariates) of prediction points
# n.model: number of the model (used only for plot labelling)
predict <- function(sp.model, coords, covars, n.model) {
  pred = spPredict(sp.model, pred.coords = coords, pred.covars = covars)

  y.hat <- rowMeans(pred$p.y.predictive.samples)
  x11()
  y.pred.surf <- mba.surf(cbind(coords, y.hat), no.X=100, no.Y=100, extend=TRUE)$xyz.est
  image(y.pred.surf, xaxs = "r", yaxs = "r", main=paste("Predicted response Model", n.model))
  points(coords, pch=1, cex=1)
  contour(y.pred.surf, add=T)
  legend(1.5,2.5, legend=c("Obs.", "Pred."), pch=c(1,19),
         cex=c(1,1), bg="white")
  return(pred)
}

# Compute MSE
mse <- function(true, pred) {
  if (length(true) != length(pred)){
    stop("Lengths don't match")
  }
  
  return(sum((true - pred)^2)/length(true))
}


####################################
###### PREDICTION (on a grid) ######
####################################

load("Prediction_Grid.RData")

##### Model 1 Grid Prediction : DEMAND ~ POPULATION + LANES + SUBWAY #####
covars = cbind(rep(1.0, length(Grid.data[,1])), as.matrix(Grid.data[,c("Block.population", "Lane.count", "Dist.metro")]))  # Add the intercept
pred.v1 = predict(sp.model.v1, Grid.data[,c("Longitude", "Latitude")], covars, n.model=1)

##### Model 4 Grid Prediction : DEMAND ~ POPULATION + LANES + SUBWAY + LANDMARKS #####
covars = cbind(rep(1.0, length(Grid.data[,1])), as.matrix(Grid.data[,c("Block.population", "Lane.count", "Dist.metro", "Landmarks")]))  # Add the intercept
pred.v4 = predict(sp.model.v4, Grid.data[,c("Longitude", "Latitude")], covars, n.model=4)


############################################
###### PREDICTION (at station points) ######
############################################
load("Test_centre.RData")

coords = as.matrix(Test_centre[,c("Longitude","Latitude")])
DEMAND = Test_centre$N.Trips

# Plot of the real observed demand
x11()
obs.surf <-
  mba.surf(cbind(coords, DEMAND), no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
points(coords)
contour(obs.surf, add=T)


##### Model 1 Stations Prediction : DEMAND ~ POPULATION + LANES + SUBWAY #####
covars = cbind(rep(1.0, length(Test_centre[,1])), as.matrix(Test_centre[,c("Block.population", "Lane.count", "Dist.metro")]))
pred.v1 = predict(sp.model.v1, Test_centre[,c("Longitude", "Latitude")], covars, n.model=1)
mse1 = mse(Test_centre$N.Trips, rowMeans(pred.v1$p.y.predictive.samples))

##### Model 2 Stations Prediction : DEMAND ~ POPULATION + LANES + SUBWAY + PROXIMITY #####
covars = cbind(rep(1.0, length(Test_centre[,1])), as.matrix(Test_centre[,c("Block.population", "Lane.count", "Dist.metro", "Proximity.score")]))
pred.v2 = predict(sp.model.v2, Test_centre[,c("Longitude", "Latitude")], covars, n.model=2)
mse2 = mse(Test_centre$N.Trips, rowMeans(pred.v2$p.y.predictive.samples))

##### Model 3 Stations Prediction : DEMAND ~ POPULATION + LANES + SUBWAY + PROXIMITY + LANDMARKS #####
covars = cbind(rep(1.0, length(Test_centre[,1])), as.matrix(Test_centre[,c("Block.population", "Lane.count", "Dist.metro", "Proximity.score", "Landmarks")]))
pred.v3 = predict(sp.model.v3, Test_centre[,c("Longitude", "Latitude")], covars, n.model=3)
mse3 = mse(Test_centre$N.Trips, rowMeans(pred.v3$p.y.predictive.samples))

##### Model 4 Stations Prediction : DEMAND ~ POPULATION + LANES + SUBWAY + LANDMARKS #####
covars = cbind(rep(1.0, length(Test_centre[,1])), as.matrix(Test_centre[,c("Block.population", "Lane.count", "Dist.metro", "Landmarks")]))
pred.v4 = predict(sp.model.v4, Test_centre[,c("Longitude", "Latitude")], covars, n.model=4)
mse4 = mse(Test_centre$N.Trips, rowMeans(pred.v4$p.y.predictive.samples))

gof[3,]=c(mse1,mse2,mse3,mse4)
gof
