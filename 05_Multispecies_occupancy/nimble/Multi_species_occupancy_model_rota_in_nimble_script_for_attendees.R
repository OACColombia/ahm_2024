# title: R script for Rota et al. multi-species occupancy model session AHM workshop July 2024
# author: Joshua Twining
# date: 2024-07-02

# install the required packages if you haven't already! 
# install.packages('nimble')
# install.packages('unmarked')
# install.packages('abind')
# install.packages('AHMbook')
# install.packages('coda')
# install.packages('ggplot2')
# install.packages('runjags')
# install.packages('ggfortify')
#   
# load the required packages
library(unmarked)
library(AHMbook)
library(abind)
library(nimble)
library(coda)
library(ggplot2)
library(runjags)
library(ggfortify)

# load the example dataset from the unmarked package
data("MesoCarnivores")

# look at raw data
# check out the dataset
lapply(MesoCarnivores, head) # look at raw data


# pull out the detection/non-detection data for red and coyote
y_rf <- MesoCarnivores$redfox
y_c <- MesoCarnivores$coyote

# put them in a  list
ylist <- list(y_rf, y_c)

# Turn the species detection/non-detection data list in a 3D array (site x occ x species)
y <- abind(ylist, along=3)

str(y)

# Get sample sizes and data dimensions 
nsites <- dim(y)[1]       # 1437 sites
nsurveys <- dim(y)[2]     # 3 surveys
nspec <- dim(y)[3]        # 2 species
ncat <- 2^nspec           # 4 possible community states

# Condense multi-species detection array to be site x survey
ycat <- apply(y, c(1,2), paste, collapse = "")

ycat

#nimble do not like categories, so we convert to numeric indexes
ycat[ycat == "00"] <- 1    # Unoccupied (abbreviated 'U')
ycat[ycat == "10"] <- 2    # Only redfox ('RF') detected
ycat[ycat == "01"] <- 3    # Only coyote ('C') detected
ycat[ycat == "11"] <- 4    # Redfox and coyote ('RFC')


# Convert each column to a numeric: this is our response variable
ycat <- apply(ycat, 2, as.numeric)

ycat

# Pull out site covariate data
siteCovs <- MesoCarnivores$sitecovs

# Prepare site covariates for occupancy
dist_cov <- standardize(siteCovs$Dist_5km) #land use disturbance distance
hdens_cov <- standardize(siteCovs$HDens_5km)
lat_cov <- standardize(siteCovs$Latitude)
long_cov <- standardize(siteCovs$Longitude)
people_cov <- standardize(siteCovs$People_site)

# Prepare observation covariate for detection
table(Trail <- MesoCarnivores$sitecovs[,'Trail'])    # summarize
Trail <- matrix(cbind(Trail, Trail, Trail), ncol = 3)
Trail_cov <- Trail

NimModel <- nimbleCode({
  # --- 'Likelihood' ---
  # (1) Basic hierarchical model: state and observation
  # Latent state model
  for(i in 1:nsites) {
    z[i] ~ dcat(psi[i,1:ncat] ) 
  }
  # Observation model
  for(i in 1:nsites) {
    for(j in 1:nsurveys) {
      y[i, j] ~ dcat(rdm[i, j, 1:ncat  , z[i] ] )
    }
  }
  # (2) Specify linear models for the natural parameters (which will be used to calculate occupancy and detection probs
  # for each category)
  # Linear models for the state model natural parameters
  # natural parameters for states  f1(RF), f2(C), and f12 (RF-C)
  for( i in 1:nsites ) {
    #for the fox
    f1[i] <- beta0RF + betaRF[1] * Dist[i] + betaRF[2] * latitude[i] + betaRF[3] * longitude[i] + betaRF[4] * Hdens[i]  
    #for the coyote
    f2[i] <- beta0C +  betaC[1] * Dist[i] + betaC[2] * latitude[i] + betaC[3] * longitude[i] + betaC[4] * Hdens[i]
    # natural parameters for states f12 (red fox and coyote co-ocurring),
    f12[i] <- beta0RFC + betaRFC[1] * Hdens[i]
    # Linear models for the observation model parameters
    # => Here we could specify detection interactions as well (but we don't)
    for(j in 1:nsurveys){
      # Baseline detection linear predictors
        #detection probability for red fox (camara on a trail or not?)
      rhoRF[i, j] <-alpha0RF + alphaRF[1] * Trail[i, j] 
        #detection probability for coyote  (camara on a trail or not?)
      rhoC[i, j] <- alpha0C +  alphaC[1] * Trail[i, j] 
      # Asymmetric interactions between the two species with the simple detection probability per species
      rhoRFC[i, j] <- rhoRF[i, j] #here we can add covariates
      rhoCRF[i, j] <- rhoC[i, j]
    }
  }
  # (3) Define the latent state vector and the observation matrices
  for( i in 1:nsites ) {
    # calculate latent probabilities for each occupancy state in dcat (this automatically scale by the numerator!)
    psi[i, 1] <- 1 #----------------------------------------| f0 = none present
    psi[i, 2] <- exp( f1[i] ) #-----------------------------| f1=RF
    psi[i, 3] <- exp( f2[i] ) #-----------------------------| f2=C
    psi[i, 4] <- exp( f1[i] + f2[i] + f12[i] ) #------------| f12=RF-C (we can remove interaction removing f12[i])
    for(j in 1:nsurveys){
      # Detection matrix [i = site, j = occasion, x = observed state (os), y = true state (ts))
      # rdm = rho detection matrix.
      # True state = 1(none present - U)
      rdm[i, j, 1, 1] <- 1 #----------------------------------| OS = U
      rdm[i, j, 2, 1] <- 0 #----------------------------------| OS = RF
      rdm[i, j, 3, 1] <- 0 #----------------------------------| OS = C
      rdm[i, j, 4, 1] <- 0 #----------------------------------| OS = RF-C
      # True state = 2 (Red fox only present - RF)
      rdm[i, j, 1, 2] <- 1 #----------------------------------| OS = U (no detect coyote)
      rdm[i, j, 2, 2] <- exp( rhoRF[i, j] ) #-----------------| OS = RF (detect red fox)
      rdm[i, j, 3, 2] <- 0 #----------------------------------| OS = C
      rdm[i, j, 4, 2] <- 0 #----------------------------------| OS = RF-C
      # True state = 3 (Coyote only present)
      rdm[i, j, 1, 3] <- 1 #----------------------------------| OS = U (no detect coyote)
      rdm[i, j, 2, 3] <- 0 #----------------------------------| OS = RF
      rdm[i, j, 3, 3] <- exp( rhoC[i, j] ) #------------------| OS = C (detect coyote)
      rdm[i, j, 4, 3] <- 0 #----------------------------------| OS = RF-C
      # True state = 4 (Red fox and coyote present - all posible)
      rdm[i, j, 1, 4] <- 1 #----------------------------------| OS = U
      rdm[i, j, 2, 4] <- exp( rhoRFC[i, j] ) #----------------| OS = RF
      rdm[i, j, 3, 4] <- exp( rhoCRF[i, j] )# ----------------| OS = C
      rdm[i, j, 4, 4] <- exp( rhoRFC[i, j] +  rhoCRF[i, j]) #-| OS = RF-C
    }
  }
  # --- Priors ---
  # First order psi intercepts
  beta0RF ~ dlogis(0,1) # fo occupancy intercepts
  beta0C ~ dlogis(0,1)
  
  # Second order interaction intercepts
  beta0RFC ~  dlogis(0,1)
  
  # First order psi coefficients
  for(fo_psi in 1:4){ # fo occupancy slopes
    betaRF[fo_psi] ~ dnorm(0, 0.1)
    betaC[fo_psi] ~ dnorm(0, 0.1)
  }
  # Second order psi coefficients
  for(so_psi in 1:1){
    betaRFC[so_psi] ~ dnorm(0, 0.1)
  }
  # First order detection intercept priors (rho)
  alpha0RF ~ dlogis(0,1) # fo detection intercepts
  alpha0C  ~ dlogis(0,1)
  
  # First order detection coefficients
  for(fo_rho in 1:1){ # fo detection slopes
    alphaRF[fo_rho] ~ dnorm(0, 0.1)
    alphaC[fo_rho] ~ dnorm(0, 0.1)
  }
})# end model

# example of plotting prior (if you wondering the difference between the values in the normal distribution in the model above and here are that nimble takes
# mean and precision for its dnorm function, while , ggdistribution takes mean and sd. Where precision is the inverse of the variance, and variance = sd^2
ggdistribution(dnorm, seq(-5, 5, 0.1), mean = 0, sd = 3.16) + theme_classic()

## # # # # # # # # # # # # # # # # # # # # # # # #  WORKSHOP CHALLENGE  # # # # # # # # # # # # # # # # 
# CHALLENGE 01 - plot the logistic prior above dlogis(0,1) 

ggdistribution(dlogis,
               seq(-5, 5, 0.1), #simulated values to plot
               location = 0, #higher probability distribution peak
               scale = 1) + #dufuse
  theme_classic()

# Parameters monitored
parameters <- c('beta0RF', 'beta0C', 'betaRF','betaC',  'beta0RFC', 'betaRFC', 'alphaRF', 'alphaC', 'alpha0RF', 'alpha0C')

# Observation data (nimble data)
Nimdata <- list(y = as.matrix(ycat))

# Provide constants (set all data and model dimensions, and provide all site- and site-by-occasion covariate data)
constants=list(nsites = nsites, nsurveys = nsurveys, ncat=ncat, #dimensions
               #and covariates
               Dist = dist_cov, Trail = Trail_cov,
               Hdens = hdens_cov, latitude = lat_cov, longitude = long_cov)


# Initial values
# get the MAXIMUM possible state across all 3 potential surveys at a site
# returns a site x species matrix
zinit <- apply(y, c(1,3), sum, na.rm = TRUE) #sum across rows and columns
zinit[zinit>1] <- 1 # make binary

# convert to a category
zcat <- apply(zinit, 1, paste, collapse = '')
zcat[zcat == "00"] <- 1 # nobody there
zcat[zcat == "10"] <- 2 # only red fox
zcat[zcat == "01"] <- 3 # only coyote
zcat[zcat == "11"] <- 4 # both red fox coyote

# make numeric again
zcat <- as.numeric(zcat)

# check it out
zcat

# set all initial values
Niminits=list(z=zcat,beta0RF=0,beta0C=0, betaRF=c(0,0,0,0),betaC=c(0,0,0,0), beta0RFC=0, betaRFC=0,
              alphaRF=0, alphaC=0, alpha0RF=0, alpha0C=0)

#thinning rate
nt = 2
nt2 = 5

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, monitors2 = parameters, thin=nt, thin2 = nt2, useConjugacy = TRUE)

# block these due to strong posterior correlation. These are the most correlated parameters, can block more.
conf$addSampler(target = c("beta0C","betaC[4]","beta0RFC", "betaRFC[1]", "beta0RF", "betaRF[4]"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

conf$addSampler(target = c("betaRF[2]","betaRF[3]"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

n.iter <- 5000

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(n.iter,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2=end.time-start.time2 # post-com

#grab the samples
mvSamples = as.matrix(Cmcmc$mvSamples)

# have a look at the chain - set iterations and burn in for this mini run
n.iter=2500
n.burn = 500

# plot it
plot(mcmc(mvSamples[n.burn:nrow(mvSamples),]))


# Package it up and run three chains
# n.chains = 3
# chains = vector("list", n.chains)
# chains2 = vector("list", n.chains)
# for(chain in 1:n.chains){
#   
#   # Parameters monitored
#   parameters <- c('beta0RF', 'beta0C', 'betaRF','betaC',  'beta0RFC', 'betaRFC', 'alphaRF', 'alphaC', 'alpha0RF', 'alpha0C')
#   
#   # Observation data (nimble data)
#   Nimdata <- list(y = as.matrix(ycat))
#   
#   # Provide constants (set all data and model dimensions, and provide all site- and site-by-occasion covariate data)
#   constants=list(nsites = nsites, nsurveys = nsurveys, ncat=ncat, Dist = dist_cov, Trail = Trail_cov,
#                  Hdens = hdens_cov, latitude = lat_cov, longitude = long_cov)
#   
#   # Initial values
#   # get the maximum possible state across all 3 potential surveys at a site
#   # returns a site x species matrix
#   zinit <- apply(y, c(1,3), sum, na.rm = TRUE)
#   zinit[zinit>1] <- 1 # make binary
#   
#   # convert to a category
#   zcat <- apply(zinit, 1, paste, collapse = '')
#   zcat[zcat == "00"] <- 1 # nobody there
#   zcat[zcat == "10"] <- 2 # only red fox
#   zcat[zcat == "01"] <- 3 # only coyote
#   zcat[zcat == "11"] <- 4 # both red fox coyote
#   
#   # make numeric again
#   zcat <- as.numeric(zcat)
#   
#   # check it out
#   zcat
#   
#   # set all initial values
#   Niminits=list(z=zcat,beta0RF=0,beta0C=0, betaRF=c(0,0,0,0),betaC=c(0,0,0,0), beta0RFC=0, betaRFC=0,
#                 alphaRF=0, alphaC=0, alpha0RF=0, alpha0C=0)
#   
#   #thinning rate
#   nt = 4
#   nt2 = 200
#   
#   # Build the model, configure the mcmc, and compile
#   start.time<-Sys.time()
#   Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
#                         inits=Niminits)
#   conf <- configureMCMC(Rmodel,monitors=parameters, monitors2 = parameters, thin=nt, thin2 = nt2, useConjugacy = TRUE)
#   
#   # block these due to strong posterior correlation. These are the most correlated parameters, can block more. 
#   conf$addSampler(target = c("beta0C","betaC[4]","beta0RFC", "betaRFC[1]", "beta0RF", "betaRF[4]"),
#                   type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)
#   
#   conf$addSampler(target = c("betaRF[2]","betaRF[3]"),
#                   type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)
#   
#   
#   # Build and compile
#   Rmcmc <- buildMCMC(conf)
#   # runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
#   Cmodel <- compileNimble(Rmodel)
#   Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
#   
#   n.iter <- 100000
#   
#   # Run the model.
#   start.time2<-Sys.time()
#   Cmcmc$run(n.iter,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
#   end.time<-Sys.time()
#   time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
#   time2=end.time-start.time2 # post-com
#   
#   #get the chains
#   mvSamples = as.matrix(Cmcmc$mvSamples)
#   mvSamples2 = as.matrix(Cmcmc$mvSamples2)
#   
#   chains[[chain]]=mvSamples
#   chains2[[chain]]=mvSamples2
#   
# }

load('05_Multispecies_occupancy/nimble/Two-species_Rota_cooccurence_model_Mesocarnivores_casestudy_100k_chains.RData')
load('05_Multispecies_occupancy/nimble/Two-species_Rota_cooccurence_model_Mesocarnivores_casestudy_100k_chains2.RData')

# combine the chains into an mcmc list
n.iter=25000
n.burn = 20000
a=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
            mcmc(chains[[2]][n.burn:n.iter,]),
            mcmc(chains[[3]][n.burn:n.iter,]))

# gelman-rubin diagnostics
gelman_a = gelman.diag(a)
gelman_a

# plot and do visual checks of traceplots for convergence
plot(a)

# look at model summary
summary(a)
#see the 2. Quantiles - do not overlap zero? important parameter. 
  #Important, retake the question (betaC[4])

## # # # # # # # # # # # # # # # # # # # # # # # #  WORKSHOP CHALLENGE  # # # # # # # # # # # # # # # # 
# CHALLENGE 2-4: 
#What is the relationship between Red fox occupancy and disturbance?
  # betaRF[1] is positive and does not overlap 0
# What is the relationship between latitude and coyote occupancy?
  # betaC[2] overlaps with 0, no evidence of relationship
# Is there an interaction between red foxes and coyotes? And if so, what is it? 
  # beta0RFC is positive and does not overlap 0


# For conducting posterior computations we will want to use some heavily thinned chains
# heavily thinned chain
n.iter=500
n.burn = 400
# colnames(mvSamples)
b=mcmc.list(mcmc(chains2[[1]][n.burn:n.iter,]),
            mcmc(chains2[[2]][n.burn:n.iter,]),
            mcmc(chains2[[3]][n.burn:n.iter,]))

b=runjags::combine.mcmc(b)


######### Producing marginal occupancy estimates ##########
# This example produces marginal psi ests for red fox as function of disturbance

# number of posterior samples
n.samples = nrow(b)

# range of disturbance covariate
r<- range(dist_cov)

# create sequence along range of disturbance
dist_data <- seq(r[1], r[2], length.out=500)

# set all other covs at mean of 0 (it was standardized before)
hdens_mean = 0
latitude_mean = 0
longitude_mean = 0

# create matrices to stick estimates in
f0 = matrix(1, n.samples, length(dist_data))
f1 = matrix(NA, n.samples, length(dist_data)) 
f2 = matrix(NA, n.samples, length(dist_data))   
f12 = matrix(NA, n.samples, length(dist_data))   

# Sample from posterior for range of disturbance values sampled in dataset 
for (i in 1:n.samples){ #every posterior sample
  for (j in 1:length(dist_data)){ #every of disturbance distance
    # create the linear predictors for first order natural parameters
      #Red Fox 
    f1[i,j] = b[,'beta0RF'][[i]] + b[,'betaRF[1]'][[i]] * dist_data[j] + 
      b[,'betaRF[2]'][[i]] * latitude_mean + b[,'betaRF[3]'][[i]] * longitude_mean + b[,'betaRF[4]'][[i]] * hdens_mean #all these means are zero
      #Coyote
    f2[i,j] = b[,'beta0C'][[i]] +  b[,'betaC[1]'][[i]] * dist_data[j] + 
      b[,'betaC[2]'][[i]] * latitude_mean + b[,'betaC[3]'][[i]] * longitude_mean + b[,'betaC[4]'][[i]] * hdens_mean #all these means are zero
    # natural parameters for states f12 (RF-C)
    f12[i,j] =  b[,'beta0RFC'][[i]] + b[,'betaRFC[1]'][[i]] * hdens_mean
  }}

# create a 3D array for each state (unoccupied, red fox only [f1], coyote only [f2], and red fox and coyote [f1+f2+f12])
# Structure should columns are posterior samples, rows are different values of covariate x, and slices are states (ncat=4)
natparams_3D <- array(NA, dim=c(n.samples, length(dist_data), ncat))

natparams_3D[,,1] <- 1
natparams_3D[,,2] <- exp(f1)
natparams_3D[,,3] <- exp(f2)
natparams_3D[,,4] <- exp(f12+ f1 + f2)

# check it by pulling out a single state (e.g. f1)
natparams_f1 <- natparams_3D[,,2]
str(natparams_f1)

# now create a 3D array where we will finish conversions to probabilities
psi <- array(NA, dim=c(n.samples, length(dist_data), ncat))

# scale exponents of each state so that sum of probabilities across categories = 1
for (i in 1:ncol(natparams_3D)){
  for (j in 1:nrow(natparams_3D)){
    psi[j,i, ] <- natparams_3D[j,i,(1:ncat)] / sum(natparams_3D[j,i,(1:ncat)])
  }}

# create 2D array to store the products of our calculations
redfoxmarg_preds <- array(NA, dim=c(n.samples, length(dist_data)))

# check structure
str(redfoxmarg_preds)

# calculate marginal psi ests (we ignore other species, so include all states where red fox occurs, that is 2 and 4)

redfoxmarg_preds <- psi[,,2] + psi[,,4]

# get posterior means
redfoxmarg_dist_pred <- colMeans(redfoxmarg_preds) 

# calculate 95% confidence intervals
redfoxmarg_dist_CIs <- apply(redfoxmarg_preds,2,quantile, c(0.025,0.975), na.rm=TRUE)

# package up the red fox marginal psi ests, with the sequenece of disturbance values
redfoxmargs <- data.frame(disturbance = dist_data, 
                          predicted = redfoxmarg_dist_pred, 
                          LCI = redfoxmarg_dist_CIs[1,],
                          UCI = redfoxmarg_dist_CIs[2,])

# plot it! 
ggplot(redfoxmargs, aes(x = disturbance, y = predicted)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.1, linetype = "dashed")+
  geom_path(size=1)+ 
  ylab("Predicted probability of occupancy")+ 
  xlab("Disturbance (scaled and centered)")+
  theme_classic()



## # # # # # # # # # # # # # # # # # # # # # # # #  WORKSHOP CHALLENGE  # # # # # # # # # # # # # # # # 

# CHALLENGE 5 - Produce marginal estimates of coyote occupancy over a range of disturbance values 

coyotemarg_preds <- psi[,,3] + psi[,,4]

# get posterior means
coyotemarg_dist_pred <- colMeans(coyotemarg_preds) 

# calculate 95% confidence intervals
coyotemarg_dist_CIs <- apply(coyotemarg_preds,2,quantile, c(0.025,0.975), na.rm=TRUE)

# package up the red fox marginal psi ests, with the sequenece of disturbance values
coyotemargs <- data.frame(disturbance = dist_data, 
                          predicted = coyotemarg_dist_pred, 
                          LCI = coyotemarg_dist_CIs[1,],
                          UCI = coyotemarg_dist_CIs[2,])

# plot it! 
ggplot(coyotemargs, aes(x = disturbance, y = predicted)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.1, linetype = "dashed")+
  geom_path(size=1)+ 
  labs(x = "Predicted probability of occupancy",
       y = "Disturbance (scaled and centered)",
       title = "For coyote")+
  theme_classic()

# CHALLENGE 6 - Produce marginal estimates of red fox and coyote occupancy over a range of latitudes   

# range of disturbance covariate
r<- range(lat_cov)

# create sequence along range of disturbance
lat_data <- seq(r[1], r[2], length.out=500)

# set all other covs at mean of 0 (it was standardized before)
hdens_mean = 0
distance_mean = 0
longitude_mean = 0

# create matrices to stick estimates in
f0 = matrix(1, n.samples, length(lat_data))
f1 = matrix(NA, n.samples, length(lat_data)) 
f2 = matrix(NA, n.samples, length(lat_data))   
f12 = matrix(NA, n.samples, length(lat_data))   

# Sample from posterior for range of latitude values sampled in dataset 
for (i in 1:n.samples){ #every posterior sample
  for (j in 1:length(lat_data)){ #every of disturbance distance
    # create the linear predictors for first order natural parameters
    #Red Fox 
    f1[i,j] = b[,'beta0RF'][[i]] + b[,'betaRF[1]'][[i]] * lat_data[j] + 
      b[,'betaRF[2]'][[i]] * distance_mean + b[,'betaRF[3]'][[i]] * longitude_mean + b[,'betaRF[4]'][[i]] * hdens_mean #all these means are zero
    #Coyote
    f2[i,j] = b[,'beta0C'][[i]] +  b[,'betaC[1]'][[i]] * lat_data[j] + 
      b[,'betaC[2]'][[i]] * distance_mean + b[,'betaC[3]'][[i]] * longitude_mean + b[,'betaC[4]'][[i]] * hdens_mean #all these means are zero
    # natural parameters for states f12 (RF-C)
    f12[i,j] =  b[,'beta0RFC'][[i]] + b[,'betaRFC[1]'][[i]] * hdens_mean
  }}

# create a 3D array for each state (unoccupied, red fox only [f1], coyote only [f2], and red fox and coyote [f1+f2+f12])
# Structure should columns are posterior samples, rows are different values of covariate x, and slices are states (ncat=4)
natparams_3D <- array(NA, dim=c(n.samples, length(lat_data), ncat))

natparams_3D[,,1] <- 1
natparams_3D[,,2] <- exp(f1)
natparams_3D[,,3] <- exp(f2)
natparams_3D[,,4] <- exp(f12+ f1 + f2)

# check it by pulling out a single state (e.g. f1)
natparams_f1 <- natparams_3D[,,2]
str(natparams_f1)

# now create a 3D array where we will finish conversions to probabilities
psi <- array(NA, dim=c(n.samples, length(lat_data), ncat))

# scale exponents of each state so that sum of probabilities across categories = 1
for (i in 1:ncol(natparams_3D)){
  for (j in 1:nrow(natparams_3D)){
    psi[j,i, ] <- natparams_3D[j,i,(1:ncat)] / sum(natparams_3D[j,i,(1:ncat)])
  }}

# create 2D array to store the products of our calculations
redfoxmarg_preds <- array(NA, dim=c(n.samples, length(lat_data)))

# check structure
str(redfoxmarg_preds)

# calculate marginal psi ests (we ignore other species, so include all states where red fox occurs, that is 2 and 4)

redfoxmarg_preds <- psi[,,2] + psi[,,4]

# get posterior means
redfoxmarg_lat_pred <- colMeans(redfoxmarg_preds) 

# calculate 95% confidence intervals
redfoxmarg_lat_CIs <- apply(redfoxmarg_preds,2,quantile, c(0.025,0.975), na.rm=TRUE)

# package up the red fox marginal psi ests, with the sequenece of disturbance values
redfoxmargs <- data.frame(latitude = lat_data, 
                          predicted = redfoxmarg_lat_pred, 
                          LCI = redfoxmarg_lat_CIs[1,],
                          UCI = redfoxmarg_lat_CIs[2,])

# plot it! 
ggplot(redfoxmargs, aes(x = latitude, y = predicted)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.1, linetype = "dashed")+
  geom_path(size=1)+ 
  labs(y = "Predicted probability of occupancy",
      x = "Latitude (scaled and centered)",
      title = "For red fox")+
  theme_classic()

#Coyote
coyotemarg_preds <- psi[,,3] + psi[,,4]

# get posterior means
coyotemarg_lat_pred <- colMeans(coyotemarg_preds) 

# calculate 95% confidence intervals
coyotemarg_lat_CIs <- apply(coyotemarg_preds,2,quantile, c(0.025,0.975), na.rm=TRUE)

# package up the red fox marginal psi ests, with the sequenece of disturbance values
coyotemargs <- data.frame(latitude = lat_data, 
                          predicted = coyotemarg_lat_pred, 
                          LCI = coyotemarg_lat_CIs[1,],
                          UCI = coyotemarg_lat_CIs[2,])

# plot it! 
ggplot(coyotemargs, aes(x = latitude, y = predicted)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.1, linetype = "dashed")+
  geom_path(size=1)+ 
  labs(y = "Predicted probability of occupancy",
       x = "Latitude (scaled and centered)",
       title = "For coyote")+
  theme_classic()

######### Producing conditional occupancy estimates ##########
# this example produces estimates of probability of coyote occupancy conditional on the presence/absence of red fox
# as a function of housing density

# number of posterior samples
n.samples = nrow(b)

# range of housing density covariate
r<- range(hdens_cov)

# create sequence along range of housing density
hdens_data <- seq(r[1], r[2], length.out=500)

# set all other covs at mean
dist_mean = 0
latitude_mean = 0
longitude_mean = 0

# create matrices to stick posterior estimates in
f0 = matrix(1, n.samples, length(hdens_data))
f1 = matrix(NA, n.samples, length(hdens_data)) 
f2 = matrix(NA, n.samples, length(hdens_data))   
f12 = matrix(NA, n.samples, length(hdens_data))   

# Sample from posterior for range of housing density values sampled in dataset 
for (i in 1:n.samples){
  for (j in 1:length(hdens_data)){
    # create the linear predictors for first order natural parameters
    f1[i,j] = b[,'beta0RF'][[i]] + b[,'betaRF[1]'][[i]] * dist_mean + b[,'betaRF[2]'][[i]] * latitude_mean + b[,'betaRF[3]'][[i]] * longitude_mean + b[,'betaRF[4]'][[i]] * hdens_data[j]
    f2[i,j] = b[,'beta0C'][[i]] +  b[,'betaC[1]'][[i]] * dist_mean + b[,'betaC[2]'][[i]] * latitude_mean + b[,'betaC[3]'][[i]] * longitude_mean + b[,'betaC[4]'][[i]] * hdens_data[j]
    # linear predictors for second order natural parameters 
    f12[i,j] = b[,'beta0RFC'][[i]] + b[,'betaRFC[1]'][[i]] * hdens_data[j]
  }}

# create a 3D array for each state (unoccupied, red fox only [f1], coyote only [f2], and red fox and coyote [f1+f2+f12])
# Structure should be columns are posterior samples, rows are different values of covariate x, and slices (3rd dimension) are states (ncat=4)
natparams_3D <- array(NA, dim=c(n.samples, length(hdens_data), ncat))

natparams_3D[,,1] <- 1
natparams_3D[,,2] <- exp(f1)
natparams_3D[,,3] <- exp(f2)
natparams_3D[,,4] <- exp(f12+ f1 + f2)

# check it by pulling out a single natural parameter
natparam_f1 <- natparams_3D[,,2]
str(natparam_f1)

# now create a 3D array where we will finish conversions to probabilities
psi <- array(NA, dim=c(n.samples, length(hdens_data), ncat))

# scale exponents of each state so that sum of probabilities across categories = 1
for (i in 1:ncol(natparams_3D)){
  for (j in 1:nrow(natparams_3D)){
    psi[j,i, ] <- natparams_3D[j,i,(1:ncat)] / sum(natparams_3D[j,i,(1:ncat)])
  }}


#  red fox probability of occupancy conditional on the absence of coyotes
red.alone.samples <- psi[,,2]/(psi[,,1]+ psi[,,2])

# calculate posterior means and 95% credible intervals
fox.alone.mean <- colMeans(red.alone.samples) 
fox.alone.CIs <- apply(red.alone.samples,2,quantile, c(0.025,0.975), na.rm=TRUE)

# package it up into df
fox.cond.alone <- data.frame(housing_density = hdens_data, 
                             predicted = fox.alone.mean, 
                             LCI = fox.alone.CIs[1,],
                             UCI = fox.alone.CIs[2,],
                             coyote = "absent")

#  red fox probability of occupancy conditional on the presence of coyotes
redfox.given.coyote.samples <- psi[,,4]/(psi[,,3]+psi[,,4])

# calculate posterior means and 95% credible intervals
redfox.given.coyote.mean <-colMeans(redfox.given.coyote.samples) 
redfox.given.coyote.CI <-apply(redfox.given.coyote.samples,2,quantile, c(0.025,0.975), na.rm=TRUE)

# package it up into df
redfox.cond.coyote <- data.frame(housing_density = hdens_data, 
                                 predicted = redfox.given.coyote.mean, 
                                 LCI = redfox.given.coyote.CI[1,],
                                 UCI = redfox.given.coyote.CI[2,],
                                 coyote = "present")
# merge the dfs
redfox.cond.psi <- rbind(fox.cond.alone, redfox.cond.coyote)

# plot the conditional estimates
ggplot(redfox.cond.psi, aes(x = housing_density, y = predicted, fill = coyote, colour = coyote))+
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.1, linetype = "dashed")+
  geom_path(size=1)+ ylab("Predicted red fox probability of occupancy")+ xlab("Housing density (scaled and centered)")+
  theme_classic()


# why the uncertainty?
hist(hdens_cov)

idx<- which(hdens_cov>2)

ycheck <- ycat[idx,]

ycheck

## # # # # # # # # # # # # # # # # # # # # # # # #  WORKSHOP CHALLENGE  # # # # # # # # # # # # # # # # 
# CHALLENGE 7 - produce occupancy estimates for coyote conditional on red fox presence/abscence as a function of housing density
coy.alone.samples <- psi[,,3]/(psi[,,1]+ psi[,,3])

# calculate posterior means and 95% credible intervals
coy.alone.mean <- colMeans(coy.alone.samples) 
coy.alone.CIs <- apply(coy.alone.samples,2,quantile, c(0.025,0.975), na.rm=TRUE)

# package it up into df
coy.cond.alone <- data.frame(housing_density = hdens_data, 
                             predicted = coy.alone.mean, 
                             LCI = coy.alone.CIs[1,],
                             UCI = coy.alone.CIs[2,],
                             fox = "absent")

#  red fox probability of occupancy conditional on the presence of coyotes
coyote.given.redfox.samples <- psi[,,4]/(psi[,,2]+psi[,,4])

# calculate posterior means and 95% credible intervals
coyote.given.redfox.mean <-colMeans(coyote.given.redfox.samples) 
coyote.given.redfox.CI <-apply(coyote.given.redfox.samples,2,quantile, c(0.025,0.975), na.rm=TRUE)

# package it up into df
coyote.cond.redfox <- data.frame(housing_density = hdens_data, 
                                 predicted = coyote.given.redfox.mean, 
                                 LCI = coyote.given.redfox.CI[1,],
                                 UCI = coyote.given.redfox.CI[2,],
                                 fox = "present")
# merge the dfs
coyote.cond.psi <- rbind(coy.cond.alone, coyote.cond.redfox)

# plot the conditional estimates
ggplot(coyote.cond.psi, aes(x = housing_density, y = predicted, fill = fox, colour = fox))+
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.1, linetype = "dashed")+
  geom_path(size=1)+ 
  labs(y = "Predicted red fox probability of occupancy", 
       x = "Housing density (scaled and centered)",
       title = "Coyote conditional on red fox")+
  theme_classic()


# CHALLENGE 7B. Red fox probability of occupancy conditional on the absence of coyotes as a function of latitude

# range of latitude covariate
r<- range(lat_cov)

# create sequence along range of housing density
lat_data <- seq(r[1], r[2], length.out=500)

# set all other covs at mean
dist_mean = 0
housing_mean = 0
longitude_mean = 0

# create matrices to stick posterior estimates in
f0 = matrix(1, n.samples, length(lat_data))
f1 = matrix(NA, n.samples, length(lat_data)) 
f2 = matrix(NA, n.samples, length(lat_data))   
f12 = matrix(NA, n.samples, length(lat_data))   

# Sample from posterior for range of housing density values sampled in dataset 
for (i in 1:n.samples){
  for (j in 1:length(lat_data)){
    # create the linear predictors for first order natural parameters
    f1[i,j] = b[,'beta0RF'][[i]] + b[,'betaRF[1]'][[i]] * dist_mean + b[,'betaRF[2]'][[i]] * lat_data[j] + b[,'betaRF[3]'][[i]] * longitude_mean + b[,'betaRF[4]'][[i]] * housing_mean
    f2[i,j] = b[,'beta0C'][[i]] +  b[,'betaC[1]'][[i]] * dist_mean + b[,'betaC[2]'][[i]] * lat_data[j] + b[,'betaC[3]'][[i]] * longitude_mean + b[,'betaC[4]'][[i]] * housing_mean
    # linear predictors for second order natural parameters 
    f12[i,j] = b[,'beta0RFC'][[i]] + b[,'betaRFC[1]'][[i]] * lat_data[j]
  }}

# create a 3D array for each state (unoccupied, red fox only [f1], coyote only [f2], and red fox and coyote [f1+f2+f12])
# Structure should be columns are posterior samples, rows are different values of covariate x, and slices (3rd dimension) are states (ncat=4)
natparams_3D <- array(NA, dim=c(n.samples, length(lat_data), ncat))

natparams_3D[,,1] <- 1
natparams_3D[,,2] <- exp(f1)
natparams_3D[,,3] <- exp(f2)
natparams_3D[,,4] <- exp(f12+ f1 + f2)

# check it by pulling out a single natural parameter
natparam_f1 <- natparams_3D[,,2]
str(natparam_f1)

# now create a 3D array where we will finish conversions to probabilities
psi <- array(NA, dim=c(n.samples, length(lat_data), ncat))

# scale exponents of each state so that sum of probabilities across categories = 1
for (i in 1:ncol(natparams_3D)){
  for (j in 1:nrow(natparams_3D)){
    psi[j,i, ] <- natparams_3D[j,i,(1:ncat)] / sum(natparams_3D[j,i,(1:ncat)])
  }}


#  red fox probability of occupancy conditional on the absence of coyotes
red.alone.samples <- psi[,,2]/(psi[,,1]+ psi[,,2])

# calculate posterior means and 95% credible intervals
fox.alone.mean <- colMeans(red.alone.samples) 
fox.alone.CIs <- apply(red.alone.samples,2,quantile, c(0.025,0.975), na.rm=TRUE)

# package it up into df
fox.cond.alone <- data.frame(latitude = lat_data, 
                             predicted = fox.alone.mean, 
                             LCI = fox.alone.CIs[1,],
                             UCI = fox.alone.CIs[2,],
                             coyote = "absent")

#  red fox probability of occupancy conditional on the presence of coyotes
redfox.given.coyote.samples <- psi[,,4]/(psi[,,3]+psi[,,4])

# calculate posterior means and 95% credible intervals
redfox.given.coyote.mean <-colMeans(redfox.given.coyote.samples) 
redfox.given.coyote.CI <-apply(redfox.given.coyote.samples,2,quantile, c(0.025,0.975), na.rm=TRUE)

# package it up into df
redfox.cond.coyote <- data.frame(latitude = lat_data, 
                                 predicted = redfox.given.coyote.mean, 
                                 LCI = redfox.given.coyote.CI[1,],
                                 UCI = redfox.given.coyote.CI[2,],
                                 coyote = "present")
# merge the dfs
redfox.cond.psi <- rbind(fox.cond.alone, redfox.cond.coyote)

# plot the conditional estimates
ggplot(redfox.cond.psi, aes(x = latitude, y = predicted, fill = coyote, colour = coyote))+
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.1, linetype = "dashed")+
  geom_path(size=1)+ 
  labs(y = "Predicted red fox probability of occupancy",
       x = "Latitude (scaled and centered)",
       title = "conditional probability of red fox through latitude")+
  theme_classic()

# CHALLENGE 8 (the final boss) - specify a three species multi-species occupancy model including coyote, bobcat, and red fox





