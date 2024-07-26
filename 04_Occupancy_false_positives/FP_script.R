# Part 1: basic Royle-Link model, and demonstration of multi-modality issue 
#   (finding the wrong mode with bad starting values)
#
# Simulate data from the Royle-Link model (occupancy with false positives)

set.seed(1) # Initialize RNGs
nsites <- 200 # number of sites (i = 1, ..., nsites=M)
nsurveys <- 7 # number of visits (j = 1, ..., nsurveys=J)
psi <- 0.6    # expected occupancy probability
p <- 0.7      # detection probability (p_11)
fp <- 0.05    # false-positive error probability (p_10)

# Simulate occupancy states and encounter histories
z <- rbinom(nsites, 1, psi) # occupancy states
y <- matrix(NA, nrow = nsites, ncol = nsurveys) # empty matrix for detections
for(i in 1:nsites){
  pr_yequals1 <- p*z[i] + fp*(1 - z[i]) # p11 + p10
  y[i,] <- rbinom(nsurveys, 1, pr_yequals1) # realized observations
}

 
# Number of false-positive detections per occasion
apply(y[z==0,]>0, 2, sum)

# [1] 8 4 3 4 1 8 7

# Number of false-negative detections per occasion
apply(y[z==1,]==0, 2, sum)
# [1] 39 32 33 39 44 34 42

# Fit the model in unmarked. Data contaminated with false positives is type 2 
#  data according to the unmarked conventions. So here we have 7 occasions of Type 2 data:

#Type 1 = perfect occupancy data
#Type 2 = some false positive detections
#type 3 = observations are assumed to include both certain detections and uncertain detections that may include false positive detections

type <- c(0, 7, 0) #columns of the y columns with false positives

# Build the unmarkedFrame using the unmarkedFrameFP constructor function
library(unmarked)
summary(umf <- unmarkedFrameOccuFP(y = y, type = type))  

# For fitting a normal occupancy model you would just do this:
umf.occ<- unmarkedFrameOccu(y=y)

#
# Using bad starting values you can see that the wrong mode of the likelihood can be found:
#

# Using these starting values will find the wrong mode.
  # Order is psi, p, fp   <- the order is important, otherwise you can't understand what's going on
largerp10 <- qlogis(c(0.5, 0.1, 0.7))  

# These starting values are consistent with our preferred constraint
largerp11 <- qlogis(c(0.5, 0.7, 0.1))


(m1 <- occuFP(detformula = ~1, # model for p_11
	FPformula = ~1,          # model for p_10
	stateformula = ~1,       # model for psi
	data = umf,              # umarkedFrameOccuFP object
	starts = largerp11) )    # add p_11 > p_10 constraint

plogis(coef(m1)) #similar to the simulated data

( m2 <- occuFP(detformula = ~1, # model for p_11
               FPformula = ~1, # model for p_10
               stateformula = ~1, # model for psi
               data = umf, # umarkedFrameOccuFP object
               starts = largerp10) ) # add p_11 < p_10 constraint

plogis(coef(m2)) #far away from simulated data: psi=0.6, p=0.7, fp=0.05!!!!

# Notice that the two models have identical AIC values, and that the values of p11 and p10 flip and psi becomes 1-psi

m1@AIC ; m2@AIC

#[1] 1550.633
#[1] 1550.633

cbind("m1" = plogis( coef(m1) ), "m2" = plogis( coef(m2)))

# You can see that ignoring false positives increases the estimate of psi. 
# In general estimate of psi is biased in this case, you would have to do a simulation to verify that

umf.occ<- unmarkedFrameOccu(y=y)
 occu(~1 ~1, data =umf.occ)
 
 #occupancy probability is overestimated
 plogis(1.06)
 
 #detection probability is underestimated
 plogis(0.277)
 
#
# Part 2: General model of Miller et al.
#
# Miller model with combination of T1 and T3 data. T3 data have "random confirmation", T1 data are normal occupancy data.
# To simulate T1 data we simulate data with false positives (as above) but then randomly declare some of them, with z = 1, to be confirmed
#
# You can think of T1 and T3 data as arising from different survey methods (e.g., camera trapping produces
#  occupancy data and T3 data might be sightings, whhich are confirmed if made by the state biologist).
# So here we will include a "method effect" on detection probability.
#

# Set parameter values of the simulation
set.seed(129) # RNG seed
nsites <- 200 # number of sites
nsurveys1 <- 3 # number of occasions with Type 1 data
nsurveys2 <- 4 # number of occasions with Type 3 data
psi <- 0.6 # expected proportion of are occupied
p <- c(0.7,0.5) # detection prob of method 1 and method 2
fp <- 0.05 # false-positive error probability (p_10)
b <- 0.2 # probability y is recorded as certain


# Simulate the occupancy states and data
z <- rbinom(nsites, 1, psi)
y <- matrix(NA, nrow = nsites, ncol = nsurveys1 + nsurveys2)
for(i in 1:nsites){
p1 <- p[1]*z[i] # certainly detection (method 1)
p2 <- p[2]*z[i] + fp*(1-z[i]) # uncertainly detection (method 2)
y[i,1:3] <- rbinom(nsurveys1, 1, p1) # simulate method 1 data
y[i,4:7] <- rbinom(nsurveys2, 1, p2) # simulate method 2 data

# Now introduce certain observations:
pr.certain <- z[i] * y[i,4:7] * b     # This looks confusing but it allows only z = 1 and y = 1 obs to be confirmed
  #if the site is occupy and is subject to false positive, flip a coin (below)
y[i, 4:7] <- y[i, 4:7] + rbinom(4, 1, pr.certain)
}

# Make a covariate to distinguish between the two methods
Method <- matrix(c(rep("1", 3), #collected with method 1
                   rep("2", 4)), #collected with method 2
                 nrow = nsites, 
                 ncol = nsurveys1 + nsurveys2, byrow = TRUE)

type <- c(nsurveys1, 0, nsurveys2)

summary(umf2 <- unmarkedFrameOccuFP(y = y, obsCovs = list(Method = Method), type = type)) 

#
# Note: starting values are not specified but , in general, you should specify starts that make sense (low fp)
#

(m2 <- occuFP(detformula = ~ -1 + Method, 
              FPformula = ~ 1, 
              Bformula = ~ 1, #indicating type3 data in unmark
              stateformula = ~ 1, 
              data = umf2) )

# Coefficients on the link (= "beta") scale
coef(m2)
 
# Coefficients on the probability (="real") scale
pred.df <- data.frame(Method = c("1", "2"))
round(rbind(
  "det" = predict(m2, type = 'det', newdata = pred.df),           #probability of detection
  "fp" = predict(m2, type = 'fp', newdata = pred.df[1,,drop=F]),  #probability of false positive
  "b" = predict(m2, type = 'b', newdata = pred.df[1,,drop=F]),    #probability of confirmation from "experts"
  "state" = predict(m2, type = 'state', newdata = pred.df[1,,drop=F])),
  3) 
 

#####
#
# Part 3:
#
# General case with all 3 data types
#
#####

# Simulation settings
set.seed(2019)  # RNG seed
nsites <- 200   # number of sites
nsurveys <- 7   # number of occasions
habitat <- rnorm(nsites) # Some (continuous) habitat descriptor
# Simulate the occupancy states and data
alpha0 <- 0   # Intercept...
alpha1 <- 1   # ... and slope of psi-habitat regression
psi <- plogis(alpha0 + alpha1*habitat) # Occupancy
z <- rbinom(nsites, 1, psi) # Latent p/a states
y <- matrix(0,nsites, nsurveys)
p <- c(0.7, 0.5) # method 2 will have a lower p
b <- 0.5 # probability that a observed positive is determined to be certain
fp <- 0.05 # False-positive prob.


# Simulate data of all 3 types. Note p differs between occ 1-2 and 3-7.
# False positives occur in occasions 3-7 but in occasion 7 there are some confirmed positives

for(i in 1:nsites){
 # Normal occupancy data
 y[i, 1:2] <- rbinom(2, 1, p[1]*z[i])
 # False-positives mixed in
 y[i, 3:6] <- rbinom(4, 1, p[2]*z[i] + fp*(1-z[i]))
 # Type 3 observations are occupancy data contaminated with false positives but then
 # we identify some of them as true
 y[i, 7] <- rbinom(1, 1, p[2]*z[i] + fp*(1-z[i]))
}

# Next we set some of the detections to confirmed positives
true.positives <- z==1 & y[,7]==1
confirmed <- (rbinom(nsites, 1, b) == 1) & true.positives
y[confirmed, 7] <- 2
# Make a covariate to distinguish between the two methods
Method <- matrix(c(rep("1", 2), rep("2", 5)), nrow = nsites, ncol = 7, byrow = TRUE)
# Type indicates a mix of all 3 data types
type <- c(2, 4, 1)

# Same covariate structure as before
siteCovs <- data.frame(habitat = habitat)
obsCovs <- list(Method = Method)
summary(umf1 <- unmarkedFrameOccuFP(y, siteCovs = siteCovs, obsCovs = obsCovs, type = type)) 

# fp starting value should be small (-1 here).

# Note: last parameter in this model is "Pcertain"
( m3 <- occuFP(detformula = ~ -1 + Method, #factor variables, two slots in "starts"
               FPformula = ~1, #
               Bformula = ~1, 
               stateformula = ~ habitat, 
               data = umf1, 
               starts=c(0, 0, 0, 0, -1, 0)) )
               # starting values order:   psi parms, p parms, fp parms, b parms.   

#Class exercise with Lynx######
library(AHMbook);
data("EurasianLynx")
str(lynx <- EurasianLynx)

#types
table(lynx$type)

# Add the columns we need for analysis in unmarked
lynx$occ.1 <- 1
lynx$occ.2 <- 2
lynx$occ.3 <- 3
lynx$sYear <- standardize(lynx$Year)

# Extract the type 1 and type 2 data separately and bind them together
lynx1 <- lynx[lynx[,"type"] == "certain", ]
lynx2 <- lynx[lynx[,"type"] == "uncertain", ]
lynx <- cbind(lynx1[,c(2,3:5)], lynx2[,3:5] ) #concatenated the two types to organize the observation events?
colnames(lynx) <- c("site.nr", "y.1", "y.2", "y.3", "y.4", "y.5", "y.6")

occ <- cbind(lynx1[,c("occ.1","occ.2","occ.3")], lynx2[, c("occ.1", "occ.2", "occ.3")])
colnames(occ) <- c("occ.1", "occ.2", "occ.3", "occ.4", "occ.5", "occ.6")
lynx <- cbind(lynx,lynx1[, c("Year", "sYear", "Cntry")])

# Make the false-positive unmarkedFrame. Be sure to indicate type!
y <- lynx[,paste0("y.", 1:6)]
siteCovs <- lynx[, c("sYear", "Year", "Cntry")]
obsCovs <- list(occ = occ)
type = c(nrow(lynx1), #type 1 - certain
         nrow(lynx2), #type 2 - uncertain
         0) 

summary(lynxumf <- unmarkedFrameOccuFP(y, 
                                       siteCovs = siteCovs, 
                                       obsCovs = obsCovs, 
                                       type = type)) 

startlynx <- qlogis(c(0.5, 0.1, 0.7)) 

( m1lynx <- occuFP(detformula = ~1,
               FPformula = ~1, 
               stateformula = ~1, 
               data = lynxumf) ) 

plogis(coef(m2)) #far away from simulated data: psi=0.6, p=0.7, fp=0.05!!!!


######
#
# Part 3B  Simulation study. What is the effect of false positives?
# Basic ideal : Just wrap the above code in a loop and fit a normal occupancy model.
#  (a methodological need: test the hypothesis that false positive rate is 0 -- this is on the boundary of the parameter space TBD)
#
######

set.seed(1) # Initialize RNGs
# PUT SEED OUTSIDE OF LOOP!

# matrix to hold the output
simout<- matrix(NA,nrow=100,ncol=2)
simout.fp<- matrix(NA,nrow=100,ncol=3)

for(simiter in 1:100){

nsites <- 200 # number of sites (i = 1, ..., nsites=M)
nsurveys <- 7 # number of visits (j = 1, ..., nsurveys=J)
psi <- 0.6    # expected occupancy probability
p <- 0.7      # detection probability (p_11)
fp <- 0.05    # false-positive error probability (p_10)

# Simulate occupancy states and encounter histories
z <- rbinom(nsites, 1, psi) # occupancy states
y <- matrix(NA, nrow = nsites, ncol = nsurveys) # empty matrix for detections
for(i in 1:nsites){
  pr_yequals1 <- p*z[i] + fp*(1 - z[i]) # p11 + p10
  y[i,] <- rbinom(nsurveys, 1, pr_yequals1) # realized observations
}

  
# Fit the model in unmarked. Data contaminated with false positives is type 2 
#  data according to the unmarked conventions. So here we have 7 occasions of T2 data:

type <- c(0, 7, 0)

# Build the unmarkedFrame
summary(umf <- unmarkedFrameOccuFP(y = y, type = type)) # not shown

# These starting values are consistent with our preferred constraint
largerp11 <- qlogis(c(0.5, 0.7, 0.1))  # Order is psi, p, fp   <- the order is important, otherwise you can't understand what's going on

 m2 <- occuFP(detformula = ~1, # model for p_11
               FPformula = ~1, # model for p_10
               stateformula = ~1, # model for psi
               data = umf, # umarkedFrameOccuFP object
               starts = largerp11)   # add p_11 < p_10 constraint
# Ignoring FPs

# Build the unmarkedFrame
 umf.occ <- unmarkedFrameOccu(y = y)  # 'type' not used
 
 m0 <- occu(~1 ~1, data = umf.occ)   # Now comma between formulas!

simout[simiter,]<- plogis(coef(m0)  )
simout.fp[simiter,]<- plogis(coef(m2))


}

# You can easily repeat this for fp = 0.01, 0.02, ..., 0.10 let's say



#  Part 4
#
#
# Bayesian analysis example
#
#
#

# Random seed and simulation settings
set.seed(129, kind = "Mersenne")
nsites <- 200 # number of sites (i = 1, ..., M)
nsurveys <- 7 # number of visits (k = 1, ..., J)
psi <- 0.6    # expected psi
p <- 0.7      # detection probability (p_11)
fp <- 0.05    # false positive error probability (p_10)

# Simulate the latent states and the data
z <- matrix(NA, nrow = nsites, ncol = 1) # empty matrix for occ states
z[1:nsites] <- rbinom(nsites, 1, psi) # occupancy states
y <- matrix(NA, nrow = nsites, ncol = nsurveys) # empty matrix for det.
for(i in 1:nsites){
  pr_yequals1 <- p*z[i] + fp*(1-z[i]) # p11 + p10
  y[i,] <- rbinom(nsurveys, 1, pr_yequals1) # realized observations
}

# Bundle data and summarize data bundle
str( bdata <- list(y = y, nsites = nrow(y), nsurveys = ncol(y)) )

# Specify model in BUGS language Royle Link model - type2
cat(file = "occufp.txt","
model {
# Priors
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  fp ~ dunif(0, 1)
# Likelihood and process model
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) # State model
    for (j in 1:nsurveys) { # Loop over replicate surveys
      y[i,j] ~ dbern(z[i]*p + (1-z[i])*fp) # Observation model
    }
  }
}
")


# Initial values
zst <- apply(y, 1, max)
inits <- function(){list(z = zst, p = 0.7, fp = 0.05)}

# Parameters monitored
params <- c("psi", "p", "fp")

# MCMC settings
na <- 1000 ; ni <- 5000 ; nt <- 1 ; nb <- 1000 ; nc <- 3

# Call JAGS (ART <1 min), assess convergence and summarize posteriors
library(jagsUI)
out1 <- jags(bdata, inits, params, "occufp.txt", n.adapt = na,
n.chains = nc,n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
 
print(out1, 3)
#very similar!


# BAD Initial values
inits <- function(){list(z = zst, p = 0.1, fp = 0.3)}

# Call JAGS (ART <1 min), assess convergence and summarize posteriors
out1b <- jags(bdata, inits, params, "occufp.txt", n.adapt = na,
n.chains = nc,n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
 
print(out1b)  # Illustrates the symmetry in the likelihood (and the wrong values are produced)


#####
#
# Class exercise break (see solution script)
#
#
#####

#
# Part 5
#
#
# Chambert et al. model for detection frequency data, contaminated with false positives
#
#

# Simulation settings
set.seed(2019, kind = "Mersenne")
nsites <- 100 # Number of sites
nsurveys <- 5 # Number of replicates/occasions
psi <- 0.7    # Occupancy
p11 <- 0.5    # Detection probability at an occupied site
p10 <- 0.05   # False detection probability
lam <- 3      # Rate of true positives from ARU
ome <- 0.50   # Rate of false positives from ARU
# Simulate true occupancy states
z <- rbinom(nsites, 1, psi)
# Define detection probability
p <- z * p11 + (1-z) * p10
# Simulate occupancy data and ARU count frequencies
yARU <- y <- K <- Q <- matrix(NA, nsites, nsurveys)
for(i in 1:nsites){
  y[i,] <- rbinom(nsurveys, 1, p[i]) # Detection/nondetection data
  K[i,] <- rpois(nsurveys, lam*z[i]) # True positive detection frequency
  Q[i,] <- rpois(nsurveys, ome)      # False-positive detection frequency
  yARU[i,] <- K[i,] + Q[i,]          # Number of ARU detections
}
# Bundle and summarize data
str( bdata <- list(y = y, yARU = yARU, nsites = nsites, nsurveys = nsurveys ))


# Specify Model A in BUGS language
cat(file = "modelA.txt","
model {
# Priors
psi ~ dunif(0, 1) # psi = Pr(Occupancy)
p10 ~ dunif(0, 1) # p10 = Pr(y = 1 | z = 0)
p11 ~ dunif(0, 1) # p11 = Pr(y = 1 | z = 1)
lam ~ dunif(0, 1000)
ome ~ dunif(0, 1000)
# Likelihood:process and observation models
for (i in 1:nsites) {
  z[i] ~ dbern(psi) # Occupancy status of site i
  p[i] <- z[i] * p11 + (1-z[i]) * p10 # false-positive detection model
  for(j in 1:nsurveys) {
    y[i,j] ~ dbern(p[i]) # Binary occupancy data
    yARU[i,j] ~ dpois(lam * z[i] + ome) # ARU detection frequency data
 }
}
}
")

# Initial values
inits <- function(){list(z = apply(y, 1, max), psi = runif(1),
p10 = runif(1, 0, 0.05), p11 = runif(1, 0.5, 0.8), lam = runif(1, 1, 2),
ome = runif(1, 0, 0.4) )}

# Parameters monitored
params <- c("psi", "p10", "p11", "lam", "ome")

# MCMC settings
na <- 1000 ; ni <- 2000 ; nt <- 1 ; nb <- 1000 ; nc <- 3

# Call JAGS (tiny ART), gauge convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "modelA.txt", n.adapt = na,
n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
 
print(out1, 3)























