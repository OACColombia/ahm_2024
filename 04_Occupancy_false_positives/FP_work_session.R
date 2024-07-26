#Lynx example

#Lynx is invading (re-establishing). This is important, there are significant trends that have to
#be accounted for or else estimates are unstable.

#1. Create the unmarked data frame using unmarkedFrameOccuFP
#2. Do a likelihood analysis of this data set using occuFP to evaluate the magnitude of the false positive probability parameter
#3. Use country and sYear as covariates on 
#4. Model forest cover on occupancy. 

library(unmarked)
library(AHMbook)
data(EurasianLynx)
str(lynx <- EurasianLynx)

# Add the columns we need for analysis in unmarked
lynx$occ.1 <- 1
lynx$occ.2 <- 2
lynx$occ.3 <- 3
lynx$sYear <- standardize(lynx$Year)

# Extract the type 1 and type 2 data separately and bind them together
lynx1 <- lynx[lynx[,"type"] == "certain", ]
lynx2 <- lynx[lynx[,"type"] == "uncertain", ]
lynx <- cbind(lynx1[,c(2,3:5)], lynx2[,3:5] ) #seven colomns
colnames(lynx) <- c("site.nr", "y.1", "y.2", "y.3", "y.4", "y.5", "y.6")

occ <- cbind(lynx1[, c("occ.1", "occ.2", "occ.3")], lynx2[, c("occ.1", "occ.2", "occ.3")])
colnames(occ) <- c("occ.1", "occ.2", "occ.3", "occ.4", "occ.5", "occ.6")
lynx <- cbind(lynx,lynx1[, c("Year", "sYear", "Cntry","forest")])

# Make the false-positive unmarkedFrame. Be sure to indicate type!
y <- lynx[,paste0("y.", 1:6)]
siteCovs <- lynx[, c("sYear", "Year", "Cntry","forest")]
obsCovs <- list(occ = occ)

summary(lynx.umf <- unmarkedFrameOccuFP(y = y,
                                        siteCovs = siteCovs, 
                                        obsCovs = obsCovs, 
                                        type = c(3, 3, 0))) #number of occasions for: certain, uncertain, confirmation

# Note: There are 4 formulas -- the 3rd formula is the formula for the "b" parameter -- probability confirmed.
# You have to fill this slot , best to use NAMES for these formulas!

starts <- c(0,0,-3) #

# Null model
mod0<- occuFP(detformula = ~ 1, 
              FPformula = ~1, 
              Bformula = ~1, #no type3 data here
              stateformula = ~1, 
              data = lynx.umf,
              starts = starts)
plogis(coef(mod0))
#maybe the citizen science data (uncertain) is not affecting that much the false positive

# Trend in occupancy
mod1<- occuFP(detformula = ~ 1, FPformula = ~1, Bformula = ~1, stateformula = ~sYear, data = lynx.umf)

# Trend everywhere
mod2<- occuFP(detformula = ~ sYear, FPformula = ~sYear, Bformula = ~1, stateformula = ~sYear, data = lynx.umf)

# Add country
mod3<- occuFP(detformula = ~ Cntry + sYear, FPformula = ~Cntry + sYear, Bformula = ~1, stateformula = ~Cntry + sYear, data = lynx.umf)

# Forest
mod4<- occuFP(detformula = ~ Cntry + sYear, FPformula = ~Cntry + sYear, Bformula = ~1, stateformula = ~Cntry + sYear + forest, data = lynx.umf)

# Country is really not needed in FP rate, probably not in detection either
mod5<- occuFP(detformula = ~ sYear, FPformula = ~ sYear, Bformula = ~1, stateformula = ~Cntry + sYear + forest, data = lynx.umf)

mods<- list("null" = mod0, "~1~1~Trend" = mod1,
   "~Trend~Trend~Trend" = mod2,
   "~Cntry+Trend~Cntry+Trend~Cntry+Trend"=mod3,
   "~Cntry+Trend~Cntry+Trend~Cntry+Trend+forest"=mod4,
   "~Trend~Trend~Cntry+Trend+forest" = mod5)
modSelFP(mods)

# A potential BUG here
modSel(fitList(mods))
 


######
#
#  Now let's do a Bayesian analysis of the Lynx data
#
######


# Bundle data and summarize data bundle
str( bdata <- list(y = as.matrix(y), nsites = nrow(y), sYear = siteCovs$sYear, forest = siteCovs$forest) )

# Specify model in BUGS language
cat(file = "occufp.txt","
model {
# Priors
psi ~ dunif(0, 1)
p ~ dunif(0, 1)
fp ~ dunif(0, 1)
# Likelihood and process model
for (i in 1:nsites) { # Loop over sites
  z[i] ~ dbern(psi) # State model
  for(j in 1:3){
    y[i,j] ~ dbern(z[i]*p) # Observation model Type 1 data

  }
  for (j in 4:6) { # Loop over replicate surveys
    y[i,j] ~ dbern(z[i]*p + (1-z[i])*fp) # Observation model Type 2 data
  }
 }
}
")


# Initial values
zst <- apply(y, 1, max,na.rm=TRUE)
inits <- function(){list(z = zst, p = 0.7, fp = 0.02)}

# Parameters monitored
params <- c("psi", "p", "fp")

# MCMC settings
na <- 400 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 6

# Call JAGS (ART <1 min), assess convergence and summarize posteriors
library(jagsUI)
lynxout <- jags(bdata, inits, params, "occufp.txt", n.adapt = na,
n.chains = nc,n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
  
#
# Now build in sYear and forest
#

# Specify model in BUGS language
cat(file = "occufp2.txt","
model {
# Priors
psi ~ dunif(0, 1)
p ~ dunif(0, 1)
fp ~ dunif(0, 1)

# Likelihood and process model
for (i in 1:nsites) { # Loop over sites
  z[i] ~ dbern(psi ) # State model
   
  for(j in 1:3){
    y[i,j] ~ dbern(z[i]*p) # Observation model Type 1 data

  }
  for (j in 4:6) { # Loop over replicate surveys
    y[i,j] ~ dbern(z[i]*p + (1-z[i])*fp) # Observation model Type 2 data
  }
 }
}
")


# Initial values
zst <- apply(y, 1, max,na.rm=TRUE)
inits <- function(){list(z = zst, p = 0.7, fp = 0.02)}

# Parameters monitored
params <- c("psi", "p", "fp")

# MCMC settings
na <- 400 ; ni <- 1000 ; nt <- 1 ; nb <- 500 ; nc <- 6

# Call JAGS (ART <1 min), assess convergence and summarize posteriors
library(jagsUI)
lynxout <- jags(bdata, inits, params, "occufp.txt", n.adapt = na,
n.chains = nc,n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
 
print(out1, 3)
