# You also need to install JAGS outside R, you can download it for
# your operating system here:
# https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/

library(jagsUI)
library(nimble)
library(coda)
library(unmarked)
library(ubms)
library(spOccupancy)
library(raster)

# Simulate an occupancy dataset
y <- matrix(NA, 50, 3)
psi <- plogis(0)
p <- plogis(0)
z <- rbinom(50, 1, psi)

for (i in 1:50){
  y[i,] <- rbinom(3, 1, z[i]*p)
}
umf <- unmarkedFrameOccu(y=y)

# Fit with unmarked
fit_unm <- occu(~1~1, umf)

# Fit with ubms
fit_ubm <- stan_occu(~1~1, umf, chains=3)

# Fit with spOccupancy
dat_sp <- list(y=umf@y, coords=matrix(runif(50*2, 0, 1), 50, 2))
fit_sp <- spPGOcc(~1, ~1, data=dat_sp, n.batch=10, batch.length=100,
                  n.chains=1, n.burn=50)

# Fit with JAGS
bugs_data <- list(M = 50, J = 3, y = umf@y)

modfile <- tempfile()
writeLines("
model{
  for (i in 1:M){
    z[i] ~ dbern(psi)
    for (j in 1:J){
      y[i,j] ~ dbern(z[i] * p)
    }
  }
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
}
", modfile)

inits <- function() list(z = rep(1, 50))
pars <- c("psi", "p")

fit_jags <- jags(bugs_data, inits, pars, modfile,
                 n.chains=2, n.iter=1000, n.burnin=500)

# Fit with NIMBLE
code <- nimbleCode({
  for (i in 1:M){
    z[i] ~ dbern(psi)
    for (j in 1:J){
      y[i,j] ~ dbern(z[i] * p)
    }
  }
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
})

fit_nimble <- nimbleMCMC(code, constants=bugs_data, inits=inits, monitors=pars,
                         nchains=2, niter=1000, nburnin=500, samplesAsCodaMCMC = TRUE)

coda::traceplot(fit_nimble)

message("Everything worked!")
