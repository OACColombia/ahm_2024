#' ---
#' title: An R Markdown document converted from "royle-nichols.ipynb"
#' output: html_document
#' ---
#' 
#' # Royle-Nichols abundance model
#' 
#' ### Ken Kellner
#' 
#' # Outline
#' 
#' 1. Introduction
#' 2. Model description
#' 3. Dataset
#' 4. Fit model with NIMBLE
#' 6. Fit model with `unmarked`
#' 7. Exercise
#' 
#' # Introduction
#' 
#' Suppose you have two sites. Site 1 (in truth) has 3 individuals of a species, and site 2 has 10. All individuals have the same probability of detection.
#' 
#' You do a detection / non-detection survey.
#' 
#' Are you more likely to detect the species at one site or the other? Which one?
#' 
#' 
#' Overall $p$ will be higher at site 2.
#' 
#' This is called **abundance-induced heterogeneity in detection**.
#' 
#' ### Exploiting this relationship to estimate abundance
#' 
#' We can use the expected relationship between abundance and overall detection to estimate true abundance at a site
#' 
#' This requires only detection / non-detection data!
#' 
#' This is the famous "Royle-Nichols" model (Royle and Nichols 2003).
#' 
#' It will be a key ingredient in a model Josh will present later.
#' 
#' ### Key equation
#' 
#' Relating individual detection probability to overall detection probability
#' 
#' Parameter                                              | Parameter equation |
#' -------------------------------------------------------|--------------------|
#' Abundance a site                                       | $N$                |
#' Probability of detecting one specific animal           | $r$                |
#' Probability of **not** detecting that animal           | $(1-r)$            |
#' Probability of detecting **no** animals                | $(1-r)^N$          |
#' Probability of detecting **at least one** animal ($p$) | $1-(1-r)^N$        |
#' 
#' $$p = 1 - (1 - r)^N$$
#' 
#' ## Assumptions/Limitations
#' 
#' * Population is closed during sampling
#' * All individuals have equal probability of detection in each $j$
#' * Individuals are detected independently
#' * "Sampled area" is unknown - so estimating density is problematic
#' * Works best with small abundances (otherwise $p$ "saturates" at 1)
#' 
#' # Model description
#' 
#' ## State process
#' 
#' **Parameters**
#' 
#' $N_i$: Latent abundance at site $i$
#' 
#' $\lambda_i$: Expected abundance at site $i$
#' 
#' **Math**
#' 
#' $$N_i \sim \mathrm{Poisson}(\lambda_i)$$
#' 
#' (or some other distribution like negative binomial)
#' 
#' ## Detection process
#' 
#' **Parameters/data**
#' 
#' $y_{ij}$: Observed detection (1) / non-detection (0) at site $i$ for repeated sample $j$
#' 
#' $r_{ij}$: Probability of detecting an individual animal at site $i$ in sample $j$
#' 
#' $p_{ij}$: Probability of detecting *at least one* animal at site $i$ in sample $j$
#' 
#' 
#' **Math**
#' 
#' $$y_{ij} \sim \mathrm{Bernoulli}(p_{ij})$$
#' 
#' $$p_{ij} = 1 - (1 - r_{ij})^{N_i} $$
#' 
#' # Example dataset
#' 
#' ## Study species
#' 
#' Mallard (*Anas platyrhynchos*)
#' 
#' ![](mallard.jpg)
#' 
#' 
#' ## Study design
#' 
#' * National breeding bird monitoring program, Switzerland 
#' * M = 239 sites
#' * J = 3 replicate counts per site
#' * Kéry, Royle and Schmid (2005)
#' 
#' ## Look at the raw data
#' 
## -----------------------------------------------------------------------------
library(unmarked)
data(mallard)
head(mallard.y)

#' 
#' Data are counts, but we'll simplify them to detection/non-detection (will allow us to do a comparison later).
#' 
## -----------------------------------------------------------------------------
ybin <- mallard.y
ybin[ybin > 0] <- 1 # truncate

#' 
#' # Royle-Nichols model in NIMBLE
#' 
#' ## Simple model with no covariates
#' 
#' ### Bundle the data for NIMBLE
#' 
## -----------------------------------------------------------------------------
library(nimble)

#' 
## -----------------------------------------------------------------------------
nimble_data <- list(y = ybin, M = 239, J = 3)

#' 
#' ### Write the model in BUGS language
#' 
#' ### State model
#' 
#' **Math**
#' 
#' $$N_i \sim \mathrm{Poisson}(\lambda)$$
#' 
#' **Code**
#' 
#' ```r
#' for (i in 1:M){
#'   N[i] ~ dpois(lambda)   
#' }
#' lambda ~ dunif(0, 100) # prior should be positive
#' ```
#' 
#' ### Detection model
#' 
#' **Math**
#' 
#' $$p_i = 1 - (1 - r)^{N_i}$$
#' $$y_{ij} \sim \mathrm{Bernoulli}(p_i)$$
#' 
#' **Code**
#' 
#' ```r
#' for (i in 1:M){
#'   p[i] <- 1 - (1 - r)^N[i]
#'   for (j in 1:J){
#'     y[i,j] ~ dbern(p[i])
#'   }
#' }
#' r ~ dunif(0,1) # prior must be between 0 and 1
#' ```
#' 
#' ### Complete model
#' 
## -----------------------------------------------------------------------------
code_null <- nimbleCode({

# State model
for (i in 1:M){
  N[i] ~ dpois(lambda)   #latent abundance per sites 
}
lambda ~ dunif(0, 100) #should be bigger than 0

# Detection model
for (i in 1:M){
  p[i] <- 1 - (1 - r)^N[i]
  for (j in 1:J){
    y[i,j] ~ dbern(p[i])
  }
}
r ~ dunif(0,1) #prior of R

sumN <- sum(N[])

})

#' 
#' ### Run model
#' 
## -----------------------------------------------------------------------------
pars <- c("lambda", "r", "sumN")

# initialize z at maximum count ever observed at a site
Ninit <- apply(ybin, 1, max, na.rm=TRUE)
Ninit[is.infinite(Ninit)] <- NA
Ninit

inits <- function() list(N = Ninit)

#' 
## -----------------------------------------------------------------------------
nimble_null <- nimbleMCMC(code_null,
                          constants = nimble_data,
                          inits = inits, 
                          monitors = pars, 
                          nchains = 3, niter = 10000, nburnin = 8000, 
                          samplesAsCodaMCMC = TRUE)

#' 
## -----------------------------------------------------------------------------
par(mfrow=c(2,2))
coda::traceplot(nimble_null)

#' 
## -----------------------------------------------------------------------------
summary(nimble_null)

#' 
#' ## Model with covariates
#' 
#' ### Bundle the data for NIMBLE
#' 
#' Including the site and obs-level covariates
#' 
## -----------------------------------------------------------------------------
head(mallard.site)
lapply(mallard.obs, head)

#' 
#' Crudely impute missing values:
#' 
## -----------------------------------------------------------------------------
date_filled <- mallard.obs$date
date_filled[is.na(date_filled)] <- mean(mallard.obs$date, na.rm=TRUE)
ivel_filled <- mallard.obs$ivel
ivel_filled[is.na(ivel_filled)] <- mean(mallard.obs$ivel, na.rm=TRUE)

#' 
## -----------------------------------------------------------------------------
nimble_data <- list(y = ybin, M = 239, J = 3,
                    ivel = ivel_filled, date = date_filled,
                    elev = mallard.site$elev, forest = mallard.site$forest)

str(nimble_data)

#' 
#' ### Add covariates to BUGS code
#' 
#' ### Updated state model
#' 
#' **Math**
#' 
#' $$N_i \sim \mathrm{Poisson}(\lambda_i)$$
#' $$\mathrm{log}(\lambda_i) = \beta_0 + \beta_{ele} \cdot \mathrm{elev}_i + \beta_{for} \cdot \mathrm{forest}_i$$
#' 
#' **Code**
#' 
#' ```r
#' for (i in 1:M){
#'   N[i] ~ dpois(lambda[i])
#'   log(lambda[i]) <- beta_0 + beta_ele * elev[i] + beta_for * forest[i]
#' }
#' beta_0 ~ dunif(-10, 10)
#' beta_ele ~ dnorm(0, 0.01)
#' beta_for ~ dnorm(0, 0.01)
#' ```
#' 
#' ### Updated detection model
#' 
#' **Math**
#' 
#' $$p_{ij} = 1 - (1 - r_{ij})^{N_i}$$
#' $$y_{ij} \sim \mathrm{Bernoulli}(p_{ij})$$
#' 
#' **Code**
#' 
#' ```r
#' for (i in 1:M){
#'   for (j in 1:J){
#'     logit(r[i,j]) <- alpha_0 + alpha_date * date[i,j] + alpha_ivel * ivel[i,j]
#'     p[i,j] <- 1 - (1 - r[i,j])^N[i]
#'     y[i,j] ~ dbern(p[i])
#'   }
#' }
#' ```
#' 
## -----------------------------------------------------------------------------
code_covs <- nimbleCode({

# State model
for (i in 1:M){
  N[i] ~ dpois(lambda[i])
  log(lambda[i]) <- beta_0 + beta_ele * elev[i] + beta_for * forest[i]
}
  #Priors for estimate lambda_{i}
beta_0 ~ dunif(-5, 5)
beta_ele ~ dnorm(0, 0.1)
beta_for ~ dnorm(0, 0.1)

# Detection model
for (i in 1:M){
  for (j in 1:J){
    logit(r[i,j]) <- alpha_0 + alpha_date * date[i,j] + alpha_ivel * ivel[i,j]
    p[i,j] <- 1 - (1 - r[i,j])^N[i]
    y[i,j] ~ dbern(p[i,j])
  }
}
  #Priors for estimate r_{ij}
alpha_0 ~ dunif(-5, 5)
alpha_date ~ dnorm(0, 0.1)
alpha_ivel ~ dnorm(0, 0.1)
    
# Fit statistics (we will understand this after)
for (i in 1:M){
  z[i] <- step(N[i] - 1)
  for (j in 1:J){
    yhat[i,j] <- z[i] * p[i,j] + 0.001 # add small value to avoid divide by zero
    pres[i,j] <- (y[i,j] - yhat[i,j])^2 / yhat[i,j] * (1 - yhat[i,j])
    
    y_new[i,j] ~ dbern(z[i] * p[i,j])
    pres_new[i,j] <- (y_new[i,j] - yhat[i,j])^2 / yhat[i,j] * (1 - yhat[i,j])
    
  }
}
fit <- sum(pres[1:M,1:J])
fit_new <- sum(pres_new[1:M, 1:J])

sumN <- sum(N[])

})

#' 
#' ### Run model
#' 
## -----------------------------------------------------------------------------
pars <- c("beta_0", "beta_ele", "beta_for", "alpha_0", "alpha_date", "alpha_ivel", "sumN",
          "fit", "fit_new")
inits <- function() list(N = Ninit, beta_0=rnorm(1), alpha_0=rnorm(1))

set.seed(123)
nimble_covs <- nimbleMCMC(code_covs,
                          constants = nimble_data,
                          inits = inits, 
                          monitors = pars, 
                          nchains = 3, niter = 10000, nburnin = 8000, 
                          samplesAsCodaMCMC = TRUE)

#' 
## -----------------------------------------------------------------------------
par(mfrow=c(3,2))
coda::traceplot(nimble_covs)

#' 
## -----------------------------------------------------------------------------
coda::gelman.diag(nimble_covs)

#' 
## -----------------------------------------------------------------------------
summary(nimble_covs)

#' 
#' **Bayesian p-value**
#' 
## -----------------------------------------------------------------------------
post_covs <- as.matrix(nimble_covs)
mean(post_covs[,"fit"] > post_covs[,"fit_new"])

#' 
#' **Posterior predictive check plot**
#' 
#' Plot `fit` vs. `fit_new`.
#' 
#' A line with intercept 0 and slope 1 should approximately bisect the point cloud.
#' 
## -----------------------------------------------------------------------------
plot(post_covs[,"fit"], post_covs[,"fit_new"])
abline(a=0, b=1, col='red')

#' 
#' # Royle-Nichols model in `unmarked`
#' 
#' Look at the data again:
#' 
## -----------------------------------------------------------------------------
library(unmarked)
data(mallard)
head(mallard.y)
head(mallard.site)

#' 
#' ## Format data for unmarked
#' 
## -----------------------------------------------------------------------------
oc <- list(ivel = ivel_filled, date = date_filled)

umf <- unmarkedFrameOccu(y = ybin, siteCovs = mallard.site, obsCovs = oc)
head(umf)

#' 
#' ## Fit a null model
#' 
#' * Use the `unmarked` function `occuRN`
#' * Also should set a value for `K`, an upper bound on possible site abundance
#' 
## -----------------------------------------------------------------------------
rn_null <- occuRN(~1~1, umf, K = 20)
summary(rn_null)

#' 
#' Estimates are on the transformed scale! To get them on the original scale, use `predict`.
#' 
## -----------------------------------------------------------------------------
predict(rn_null, type='state')[1,]

#' 
## -----------------------------------------------------------------------------
predict(rn_null, type='det')[1,]

#' 
#' ## Exercise: fit model with covariates
#'
#' Add effect of elevation and forest to abundance, and ivel and date to detection.
#' 
#' Just need to change the formulas.
#' 



#' 
#' ## Model goodness-of-fit in unmarked
#' 
#' We will use the MacKenzie-Bailey chi-square test, in `AICcmodavg` package
#' 
## -----------------------------------------------------------------------------
library(AICcmodavg)
mb.gof.test(rn_covs, nsim = 30)

#' 
#' Compare to NIMBLE:
#' 
## -----------------------------------------------------------------------------
nim <- apply(post_covs[,c(4:6, 1:3)], 2, mean)
unm <- coef(rn_covs)

cbind(nimble=nim, unmarked=unm)

#' 
## -----------------------------------------------------------------------------
par(mfrow=c(2,2))
plotEffects(rn_covs, type = "state", covariate = "elev")
plotEffects(rn_covs, type = "state", covariate = "forest")
plotEffects(rn_covs, type = "det", covariate = "date")
plotEffects(rn_covs, type = "det", covariate = "ivel")

#' 
#' ## Compare R-N and N-mixture estimates
#' 
## -----------------------------------------------------------------------------
umf_nmix <- unmarkedFramePCount(y = mallard.y, siteCovs = mallard.site, obsCovs = oc)
nmix_covs <- pcount(~date+ivel~elev+forest, umf_nmix, K = 20)

#' 
#' ### Compare abundance model parameter estimates
#' 
## -----------------------------------------------------------------------------
cbind(RN = coef(rn_covs), Nmix = coef(nmix_covs))

#' 
#' ### Compare overall abundance estimates
#' 
## -----------------------------------------------------------------------------
cbind(RN = sum(bup(ranef(rn_covs))),
      Nmix = sum(bup(ranef(nmix_covs))))

#' 
