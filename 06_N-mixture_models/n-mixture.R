#' ---
#' title: An R Markdown document converted from "n-mixture.ipynb"
#' output: html_document
#' ---
#' 
#' # Static N-mixture models
#' 
#' ### Ken Kellner
#' 
#' # Outline
#' 
#' 1. Background
#' 2. Model description
#' 3. Dataset
#' 4. Fit model with NIMBLE
#' 5. Exercise
#' 6. Fit model with `unmarked`
#' 7. Exercise
#' 
#' # Background
#' 
#' How to estimate abundance without marking animals?
#' 
#' ### Repeated samples!
#' 
#' In an occupancy model:
#' * We record repeated detection/non-detection samples (0/1)
#' * We model the unobserved occupancy state at a site (0/1)
#' 
#' Replace detection with counts and occupancy with abundance:
#' * We record repeated *count* samples (0, 1, 2, ...)
#' * We model the unobserved *abundance* at a site (0, 1, 2, ...)
#' 
#' ### N-mixture model
#' 
#' This is the N-mixture model of Royle (2004)
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
#' You could also use other distributions, such as:
#' 
#' * Negative binomial
#' * Zero-inflated Poisson
#' 
#' ## Detection process
#' 
#' **Parameters/data**
#' 
#' $y_{ij}$: Observed count at site $i$ for repeated sample $j$
#' 
#' $p_{ij}$: Probability of detecting an individual at site $i$ in sample $j$
#' 
#' **Math**
#' 
#' $$y_{ij} \sim \mathrm{Binomial}(p_{ij}, N_i)$$
#' 
#' ## Assumptions/Limitations
#' 
#' * Population is closed during sampling
#' * All individuals have equal probability of detection in each $j$
#' * Individuals are detected independently
#' * "Sampled area" is unknown - so estimating density is problematic
#' * Works best with small-medium counts; large counts can be very hard to fit
#' * See AHM Vol 1, pg 248
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
#' # N-mixture model in NIMBLE
#' 
#' ## Simple model with no covariates
#' 
#' ### Bundle the data for NIMBLE
#' 
## -----------------------------------------------------------------------------
library(nimble)

#' 
## -----------------------------------------------------------------------------
nimble_data <- list(y = mallard.y, M = 239, J = 3)

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
#' $$y_{ij} \sim \mathrm{Binomial}(p,N_i)$$
#' 
#' **Code**
#' 
#' ```r
#' for (i in 1:M){
#'   for (j in 1:J){
#'     y[i,j] ~ dbinom(p, N[i])
#'   }
#' }
#' p ~ dunif(0,1) # prior must be between 0 and 1
#' ```
#' 
#' ### Complete model
#' 
## -----------------------------------------------------------------------------
code_null <- nimbleCode({

# State model
for (i in 1:M){
  N[i] ~ dpois(lambda)   
}
lambda ~ dunif(0, 100)

# Detection model
for (i in 1:M){
  for (j in 1:J){
    y[i,j] ~ dbinom(p, N[i])
  }
}
p ~ dunif(0,1) #prior between 0-1

sumN <- sum(N[])

})

#' 
#' ### Run model
#' 
## -----------------------------------------------------------------------------
pars <- c("lambda", "p", "sumN")

# initialize z at maximum count ever observed at a site
Ninit <- apply(mallard.y, 1, max, na.rm=TRUE) + 1
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
#' Including the site-level covariates
#' 
## -----------------------------------------------------------------------------
head(mallard.site)

# already standardized?
round(colMeans(mallard.site), 2)
round(apply(mallard.site, 2, sd), 2)

#' 
## -----------------------------------------------------------------------------
nimble_data <- list(y = mallard.y, M = 239, J = 3,
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
## -----------------------------------------------------------------------------
code_covs <- nimbleCode({

# State model
for (i in 1:M){
  N[i] ~ dpois(lambda[i])
  log(lambda[i]) <- beta_0 + beta_ele * elev[i] + beta_for * forest[i]
}
beta_0 ~ dunif(-10, 10)
beta_ele ~ dnorm(0, 0.01)
beta_for ~ dnorm(0, 0.01)

# Detection model
for (i in 1:M){
  for (j in 1:J){
    y[i,j] ~ dbinom(p, N[i])
  }
}
p ~ dunif(0,1)

sumN <- sum(N[])

})

#' 
#' ### Run model
#' 
## -----------------------------------------------------------------------------
pars <- c("beta_0", "beta_ele", "beta_for", #these are associated with lambda
          "p", "sumN")
inits <- function() list(N = Ninit, beta_0=rnorm(1), alpha_0=rnorm(1))

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
#' ## Add goodness-of-fit test
#' 
#' **In each MCMC iteration**:
#' 
#' 1. Calculate a fit statistic `fit` (Pearson's $\chi^2$) on the real data
#' 2. Simulate a new dataset based on the fitted model
#' 3. Calculate the fit statistic `fit_new` for the new dataset
#' 
#' **After the run**:
#' 
#' * Calculate the proportion of iterations where `fit < fit_new` = Bayesian p-value
#' * If model fits data well, should be ~ 0.5
#' 
#' **Math**
#' 
#' $$\hat{y}_{ij} = N_i \cdot p $$
#' $$\mathrm{fit}_{\chi^2} = \sum{\frac{(y_{ij} - \hat{y}_{ij})^2}{\hat{y}_{ij}}}$$
#' 
#' **Code**
#' 
#' ```r
#' # real data
#' for (i in 1:M){
#'   for (j in 1:J){
#'     yhat[i,j] = N[i] * p
#'     chi2[i,j] = (y[i,j] - yhat[i,j])^2 / yhat[i,j]
#'   }
#' }
#' fit = sum(chi2[1:M, 1:J])
#' ```
#' 
#' **Simulate new dataset**
#' 
#' $$y^{new}_{ij} \sim \mathrm{Binomial}(p, N_i)$$
#' 
#' **Calculate fit for new dataset**
#' 
#' $$\hat{y}_{ij} = N_i \cdot p $$
#' $$\mathrm{fit}^{new}_{\chi^2} = \sum{\frac{(y^{new}_{ij} - \hat{y}_{ij})^2}{\hat{y}_{ij}}}$$
#' 
#' **Code**
#' 
#' ```r
#' # simulated data
#' for (i in 1:M){
#'   for (j in 1:J){
#'     y_new[i,j] ~ dbinom(p, N[i]) # simulate new datapoint
#'     yhat[i,j] = N[i] * p
#'     chi2_new[i,j] = (y_new[i,j] - yhat[i,j])^2 / yhat[i,j]
#'   }
#' }
#' fit_new = sum(chi2_new[1:M, 1:J]) #indexed in nimble!
#' ```
#' 
## -----------------------------------------------------------------------------
code_covs <- nimbleCode({

# State model
for (i in 1:M){
  N[i] ~ dpois(lambda[i])
  log(lambda[i]) <- beta_0 + beta_ele * elev[i] + beta_for * forest[i]
}
beta_0 ~ dunif(-10, 10)
beta_ele ~ dnorm(0, 0.01)
beta_for ~ dnorm(0, 0.01)

# Detection model
for (i in 1:M){
  for (j in 1:J){
    y[i,j] ~ dbinom(p, N[i])
  }
}
p ~ dbeta(1,1)

# Fit statistic for real data
for (i in 1:M){
  for (j in 1:J){
    yhat[i,j] <- N[i] * p + 0.001 # add small value to avoid divide by zero
    chi2[i,j] <- (y[i,j] - yhat[i,j])^2 / yhat[i,j]
  }
}
fit <- sum(chi2[1:M,1:J])

# Fit statistic for simulated data
for (i in 1:M){
  for (j in 1:J){
    y_new[i,j] ~ dbinom(p, N[i]) # simulate new datapoint
    chi2_new[i,j] <- (y_new[i,j] - yhat[i,j])^2 / yhat[i,j]
  }
}
fit_new <- sum(chi2_new[1:M,1:J])

sumN <- sum(N[])

})

#' 
## -----------------------------------------------------------------------------
pars <- c("beta_0", "beta_ele", "beta_for", "p", "fit", "fit_new", "sumN")
inits <- function() list(N = Ninit, beta_0=rnorm(1), alpha_0=rnorm(1))

nimble_covs <- nimbleMCMC(code_covs,
                          constants = nimble_data,
                          inits = inits, 
                          monitors = pars, 
                          nchains = 3, niter = 10000, nburnin = 8000, 
                          samplesAsCodaMCMC = TRUE)

#' 
## -----------------------------------------------------------------------------
par(mfrow=c(3,3))
coda::traceplot(nimble_covs)

#' 
#' **Bayesian p-value**
#' 
## -----------------------------------------------------------------------------
post <- as.matrix(nimble_covs)
mean(post[,"fit"] > post[,"fit_new"])

#' 
#' **Posterior predictive check plot**
#' 
#' Plot `fit` vs. `fit_new`.
#' 
#' A line with intercept 0 and slope 1 should approximately bisect the point cloud.
#' 
## -----------------------------------------------------------------------------
plot(post[,"fit"], post[,"fit_new"])
abline(a=0, b=1)

#' 
#' Based on this test, the model appears to fit the data poorly.
#' 
#' Potential solutions:
#' 
#' * Has the model converged? (maybe not in this case)
#' * Are key covariates missing from the model?
#' * Other random effects needed?
#' * Change assumed distribution for $N$ (e.g. negative binomial, zero-inflated Poisson)
#' * Check to make sure your priors are not unduly influencing the results (if you are using "uninformative" priors)
#' 
#' More research should be done in this area to develop better GoF tests.
#' 
#' # Exercise
#' 
#' Use NIMBLE to fit the following model.
#' 
#' * Abundance: elevation + forest + elevation<sup>2</sup>
#' * Detection: ivel (see dataset `mallard.obs`)
#' Line 505
#' 
#' Note there are missing values in the observation covariate `ivel`! See the dynamic occupancy model code for one way to handle this in NIMBLE. Hint: you can model the missing values of `ivel`.
#' 
#' Don't forget: `p` now becomes indexed by `i` and `j`
#' 
#' Check goodness-of-fit.
#' 
## -----------------------------------------------------------------------------

head(mallard.site)
str(mallard.obs)


#' 
#' ## Other distributions for N
#' 
#' * Frequently with count data, a Poisson distribution may not be reasonable because mean != variance.
#' * Specifically, variance > mean (overdispersion)
#' 
#' Solution: add another parameter
#' 
#' * Negative binomial
#' * Zero-inflation
#' * Poisson-lognormal (obs-level random effect)
#' 
#' ### Observation-level random effect
#' 
#' (also called Poisson-lognormal model)
#' 
#' **Math**
#' 
#' $$N_i \sim \mathrm{Poisson}(\lambda_i)$$
#' $$\mathrm{log}(\lambda_i) = \beta_0 + \beta_{ele} \cdot \mathrm{elev}_i + \beta_{for} \cdot \mathrm{forest}_i + \epsilon_i$$
#' $$\epsilon_i \sim \mathrm{Normal}(0, \sigma^2)$$
#' 
#' **Code**
#' 
#' ```r
#' for (i in 1:M){
#'   N[i] ~ dpois(lambda[i])
#'   log(lambda[i]) <- beta_0 + beta_ele * elev[i] + beta_for * forest[i] +
#'                     beta_ele2 * elev[i] * elev[i] + eps[i]
#'   eps[i] ~ dnorm(0, tau)
#' }
#' beta_0 ~ dunif(-10, 10)
#' beta_ele ~ dnorm(0, 0.01)
#' beta_for ~ dnorm(0, 0.01)
#' beta_ele2 ~ dnorm(0, 0.01)
#' tau <- pow(sigma, -2)
#' sigma ~ dunif(0, 10)
#' ```
#' 
## -----------------------------------------------------------------------------
code_pln <- nimbleCode({

# State model
for (i in 1:M){
  N[i] ~ dpois(lambda[i])
  log(lambda[i]) <- beta_0 + beta_ele * elev[i] + beta_for * forest[i] +
                    beta_ele2 * elev[i] * elev[i] + eps[i]
  eps[i] ~ dnorm(0, tau) # random effect
}
beta_0 ~ dunif(-10, 10)
beta_ele ~ dnorm(0, 0.01)
beta_for ~ dnorm(0, 0.01)
beta_ele2 ~ dnorm(0, 0.01)
tau <- pow(sigma, -2)
sigma ~ dunif(0, 10)

# Model missing values in ivel
for (i in 1:M){
  for (j in 1:J){
   ivel[i,j] ~ dnorm(ivel_mean, ivel_tau)
  }
}
ivel_mean ~ dnorm(0, 0.01)
ivel_tau <- pow(ivel_sd, -2)
ivel_sd ~ dunif(0, 10)

# Detection model
for (i in 1:M){
  for (j in 1:J){
    logit(p[i,j]) <- alpha_0 + alpha_ivel * ivel[i,j]
    y[i,j] ~ dbinom(p[i,j], N[i])
  }
}
alpha_0 ~ dunif(-10,10)
alpha_ivel ~ dnorm(0, 0.01)

# Fit statistic for real data
for (i in 1:M){
  for (j in 1:J){
    yhat[i,j] <- N[i] * p[i,j] + 0.001 # add small value to avoid divide by zero
    chi2[i,j] <- (y[i,j] - yhat[i,j])^2 / yhat[i,j]
  }
}
fit <- sum(chi2[1:M, 1:J])

# Fit statistic for simulated data
for (i in 1:M){
  for (j in 1:J){
    y_new[i,j] ~ dbinom(p[i,j], N[i]) # simulate new datapoint
    chi2_new[i,j] <- (y_new[i,j] - yhat[i,j])^2 / yhat[i,j]
  }
}
fit_new <- sum(chi2_new[1:M, 1:J])

sumN <- sum(N[])

})

#' 
## -----------------------------------------------------------------------------
pars <- c("beta_0", "beta_ele", "beta_for", "beta_ele2", "alpha_0", "alpha_ivel", "sigma", "fit", "fit_new", "sumN")
inits <- function() list(N = Ninit, beta_0=rnorm(1), alpha_0 = rnorm(1), sigma=rlnorm(1))

nimble_pln <- nimbleMCMC(code_pln,
                         constants = nimble_data,
                         inits = inits, 
                         monitors = pars, 
                         nchains = 3, niter = 10000, nburnin = 8000, 
                         samplesAsCodaMCMC = TRUE) 

#' 
## -----------------------------------------------------------------------------
par(mfrow=c(3,3))
coda::traceplot(nimble_pln)

#' 
## -----------------------------------------------------------------------------
summary(nimble_pln)

#' 
## -----------------------------------------------------------------------------
post <- as.matrix(nimble_pln)
mean(post[,"fit"]>post[,"fit_new"])

#' 
## -----------------------------------------------------------------------------
plot(post[,"fit"], post[,"fit_new"])
abline(a=0, b=1, col='red')

#' 
#' ### Comparison effect of elevation
#' 
## -----------------------------------------------------------------------------
post_covs <- as.matrix(nimble_covs)
post_pln <- as.matrix(nimble_pln)

est <- c(mean(post_covs[,"beta_ele"]), mean(post_pln[,"beta_ele"]))
lower <- c(quantile(post_covs[,"beta_ele"], 0.025), quantile(post_pln[,"beta_ele"], 0.025))
upper <- c(quantile(post_covs[,"beta_ele"], 0.975), quantile(post_pln[,"beta_ele"], 0.975))


plot(1:2, est, pch=19, xaxt='n', xlim=c(0.5, 2.5),
     ylim=c(-2.5, -0.5), xlab="", ylab="beta_ele")
axis(1, at=1:2, labels=c("Pois","PLN"))
segments(1:2, lower, 1:2, upper)

#' 
#' ### Comparison total N
#' 
## -----------------------------------------------------------------------------
est <- c(mean(post_covs[,"sumN"]), mean(post_pln[,"sumN"]))
lower <- c(quantile(post_covs[,"sumN"], 0.025), quantile(post_pln[,"sumN"], 0.025))
upper <- c(quantile(post_covs[,"sumN"], 0.975), quantile(post_pln[,"sumN"], 0.975))

plot(1:2, est, pch=19, xaxt='n', xlim=c(0.5, 2.5),
     ylim=c(70, 190), xlab="", ylab="sumN")
axis(1, at=1:2, labels=c("Pois","PLN"))
segments(1:2, lower, 1:2, upper)

#' 
#' # N-mixture model in `unmarked`
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
umf <- unmarkedFramePCount(y = mallard.y, siteCovs = mallard.site)
head(umf)

#' 
#' ## Fit a null model
#' 
#' * Use the `unmarked` function `pcount`
#' * `pcount` fits the N-mixture model: does not need to be point-count data!
#' * Also should set a value for `K`, an upper bound on possible site abundance
#' 
## -----------------------------------------------------------------------------
nmix_null <- pcount(~1~1, umf, K = 20) #K is integrated
summary(nmix_null)

#' 
#' Estimates are on the transformed scale! To get them on the original scale, use `predict`.
#' 
## -----------------------------------------------------------------------------
predict(nmix_null, type='state')[1,]

#' 
## -----------------------------------------------------------------------------
predict(nmix_null, type='det')[1,]

#' 
#' ### Impact of K value setting
#' 
## -----------------------------------------------------------------------------
max(umf@y, na.rm=TRUE)

#' 
## -----------------------------------------------------------------------------
test_K <- function(K){
  fit <- suppressWarnings(pcount(~1~1, umf, K=K))
  exp(coef(fit)[1])
}
kvals <- c(13,14,15,20,25,30) #starting with value > max(umf@y, na.rm=TRUE)
lams <- sapply(kvals, test_K)
plot(kvals, lams, type='o', pch=19, xlim=c(10, 30),
     xlab="K", ylab="lambda estimate")
abline(v=12, col='red', lty=2)

#' 
#' ## Fit model with covariates
#' 
#' Just need to change the formulas
#' 
## -----------------------------------------------------------------------------
nmix_covs <- pcount(~1~elev+forest, umf, K = 20)
nmix_covs

#' 
#' The estimates are a little different from NIMBLE, but estimates in both cases are quite uncertain
#' 
## -----------------------------------------------------------------------------
nim <- apply(post_covs[,c(1:3,6)], 2, mean)
unm <- c(coef(nmix_covs)[1:3], plogis(coef(nmix_covs)[4]))

cbind(nimble=nim, unmarked=unm)

#' 
#' ## Plot covariate effects
#' 
#' 1. Make sequence of covariate values from min to max
#' 2. Insert into `newdata`, holding other covariates at mean or median
#' 2. Use `predict` to generate corresponding abundance values
#' 3. Plot
#' 
## -----------------------------------------------------------------------------
# sequence of covariate values
elev_rng <- range(mallard.site$elev)                 #this is going to vary
elev_seq <- seq(elev_rng[1], elev_rng[2], length.out=100)
# insert into newdata
nd <- data.frame(elev=elev_seq,                      #this is going to vary
                 forest=median(mallard.site$forest)) #pick value arbitrariely, fixed
# predict
pr <- predict(nmix_covs, type='state', newdata=nd)
head(pr)

#' 
## -----------------------------------------------------------------------------
plot(elev_seq, pr$Predicted, type='l')

#' 
## -----------------------------------------------------------------------------
# add CIs
plot(elev_seq, pr$Predicted, type='l', ylim=c(0, 1.5))
polygon(c(elev_seq, rev(elev_seq)), c(pr$lower, rev(pr$upper)), col='lightgrey', border=NA)
lines(elev_seq, pr$Predicted)

#' 
## -----------------------------------------------------------------------------
plotEffects(nmix_covs, type = "state", covariate = "elev")

#Empirical Bayes value 




#' 
#' ## Model goodness-of-fit in unmarked
#' 
#' We will use tools in `AICcmodavg` package
#' 
#' **Process:**
#' 
#' 1. Calculate fit metric (Pearson's $\chi^2$) using model and real dataset
#' 2. Simulate $n$ new datasets from fitted model
#' 3. Calculate fit metric for each simulated dataset
#' 4. Compare distribution of metrics from simulated datasets with real dataset
#' 
#' $$\mathrm{Pearson's \chi^2} = \sum{\frac{(y_{ij} - \hat{y}_{ij})^2}{\hat{y}_{ij}}}$$
#' 
#' $$\hat{y}_{ij} = \lambda_i \cdot p $$
#' 
#' Calculate manually in R:
#' 
## -----------------------------------------------------------------------------
yhat <- fitted(nmix_covs)
sum((umf@y - yhat)^2 / yhat, na.rm=TRUE)

#' 
## -----------------------------------------------------------------------------
library(AICcmodavg)
Nmix.gof.test(nmix_covs, nsim = 30)

#' 
#### Try negative binomial ####
#' 
## -----------------------------------------------------------------------------
nmix_covs_nb <- pcount(~1~elev+forest, umf, K = 20, mixture = "NB")
summary(nmix_covs_nb)

#' 
## -----------------------------------------------------------------------------
Nmix.gof.test(nmix_covs_nb, nsim = 30) #with nsim = 1000

sum(bup(ranef(nmix_covs_nb)))
#104 individuals

#' 
# ### Zero-inflated Poisson ####
#' 
## -----------------------------------------------------------------------------
nmix_covs_zip <- pcount(~1~elev+forest, umf, K = 20, mixture = "ZIP")
summary(nmix_covs_zip)

#' 
## -----------------------------------------------------------------------------
Nmix.gof.test(nmix_covs_zip, nsim = 30)

#' 
#' ### PLN
#' 
## -----------------------------------------------------------------------------
umf@siteCovs$id <- factor(1:numSites(umf)) #Poisson log-normal require an id vlaue
nmix_covs_pln <- pcount(~1~elev+forest+(1|id), #random effect
                        umf, K = 20) 
nmix_covs_pln

#' 
## -----------------------------------------------------------------------------
Nmix.gof.test(nmix_covs_pln, nsim = 30)

#' 
#' ### Estimate of total N
#' 
## -----------------------------------------------------------------------------
sum(bup(ranef(nmix_covs_pln)))

#' 
#' ## Model selection
#' 
## -----------------------------------------------------------------------------
fl <- fitList(null=nmix_null, 
              covs=nmix_covs, 
              nb=nmix_covs_nb, 
              zip=nmix_covs_zip)
modSel(fl)

#' 
#' # Exercise
#' 
#' Use unmarked to fit the following model.
#' 
#' * Abundance: elevation + forest + elevation<sup>2</sup>
#' * Detection: ivel (see dataset `mallard.obs`, you will need to modify the `unmarkedFrame`)
#' 
#' Check goodness-of-fit.
#' 
#' `I(elev^2)`
#' 
## -----------------------------------------------------------------------------

umf_obs <- unmarkedFramePCount(y = mallard.y,
                               siteCovs = mallard.site, 
                               obsCovs = mallard.obs)

nmix_obs <- pcount(~ivel~elev+forest+I(elev^2),
                   umf_obs)

nmix_obs

#' 
#' ## Fit the PLN model with `ubms`
#' 
## -----------------------------------------------------------------------------
library(ubms)
nmix_stan <- stan_pcount(~1~elev+forest+(1|id), umf, K = 20, chains = 3, iter = 300, refresh = 0)
nmix_stan

#' 
#' ### Get the posterior of latent N
#' 
## -----------------------------------------------------------------------------
Nest <- posterior_predict(nmix_stan, "z") #get posterior distribution of N
dim(Nest)

#' 
## -----------------------------------------------------------------------------
hist(apply(Nest, 1, sum, na.rm=TRUE), main="sumN", xlab="sumN")

#' 


# ### Prior sensitivity analysis ####
#' 
#' Demonstration of how prior settings can affect coefficient estimates
#' 
#' Set the `prior_intercept_state` argument in `stan_pcount`
#' 
## -----------------------------------------------------------------------------
vary_prior <- function(prior){
    suppressWarnings(stan_pcount(~1~elev+forest, umf, K = 20, chains = 3, iter = 300, refresh = 0,
                prior_intercept_state = prior))
}

priors <- list(p5=normal(0, 5), p1=normal(0, 1), p25=normal(0, 2.5), p05=normal(0, 0.5), p01=normal(0, 0.1))

mods <- lapply(priors, vary_prior)

#' 
## -----------------------------------------------------------------------------
stats <- t(sapply(mods, function(x) as.numeric(summary(x, 'state')[1,c(1,4,8)])))
                  
plot(1:5, stats[,1], ylim=c(-2.5, 0), pch=19, xaxt='n',
     xlab="", ylab="Intercept estimate")
axis(1, at=1:5, labels=c("N(0,5)","N(0,2.5)","N(0,1)", "N(0,0.5)","N(0,0.1)"))
segments(1:5, stats[,2], 1:5, stats[,3])

#' 
