#' ---
#' title: An R Markdown document converted from "dynamic-occupancy.ipynb"
#' output: html_document
#' ---
#' 
#' # Dynamic occupancy models
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
#' Site occupancy can change over time (depending on your definition of time!)
#' 
#' Examples:
#' * Anthropogenic land use change affects species distribution
#' * Spread of an invasive species
#' * Reintroduction of endangered species
#' * Climate change impacts reduces/increases species range
#' 
#' We want to understand what drives these patterns of change, while accounting for imperfect detection.
#' 
#' **Solution:**
#' 
#' "Dynamic" or multi-season occupancy model (MacKenzie et al. 2003 *Ecology*)
#' 
#' **Model overview diagram**
#' 
#' **Review Concepts**
#' 
#' * Sites have initial ocupancy
#' * Occupancy at later states depend on previous state and some transition probability (Markov)
#' 
#' * Assumptions same as single-season model
#' * Most importantly closure within each primary period
#' * Don't necessarily need balance data
#' 
#' **Review parameters of interest**
#' 
#' State submodel:
#' 
#' * Initial probability of occupancy: $\psi$
#' * Probability occupied sites become unoccupied (local *extinction*): $\epsilon$
#' * Probability occupied sites *stay* occupied (persistance, 1 - $\epsilon$): $\phi$
#' * Probability unoccupied sites become occupied (*colonization*): $\gamma$
#' 
#' Detection submodel:
#' 
#' * Probability of detection: $p$
#' 
#' # Model description
#' 
#' ### State model at $t=1$
#' 
#' **Parameters**:
#' 
#' $z_{i, 1}$: Initial occupancy state of site $i$
#' 
#' $\psi$: Initial occupancy probability at time $t=1$
#' 
#' **Math**:
#' 
#' $$
#' z_{i,1} \sim \mathrm{Bernoulli}(\psi)
#' $$
#' 
#' ### State model for $t > 1$
#' 
#' **Parameters**:
#' 
#' $z_{i,t}$: Occupancy state of site $i$ in time $t$
#' 
#' $\epsilon$: Probability of extinction, i.e. probability site occupied in time $t-1$ is unoccupied in time $t$
#' 
#' $\phi$: Probability of persistence; $1-\epsilon$
#' 
#' $\gamma$: Probability of colonization; i.e. probability non-occupied site becomes occupied in time $t$
#' 
#' **Math**:
#' 
#' $$
#' z_{i,t} \sim \begin{cases}
#' \mathrm{Bernoulli}(1 - \epsilon), & \mathrm{if} z_{i, t-1} = 1 \\
#' \mathrm{Bernoulli}(\gamma), & \mathrm{if} z_{i, t-1} = 0
#' \end{cases}
#' $$
#' 
#' **Alternatively**:
#' 
#' $$
#' z_{i,t} \sim \mathrm{Bernoulli}(z_{i,t-1}\cdot(1-\epsilon) + (1-z_{i,t-1})\cdot\gamma)
#' $$
#' 
#' When $z_{i,t-1} = 1$:
#' 
#' $$
#' z_{i,t} \sim \mathrm{Bernoulli}(1\cdot(1-\epsilon) + 0\cdot\gamma) = 
#' $$
#' 
#' $$z_{i,t} \sim \mathrm{Bernoulli}(1-\epsilon)$$
#' 
#' When $z_{i,t-1} = 0$:
#' 
#' $$
#' z_{i,t} \sim \mathrm{Bernoulli}(0\cdot(1-\epsilon) + 1\cdot\gamma) = 
#' $$
#' 
#' $$z_{i,t} \sim \mathrm{Bernoulli}(\gamma)$$
#' 
#' ### Detection submodel
#' 
#' **Parameters**:
#' 
#' $p$: Probability of detection
#' 
#' $z_{i,t}$: Occupancy state of site $i$ in time $t$
#' 
#' **Data**:
#' 
#' $y_{i,j,t}$: Detection/non-detection data at site $i$, time $t$, occasion $j$
#' 
#' **Math**:
#' 
#' $$
#' y_{i,j,t} \sim \mathrm{Bernoulli}(p \cdot z_{i,t})
#' $$
#' 
#' When $z_{i,t} = 1$:
#' 
#' $$y_{i,j,t} \sim \mathrm{Bernoulli}(p\cdot1)$$
#' 
#' When $z_{i,t} = 0$:
#' 
#' $$y_{i,j,t} \sim \mathrm{Bernoulli}(p\cdot0) = 0$$
#' 
#' ### Data needed
#' 
#' * Detection/non-detection data from multiple primary periods ($t$)
#' * As with single-season occupancy, need repeated samples ($j$) at at least some sites in some primary periods
#' 
#' ### To learn more:
#' 
#' * MacKenzie, D. I., Nichols, J. D., Hines, J. E., Knutson, M. G., & Franklin, A. B. (2003). Estimating site occupancy, colonization, and local extinction when a species is detected imperfectly. Ecology, 84(8), 2200-2207.
#' 
#' * Royle, J. A., & Kéry, M. (2007). A Bayesian state‐space formulation of dynamic occupancy models. Ecology, 88(7), 1813-1823.
#' 
#' * AHM vol 2, Chapter 4
#' 
#' # Example dataset
#' 
#' *Loxia curvirostra*, European (common) crossbill
#' 
#' ![](Loxia_curvirostra.jpg)
#' 
#' ## Study design
#' 
#' * 267 quadrats
#' * 3 surveys per year
#' * 9 years
#' 
#' ## Look at the raw data
#' 
#' It comes with `unmarked`.
#' 
## -----------------------------------------------------------------------------
library(unmarked)
data(crossbill)
head(crossbill)

#' 
#' # Dynamic occupancy model in NIMBLE
#' 
#' ## Simple model with no covariates
#' 
#' ### Format the y-matrix for NIMBLE
#' 
#' First, select only the columns corresponding to detection/non-detection data
#' 
## -----------------------------------------------------------------------------
ycols <- grepl("det", names(crossbill)) #in base R!!! grepl(pattern, x)
head(ycols)
y <- as.matrix(crossbill[,ycols])
head(y)

#' 
#' The y-matrix is in "wide" format, we want to convert it to a M x J x T array
#' 
#' M: Sites
#' 
#' J: Sampling occasions
#' 
#' T: Time periods
#' 
## -----------------------------------------------------------------------------
M <- nrow(y)
J <- 3
T <- ncol(y) / J

#' 
## -----------------------------------------------------------------------------
y_array <- array(NA, dim = c(M, J, T))
period <- rep(1:T, each = J)
period

#' 
## -----------------------------------------------------------------------------
#iterate to fill the y_array
for (t in 1:T){
    y_array[1:M,1:J,t] <- y[, which(period == t)]
}

print(y_array[1:5,,1:3])
#1-3 years

#' 
#' ### Bundle the data for NIMBLE
#' 
## -----------------------------------------------------------------------------
library(nimble)

#' 
## -----------------------------------------------------------------------------
nimble_data <- list(M = M, J = J, T = T,
                  y = y_array
                 )
str(nimble_data)

#' 
#' ### Write the model in BUGS language
#' 
#' **Initial occupancy state**
#' 
#' Math:
#' 
#' $$
#' z_{i,1} \sim \mathrm{Bernoulli}(\psi)
#' $$
#' 
#' BUGS:
#' 
#' ```r
#' for (i in 1:M){
#'     z[i, 1] ~ dbern(psi)
#' }
#' 
#' (notice z is a matrix! and the "1" is the first period)
#' ```
#' 
#' **Occupancy state at time t**
#' 
#' Math:
#' 
#' at $t = 1$:
#' 
#' $$
#' z_{i,1} \sim \mathrm{Bernoulli}(\psi)
#' $$
#' 
#' at $t > 1$:
#' 
#' 
#' $$ z_{i,t} \sim \mathrm{Bernoulli}(z_{i, t-1}\cdot(1-\epsilon) + (1-z_{i, t-1})\cdot\gamma) $$
#' 
#' BUGS:
#' 
#' ```r
#' for (i in 1:M){
#'     
#'     z[i, 1] ~ dbern(psi)
#'     
#'     for (t in 2:T){
#'         z[i, t] ~ dbern(z[i, t-1] * (1 - eps) + (1 - z[i, t-1]) * gam)
#'     }
#' }
#' ```
#' 
#' **Detection process**
#' 
#' Math:
#' 
#' $$ y_{i, j, t} \sim \mathrm{Bernoulli}(p \cdot z_{i, t}) $$
#' 
#' BUGS:
#' 
#' ```r
#' # Detection - three dimensions in y_array
#' for (i in 1:M){
#'     for (t in 1:T){
#'         for (j in 1:J){
#'             y[i,j,t] ~ dbern(p * z[i, t])
#'         }
#'     }
#' }
#' ```
#' 
#' ### Complete model code
#' 
## -----------------------------------------------------------------------------
code_null <- nimbleCode({

# State
for (i in 1:M){
    
    z[i, 1] ~ dbern(psi)

    for (t in 2:T){
        z[i, t] ~ dbern(z[i, t-1] * (1 - eps) + (1 - z[i, t-1]) * gam)
    }
}

# Detection
for (i in 1:M){
    for (t in 1:T){
        for (j in 1:J){
            y[i,j,t] ~ dbern(p * z[i, t])
        }
    }
}

# Priors
psi ~ dunif(0, 1)
eps ~ dunif(0, 1)
gam ~ dunif(0, 1)
p ~ dunif(0, 1)

})

#' 
#' ### Other inputs
#' 
## -----------------------------------------------------------------------------
pars <- c("psi", "eps", "gam", "p")

#' 
#' As with single-season occupancy, need to initialize the latent state reasonably.
#' 
## -----------------------------------------------------------------------------
inits <- function() list(z = matrix(1, M, T))

#' 
#' ### Run model
#' 
## -----------------------------------------------------------------------------
nimble_null <- nimbleMCMC(code = code_null,
                  constants = nimble_data, 
                  inits = inits, 
                  monitors = pars, 
                  nchains = 3, niter = 3000, nburnin = 1000,
                  samplesAsCodaMCMC = TRUE)

#' 
## -----------------------------------------------------------------------------
par(mfrow=c(2,2))
coda::traceplot(nimble_null)
coda::gelman.diag(nimble_null)

#' 
## -----------------------------------------------------------------------------
summary(nimble_null)

#' 
#' ### Plot the posterior of `psi`
#' 
## -----------------------------------------------------------------------------
mat <- as.matrix(nimble_null)
dim(mat)
head(mat)
hist(mat[,"psi"])
abline(v=mean(mat[,"psi"]), col='red')

#' 
#' ### Plot 95% credible intervals
#' 
## -----------------------------------------------------------------------------
psi_est <- mean(mat[,"psi"])
psi_ci <- quantile(mat[,"psi"], c(0.025, 0.975))

plot(1, psi_est, cex=2, pch=19, xaxt='n', xlab="psi")
segments(1, psi_ci[1], 1, psi_ci[2])

#' 
#' ## Model with covariates
#' 
#' ### Site covariates
#' 
## -----------------------------------------------------------------------------
head(crossbill)

#' 
#' ### Format the detection covariate matrix for JAGS
#' 
#' Select only the columns corresponding to date data
#' 
#' Select text based on some pattern:
#' 
#' ```
#' grepl
#' ```
#' 
#' Incredibly powerful tool for finding things in R, definitely worth learning to use (along with regular expressions)
#' 
## -----------------------------------------------------------------------------
datecols <- grepl("date", names(crossbill))
date <- as.matrix(crossbill[,datecols])
head(date)

#' 
#' Convert to an array, the same way we did with the y-matrix
#' 
## -----------------------------------------------------------------------------
date_array <- array(NA, dim = c(M, J, T))
period <- rep(1:T, each = J)

for (t in 1:T){
    date_array[,,t] <- date[, which(period == t)]
}
print(date_array[1:5,,1:3])

#' 
#' ### Bundle the data for NIMBLE
#' 
#' Including the site-level covariates, which we standardize
#' 
## -----------------------------------------------------------------------------
nimble_data <- list(M = M, J = J, T = T,
                  y = y_array,
                  date = (date_array - mean(date_array, na.rm=TRUE)) / sd(date_array, na.rm=TRUE), #standardized
                  ele = (crossbill$ele - mean(crossbill$ele)) / sd(crossbill$ele),
                  forest = (crossbill$forest - mean(crossbill$forest)) / sd(crossbill$forest)
                 )
str(nimble_data)

#' 
#' ### Add covariates to BUGS code
#' 
#' Initial occupancy: Elevation
#' 
#' Colonization: Forest
#' 
#' Extinction: None
#' 
#' Detection: Date
#' 
#' **Initial occupancy**:
#' 
#' $$ \mathrm{logit}(\psi_{i,1}) = \alpha_0 + \alpha_{ele} \cdot \mathrm{ele}_i $$
#' 
#' ```r
#' for (i in 1:M){
#'   logit(psi[i,1]) <- alpha_0 + alpha_ele * ele[i]
#' }
#' # psi priors
#' alpha_0 ~ dlogis(0, 1)
#' alpha_ele ~ dlogis(0, 1)
#' ```
#' 
#' **Why `dlogis`**?
#' 
#' Looks better on the original (probability) scale.
#' 
## -----------------------------------------------------------------------------
# Typical normal prior
hist(plogis(rnorm(1000, 0, 10))) #like a smile, a little bad

#' 
## -----------------------------------------------------------------------------
# dlogis prior
hist(plogis(rlogis(1000, 0, 1))) #a little more uniform

#' 
#' **Add linear predictor for `gamma`**:
#' 
#' $$ \mathrm{logit}(\gamma_{i,t}) = \beta_0 + \beta_{forest} \cdot \mathrm{forest}_{i} $$
#' 
#' ```r
#' for (i in 1:M){
#'   logit(psi[i,1]) <- alpha_0 + alpha_ele * ele[i]
#'   for (t in 2:T){
#'     logit(gam[i,t]) <- beta_0 + beta_forest * forest[i]
#'   }
#' }
#' 
#' # psi priors
#' alpha_0 ~ dlogis(0, 1)
#' alpha_ele ~ dlogis(0, 1)
#' # gamma priors
#' beta_0 ~ dlogis(0, 1)
#' beta_ele ~ dlogis(0, 1)
#' # eps priors
#' eps ~ dunif(0, 1)
#' ```
#' 
#' 
#' * Extinction (`eps`) has no covariates, so it stays the same as the old model
#' * Note that `t` starts with 2
#' * We could also start `t` at 1 and go to `T-1` instead
#' 
#' **Detection parameters**:
#' 
#' $$ \mathrm{logit}(p_{i,j,t}) = \rho_0 + \rho_{date} \cdot \mathrm{date}_{i,j,t} $$
#' 
#' * Now we need to also loop over the $J$ occasions in each of the $T$ years
#' * Note here that $t$ starts with 1 instead of 2.
#' 
#' ```r
#' for (i in 1:M){
#'   for (t in 1:T){
#'     for (j in 1:J){
#'       logit(p[i,j,t]) <- rho_0 + rho_date * date[i,j,t]
#'     }
#'   }  
#' }
#' 
#' # detection priors
#' rho_0 ~ dlogis(0, 1)
#' rho_date ~ dlogis(0, 1)
#' ```
#' 
#' ### Complete BUGS model
#' 
## -----------------------------------------------------------------------------
code_covs <- nimbleCode({

# State submodel
for (i in 1:M){
    
    z[i, 1] ~ dbern(psi[i,1])

    for (t in 2:T){
        z[i, t] ~ dbern(z[i, t-1] * (1 - eps) + (1 - z[i, t-1]) * gam[i,t])
    }
}

# ----------------------

# Calculate psi, gamma, eps
for (i in 1:M){

    logit(psi[i,1]) <- alpha_0 + alpha_ele * ele[i]

    for (t in 2:T){
        logit(gam[i,t]) <- beta_0 + beta_forest * forest[i]
    }
}

# psi priors
alpha_0 ~ dlogis(0, 1)
alpha_ele ~ dlogis(0, 1)
# gamma priors
beta_0 ~ dlogis(0, 1)
beta_forest ~ dlogis(0, 1)
# eps priors
eps ~ dunif(0, 1)

#--------------------------

# Detection submodel
for (i in 1:M){
    for (t in 1:T){
        for (j in 1:J){
            y[i,j,t] ~ dbern(p[i,j,t] * z[i, t])
        }
    }
}

# -----------------------

# Calculate p
for (i in 1:M){
    for (t in 1:T){
        for (j in 1:J){
            logit(p[i,j,t]) <- rho_0 + rho_date * date[i,j,t]
        }
    }  
}
rho_0 ~ dlogis(0, 1)
rho_date ~ dlogis(0, 1)

})

#' 
#' ### Run model
#' 
## -----------------------------------------------------------------------------
pars <- c("alpha_0", "alpha_ele", "beta_0", "beta_forest",
          "rho_0", "rho_date")

inits <- function() list(z = matrix(1, M, T))

nimble_covs <- nimbleMCMC(code = code_covs,
                          constants = nimble_data, 
                          inits = inits, 
                          monitors = pars, 
                          nchains = 3, niter = 1000, nburnin = 500, 
                          samplesAsCodaMCMC=TRUE)

#' 
## -----------------------------------------------------------------------------

par(mfrow=c(4,2))
coda::traceplot(nimble_covs)

#' 
#' Any guesses about the problem?
#' 
#' ### Missing values
#' 
#' (sometimes the result of unequal sample sizes)
#' 
#' A very common problem in multi-year datasets!
#' 
## -----------------------------------------------------------------------------
nimble_data$date[43,3,1]

#' 
## -----------------------------------------------------------------------------
sum(is.na(nimble_data$date))

#' 
#' Solutions:
#' 
#' * Change the indexing to work around the missing values
#' * Impute the missing values
#' 
#' ### Impute missing dates
#' 
#' Missing date values are simulated as coming from a normal distribution with mean and SD of the original dataset
#' 
#' ```r
#' # Impute missing dates
#' for (i in 1:M){
#'     for (t in 1:T){
#'         for (j in 1:J){
#'             date[i,j,t] ~ dnorm(mu_date, tau_date)
#'         }
#'     }  
#' }
#' mu_date ~ dnorm(0, 0.01)
#' sigma_date ~ dunif(0, 10)
#' tau_date <- pow(sigma_date, -2) #this is precision = 1/sigmasqr
#' ```
#' 
#' It's OK that `date` has a mixture of data and missing values!
#' 
#' One way to think about it:
#' 
#' * When `date[i,j,k]` is not `NA`: information flows from left to right
#' * When `date[i,j,k]` is `NA`: information flows from right to left
#' 
#' ### Complete model .... again
#' 
## -----------------------------------------------------------------------------
code_covs <- nimbleCode({

# State submodel
for (i in 1:M){
    
    z[i, 1] ~ dbern(psi[i,1])

    for (t in 2:T){
        z[i, t] ~ dbern(z[i, t-1] * (1 - eps) + (1 - z[i, t-1]) * gam[i,t])
    }
}

# ----------------------

# Calculate psi, gamma, eps
for (i in 1:M){

    logit(psi[i,1]) <- alpha_0 + alpha_ele * ele[i]

    for (t in 2:T){
        logit(gam[i,t]) <- beta_0 + beta_forest * forest[i]
    }
}

# psi priors
alpha_0 ~ dlogis(0, 1)
alpha_ele ~ dlogis(0, 1)
# gamma priors
beta_0 ~ dlogis(0, 1)
beta_forest ~ dlogis(0, 1)
# eps priors
eps ~ dunif(0, 1)

# -----------------------

# Detection submodel
for (i in 1:M){
    for (t in 1:T){
        for (j in 1:J){
            y[i,j,t] ~ dbern(p[i,j,t] * z[i, t])
        }
    }
}

# -----------------------

# Impute missing dates
for (i in 1:M){
    for (t in 1:T){
        for (j in 1:J){
            date[i,j,t] ~ dnorm(mu_date, tau_date)
        }
    }  
}
mu_date ~ dnorm(0, 0.01)
sigma_date ~ dunif(0, 10)
tau_date <- pow(sigma_date, -2)

# -----------------------

# Calculate p
for (i in 1:M){
    for (t in 1:T){
        for (j in 1:J){
            logit(p[i,j,t]) <- rho_0 + rho_date * date[i,j,t]
        }
    }  
}
rho_0 ~ dlogis(0, 1)
rho_date ~ dlogis(0, 1)

})

#' 
## -----------------------------------------------------------------------------
pars <- c("alpha_0", "alpha_ele", "beta_0", "beta_forest",
          "rho_0", "rho_date", "mu_date", "sigma_date")

inits <- function() list(z = matrix(1, M, T))

#' 
## -----------------------------------------------------------------------------
nimble_covs <- nimbleMCMC(code = code_covs,
                          constants = nimble_data, 
                          inits = inits, 
                          monitors = pars, 
                          nchains = 3, niter = 3000, nburnin = 1000, 
                          samplesAsCodaMCMC=TRUE)

#' 
## -----------------------------------------------------------------------------
#saveRDS(nimble_covs, "nimble_covs.Rds")
#nimble_covs <- readRDS('nimble_covs.Rds')

#' 
## -----------------------------------------------------------------------------
par(mfrow=c(4,2))
coda::traceplot(nimble_covs)
coda::gelman.diag(nimble_covs)

#' 
## -----------------------------------------------------------------------------
options(scipen=5)
summary(nimble_covs, digits=2) #see mu_date close to 0, sigma_date close to 1!

#' 
#' ### Plot effect of elevation on `psi`
#' 
#' 1. Generate sequence of elevation values
#' 2. Standardize them (the original values, no the sequence!)
#' 3. Plug them into the linear predictor for `psi`
#' 4. Plot results
#' 
## -----------------------------------------------------------------------------
elev_rng <- range(crossbill$ele)
elev_seq <- seq(elev_rng[1], elev_rng[2], length.out=100) # 1.
elev_std <- (elev_seq - mean(crossbill$ele)) / sd(crossbill$ele) # 2.
head(elev_std)

#' 
#' $$ \mathrm{logit}(\psi_{i,1}) = \alpha_0 + \alpha_{ele} \cdot \mathrm{ele}_i $$
#' 
#' ```r
#' logit(psi[i,1]) <- alpha_0 + alpha_ele * ele[i]
#' ```
#' 
## -----------------------------------------------------------------------------
post <- as.matrix(nimble_covs)
alpha_0 <- mean(post[,"alpha_0"])
alpha_ele <- mean(post[,"alpha_ele"])
psi_est <- alpha_0 + alpha_ele * elev_std # 3.(Vectorized)
psi_est <- plogis(psi_est)
head(psi_est)

#' 
## -----------------------------------------------------------------------------
plot(elev_seq, psi_est, type='l', xlab="Elevation", ylab="psi") # 4.

#' 
#' **What about the 95% credible interval?**
#' 
#' We need to calculate `psi` for *each saved iteration* of parameters in the NIMBLE output (not just the posterior mean).
#' 
## -----------------------------------------------------------------------------
nsamples <- nrow(mat)
nsamples

#' 
## -----------------------------------------------------------------------------
psi_post <- matrix(NA, nrow=nsamples, ncol=100)
for (i in 1:100){
  psi_post[,i] <- post[,"alpha_0"] + post[,"alpha_ele"] * elev_std[i]
}
psi_post <- plogis(psi_post)
psi_post[1:5,1:5]

#' 
## -----------------------------------------------------------------------------
psi_lower <- apply(psi_post, 2, quantile, 0.025) # apply quantile over the columns (2) lower bound of 95% ci
psi_upper <- apply(psi_post, 2, quantile, 0.975) # upper bound of 95% ci

head(psi_lower)
head(psi_upper)

#' 
#' **Plot in base R**
#' 
## -----------------------------------------------------------------------------
plot(elev_seq, psi_est, type='l', xlab="Elevation", ylab="psi", ylim=c(0,0.8))
polygon(c(elev_seq, rev(elev_seq)),
        c(psi_lower, rev(psi_upper)), 
        col='grey90', border=NA)
lines(elev_seq, psi_est)

#' 
#' **Plot in ggplot2**
#' 
## -----------------------------------------------------------------------------
library(ggplot2)

plot_data <- data.frame(elev=elev_seq, psi=psi_est, lower=psi_lower, upper=psi_upper)
head(plot_data)

#' 
## -----------------------------------------------------------------------------
ggplot(data=plot_data, aes(x=elev)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  geom_line(aes(y=psi)) +
  theme_bw(base_size=18) +
  labs(y = expression(psi))+
  theme(panel.grid=element_blank())

#' 
#' ### Estimate occupancy probability for years 2...T
#' 
#' Right now we only estimate $\psi$ at $t=1$!
#' 
#' What about for other times?
#' 
#' * For this we need a *derived* parameter, that is a function of other parameters already in the model
#' 
#' * $\psi_{t}$ at times $t>1$ is a function of occupancy in the previous time step, extinction, and colonization
#' 
#' $$\psi_{i,t} = \psi_{i,t-1} \cdot (1-\epsilon_{i,t}) + (1 - \psi_{i,t}) * \gamma_{i,t} $$
#' 
#' $$\psi_{i,t} = \psi_{i,t-1} \cdot (1-\epsilon_{i,t}) + (1 - \psi_{i,t}) * \gamma_{i,t} $$
#' 
#' This math is very similar to BUGS code we already have in our model:
#' 
#' ```r
#' for (t in 2:T){
#'     z[i, t] ~ dbern(z[i, t-1] * (1 - eps) + (1 - z[i, t-1]) * gam[i,t])
#' }
#' ```
#' 
#' We just need to replace `z` with `psi`.
#' 
#' ```r
#' # Derive psi at t>1
#' for (i in 1:M){
#'     for (t in 2:T){
#'         psi[i, t] <- psi[i, t-1] * (1 - eps) + (1 - psi[i, t-1]) * gam[i,t]
#'     }
#' }
#' ```
#' 
#' * Note that `psi` is now indexed by both `i` and `t`
#' 
#' ### Model goodness-of-fit
#' 
#' 1. Calculate some test statistic on data
#' 2. Simulate new data based on the model
#' 3. Calculate the test statistic on simulated data
#' 4. Compare
#' 
#' One possible statistic is **sum of squared Pearson residuals (SSE)**.
#' 
#' $$ r = \frac{obs - expect}{\sqrt{(Var(expect))}} $$
#' 
#' $$ r^2 = \frac{(obs - expect)^2}{Var(expect)} $$
#' 
#' obs: $y_{ijt}$
#' 
#' expect: $\hat{y}_{ijt} = \psi_{it} \cdot p_{ijt}$
#' 
#' Var(expect): $\hat{y}_{ijt} \cdot (1 - \hat{y}_{ijt})$
#' 
#' $$ \sum r^2 =  \sum{\frac{(y_{ijt} - \hat{y}_{ij})^2}{\hat{y}_{ijt}  \cdot (1 - \hat{y}_{ijt})}} $$
#' 
#' $$\hat{y}_{ijt} = \psi_{it} \cdot p_{ijt} $$
#' 
#' 
#' Now that we have `psi` for all years we can calculate this.
#' 
#' ```r
#' # Goodness-of-fit
#' for (i in 1:M){
#'     for (t in 1:T){
#'         for (j in 1:J){
#'             yhat[i,j,t] <- psi[i,t] * p[i,j,t]
#'             # Pearson residual
#'             pres[i,j,t] <- (y[i,j,t] - yhat[i,j,t])^2 / (yhat[i,j,t] * (1-yhat[i,j,t]))
#'             # Simulate new data
#'             y_new[i,j,t] ~ dbern(p[i,j,t] * z[i, t])
#'             # Pearson residual for new data
#'             pres_new[i,j,t] <- (y_new[i,j,t] - yhat[i,j,t])^2 / (yhat[i,j,t] * (1-yhat[i,j,t]))
#'         }
#'     }  
#' }
#' sse <- sum(pres[1:M, 1:J, 1:T]) #require the dimmensions of the array
#' sse_new <- sum(pres_new[1:M, 1:J, 1:T])
#' ```
#' 
## -----------------------------------------------------------------------------
code_covs <- nimbleCode({

# State submodel
for (i in 1:M){
    
    z[i, 1] ~ dbern(psi[i,1])

    for (t in 2:T){
        z[i, t] ~ dbern(z[i, t-1] * (1 - eps) + (1 - z[i, t-1]) * gam[i,t])
    }
}

# ----------------------

# Calculate psi, gamma, eps
for (i in 1:M){

    logit(psi[i,1]) <- alpha_0 + alpha_ele * ele[i]

    for (t in 2:T){
        logit(gam[i,t]) <- beta_0 + beta_forest * forest[i]
    }
}

# psi priors
alpha_0 ~ dlogis(0, 1)
alpha_ele ~ dlogis(0, 1)
# gamma priors
beta_0 ~ dlogis(0, 1)
beta_forest ~ dlogis(0, 1)
# eps priors
eps ~ dunif(0, 1)

# -----------------------

# Detection submodel
for (i in 1:M){
    for (t in 1:T){
        for (j in 1:J){
            y[i,j,t] ~ dbern(p[i,j,t] * z[i, t])
        }
    }
}
# -----------------------

# Impute missing dates
for (i in 1:M){
    for (t in 1:T){
        for (j in 1:J){
            date[i,j,t] ~ dnorm(mu_date, tau_date)
        }
    }  
}
mu_date ~ dnorm(0, 0.01)
sigma_date ~ dunif(0, 10)
tau_date <- pow(sigma_date, -2)

# -----------------------

# Calculate p
for (i in 1:M){
    for (t in 1:T){
        for (j in 1:J){
            logit(p[i,j,t]) <- rho_0 + rho_date * date[i,j,t]
        }
    }  
}
rho_0 ~ dlogis(0, 1)
rho_date ~ dlogis(0, 1)

# -----------------------

# Derive psi at t>1
for (i in 1:M){
    for (t in 2:T){
        psi[i, t] <- psi[i, t-1] * (1 - eps) + (1 - psi[i, t-1]) * gam[i,t]
    }
}
    
# Goodness-of-fit
for (i in 1:M){
    for (t in 1:T){
        for (j in 1:J){
            yhat[i,j,t] <- psi[i,t] * p[i,j,t]
            pres[i,j,t] <- (y[i,j,t] - yhat[i,j,t])^2 / (yhat[i,j,t] * (1-yhat[i,j,t]))
            y_new[i,j,t] ~ dbern(p[i,j,t] * z[i, t])
            pres_new[i,j,t] <- (y_new[i,j,t] - yhat[i,j,t])^2 / (yhat[i,j,t] * (1-yhat[i,j,t]))
        }
    }  
}

# Sum up the residuals
# nimble can only sum over max 2 dimensions at a time
for (i in 1:M){
    sse_temp[i] <- sum(pres[i, 1:J, 1:T])
    sse_new_temp[i] <- sum(pres_new[i, 1:J, 1:T])
}
sse <- sum(sse_temp[1:M])
sse_new <- sum(sse_new_temp[1:M])

})

#' 
## -----------------------------------------------------------------------------
pars <- c("alpha_0", "alpha_ele", "beta_0", "beta_forest",
          "rho_0", "rho_date", "mu_date", "sigma_date", "psi",
          "sse", "sse_new")

inits <- function() list(z = matrix(1, M, T))

#' 
## -----------------------------------------------------------------------------
nimble_covs <- nimbleMCMC(code = code_covs,
                          constants = nimble_data, 
                          inits = inits, 
                          monitors = pars, 
                          nchains = 3, niter = 3000, nburnin = 1000, 
                          samplesAsCodaMCMC=TRUE)

#' 
## -----------------------------------------------------------------------------
#saveRDS(nimble_covs, "nimble_covs_gof.Rds")
#nimble_covs <- readRDS("nimble_covs_gof.Rds")

#' 
#' ### Posterior predictive check
#' 
#' * Plot posterior of SSE for real data vs. SSE for simulated data
#' * If model fits well, SSE for simulated data should be bigger about half the time, and vice-versa
#' * The ratio is the Bayesian p-value and should be about 0.5
#' 
## -----------------------------------------------------------------------------
post <- as.matrix(nimble_covs)
mean(post[,"sse"] > post[,"sse_new"])

#' 
## -----------------------------------------------------------------------------
plot(post[,"sse"], post[,"sse_new"])
abline(a=0, b=1, col='red')

#' 
#' Fit is not so good.
#' 
#' ### Plot `psi` over time for sites 1-5
#' 
#' `psi` in the NIMBLE output is in "wide" format
#' 
## -----------------------------------------------------------------------------
head(post)

#' 
#' We want `psi` in a matrix format (sites by time) so we can loop over it
#' 
#' ```r
#' psi_est <- psi_lower <- psi_upper <- matrix(NA, nrow=5, ncol=9)
#' for (t in 1:9){
#'   for (i in 1:5){
#'     psi_it <- ????
#'     psi_est[i,t] <- mean(psi_it)
#'     psi_lower[i,t] <- quantile(psi_it, 0.025)
#'     psi_upper[i,t] <- quantile(psi_it, 0.975)
#'   }
#' }
#' ```
#' 
#' How to reorganize the `psi` columns, i.e. convert between wide and "matrix" format?
#' 
## -----------------------------------------------------------------------------
paste0("psi[",1,", ", 1,"]")

#' 
## -----------------------------------------------------------------------------
psi_est <- psi_lower <- psi_upper <- matrix(NA, nrow=5, ncol=9)
for (t in 1:9){
  for (i in 1:5){
    psi_name <- paste0("psi[",i,", ",t,"]")
    psi_it <- post[,psi_name]
    psi_est[i,t] <- mean(psi_it)
    psi_lower[i,t] <- quantile(psi_it, 0.025)
    psi_upper[i,t] <- quantile(psi_it, 0.975)
  }
}

#' 
## -----------------------------------------------------------------------------
cols <- rainbow(5)
plot(1:9, psi_est[1,1:9], ylim=c(0,0.5), type='o', pch=19, 
     col = cols[1],
     ylab="psi", xlab="time")
for (i in 2:5){
  lines(1:9, psi_est[i,1:9], type='o', pch=19, col=cols[i])    
}

#' 
#' ### Exercise
#' 
#' Use NIMBLE (or JAGS if you want) to fit the following model:
#' 
#' * initial occupancy: elevation + elevation^2 + forest
#' * colonization: forest
#' * extinction: elevation
#' * detection: date
#' 
#' Also try plotting the effect of elevation on extinction.
#' 
#' Probably remove the parts calculating `psi` and `sse`, to speed things up.
#' 
#' **NIMBLE model**
#' 
## -----------------------------------------------------------------------------




#' 
#' # Dynamic occupancy model in `unmarked`
#' 
#' Look at the data again:
#' 
## -----------------------------------------------------------------------------
head(crossbill)

#' 
#' Grab the detection/non-detection data.
#' 
#' This time we want to keep it in "wide" format (M x TJ), because that's what `unmarked` uses.
#' 
## -----------------------------------------------------------------------------
ycols <- grepl("det", names(crossbill))
y <- crossbill[,ycols]
names(y)
head(y)

#' 
#' Get the site covariates:
#' 
## -----------------------------------------------------------------------------
site_covs <- crossbill[,1:4]
head(site_covs)

#' 
#' And the observation covariate (date):
#' 
## -----------------------------------------------------------------------------
datecols <- grepl("date", names(crossbill))
date <- crossbill[,datecols]
names(date)

#' 
#' **Make the `unmarkedFrame`**
#' 
#' Special unmarked frame type: `unmarkedMultFrame`
#' 
#' We need to specify the number of seasons (`numPrimary`)
#' 
## -----------------------------------------------------------------------------
umf <- unmarkedMultFrame(y = y, 
                         siteCovs = site_covs, 
                         obsCovs = list(date = date), 
                         numPrimary = 9)
summary(umf)

#' 
#' ### Fit a null model
#' 
#' We will use the `colext` function.
#' 
#' We need to specify R formulas for four parameters:
#' 
#' * $\psi$, `psiformula`
#' * $\epsilon$, `epsilonformula`
#' * $\gamma$, `gammaformula`
#' * $p$, `pformula`
#' 
#' For the null model, all four will be intercept-only, or `~1`.
#' 
## -----------------------------------------------------------------------------
mod_null <- colext(psiformula = ~1, 
                   epsilonformula = ~1, 
                   gammaformula = ~1, 
                   pformula = ~1, 
                   data = umf)
summary(mod_null)

#' 
#' ### Get parameters on the probability scale
#' 
## -----------------------------------------------------------------------------
plogis(-0.795) #for psi

#' 
#' Or use the `backTransform` function:
#' 
## -----------------------------------------------------------------------------
backTransform(mod_null, 'psi')

#' 
## -----------------------------------------------------------------------------
backTransform(mod_null, 'det')

#' 
## -----------------------------------------------------------------------------
confint(backTransform(mod_null, 'psi'))

#' 
#' Or use the `predict` function:
#' 
## -----------------------------------------------------------------------------
head(predict(mod_null, 'psi'))

#' 
#' ### Comparison with NIMBLE estimates
#' 
## -----------------------------------------------------------------------------
cbind(NIMBLE=apply(as.matrix(nimble_null), 2, mean)[c(4,2,1,3)],
      unmarked=plogis(coef(mod_null)))

#' 
#' ## Model with covariates
#' 
#' * $\psi$: elevation
#' * $\epsilon$: intercept-only
#' * $\gamma$: forest
#' * $p$: date
#' 
#' (same as with NIMBLE)
#' 
## -----------------------------------------------------------------------------
mod_covs <- colext(psiformula = ~scale(ele), 
                   epsilonformula = ~1, 
                   gammaformula = ~scale(forest), 
                   pformula = ~scale(date), 
                   data = umf)
summary(mod_covs)

#' 
#' **Get estimates on probability scale**
#' 
## -----------------------------------------------------------------------------
head(predict(mod_covs, "psi"))

#' 
## -----------------------------------------------------------------------------
head(predict(mod_covs, "det"))

predict(mod_covs, "det") 

#' 
#' **Compare with NIMBLE**
#' 
## -----------------------------------------------------------------------------
cbind(NIMBLE=apply(as.matrix(nimble_covs[,c(1:4, 2409:2410)]), 2, mean),
      unmarked=coef(mod_covs)[c(1:4,6:7)])

#' 
#' ### Model goodness-of-fit
#' 
#' Procedure applying parameteric bootstrap:
#' 
#' 1. Calculate some fit statistic for model fit to real data
#' 2. Simulate $n$ datsets from fitted model
#' 3. Calculate fit statistics for simulated datasets
#' 4. Compare real fit stat to distribution of simulated fit stats
#' 
#' If model fits well, well they should be similar.
#' 
#' One possible statistic is **sum of squared Pearson residuals (SSE)**.
#' 
#' $$ \sum{\frac{(y_{ijt} - \hat{y}_{ij})^2}{\hat{y}_{ijt}  \cdot (1 - \hat{y}_{ijt})}} $$
#' 
#' $$\hat{y}_{ijt} = \psi_{it} \cdot p_{ijt} $$
#' 
## -----------------------------------------------------------------------------
pearson2 <- function(fit, ...){
  ft <- fitted(fit) # yhat
  pr <- (fit@data@y - ft)^2 / (ft * (1-ft))
  sum(pr, na.rm=TRUE)
}

#' 
#' We can do parametric bootstrap with the `parboot` function.
#' 
## -----------------------------------------------------------------------------
gof <- suppressWarnings(parboot(mod_covs, statistic=pearson2, nsim=30))

#' 
## -----------------------------------------------------------------------------
gof
plot(gof)

#' 
#' Fit is not so good.
#' 
#' Alternative statistic: **MacKenzie-Bailey chi-square** (MacKenzie & Bailey 2004, JABES)
#' 
#' Test is based on comparing expected vs. observed frequency of specific detection histories (001, 110, 101, etc.)
#' 
## -----------------------------------------------------------------------------
library(AICcmodavg)
mb_test <- suppressWarnings(mb.gof.test(mod_covs, nsim=30))

#' 
## -----------------------------------------------------------------------------
mb_test

#' 
#' **Thoughts on GOF tests**
#' 
#' * Poor fit, possibly missing key covariates or need a random effect
#' * This approach to calculating $\hat{y}$, residuals (and thus SSE) is a bit crude
#' * M-B test struggles when there are a lot of repeated samples or a lot of `NA`s
#' * More work needed in this area
#' 
#' ### Plot covariate effects
#' 
#' To do this manually, use `predict`.
#' 
## -----------------------------------------------------------------------------
elev_rng <- range(crossbill$ele)
elev_seq <- seq(elev_rng[1], elev_rng[2], length.out=100) # NOT STANDARDIZED!!!
head(elev_seq)

#' 
## -----------------------------------------------------------------------------
nd <- data.frame(ele = elev_seq)
pr <- predict(mod_covs, type = "psi", newdata = nd, appendData = TRUE)
head(pr)

#' 
## -----------------------------------------------------------------------------
ggplot(data = pr, aes(x = ele)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) +
  geom_line(aes(y = Predicted)) +
  theme_bw(base_size=24) +
  labs(y = expression(psi~" (Predicted)"),
       x = "Elevation")+
  theme(panel.grid = element_blank())

#' 
## -----------------------------------------------------------------------------
plotEffects(mod_covs, type="psi", "ele")

#' 
## -----------------------------------------------------------------------------
plotEffects(mod_covs, type="det", "date")

#' 
#' ### Get $\psi$ for each time $t$
#' 
#' Use the `projected` function
#' 
#' Mean for each $t$ across sites:
#' 
## -----------------------------------------------------------------------------
unmarked::projected(mod_covs)

#' 
#' For each site separately:
#' 
## -----------------------------------------------------------------------------
proj <- unmarked::projected(mod_covs, mean = FALSE)
dim(proj)

#' 
## -----------------------------------------------------------------------------
proj <- t(proj[2,,])  # just "occupied" probability
proj[1:5,]

#' 
#' **Plot for first 5 sites**
#' 
## -----------------------------------------------------------------------------
cols <- rainbow(5)
plot(1:9, proj[1,1:9], ylim=c(0,0.5), type='o', pch=19, 
     col = cols[1],
     ylab="psi", xlab="time")
for (i in 2:5){
  lines(1:9, proj[i,1:9], type='o', pch=19, col=cols[i])    
}

#' 
#' **Getting CIs on these estimates**
#' 
#' Use `parboot` function for parametric bootstrapping
#' 
#' Steps:
#' 
#' * Simulate $n$ datasets from model
#' * Apply `statistic` function to each dataset
#' * Return results
#' 
## -----------------------------------------------------------------------------
proj_output <- function(mod){
  unmarked::projected(mod)[2,] # just save "occupied" only
}
pb <- suppressWarnings(parboot(mod_covs, statistic=proj_output, nsim=30))

#' 
## -----------------------------------------------------------------------------
pb

#' 
## -----------------------------------------------------------------------------
dim(pb@t.star)
head(pb@t.star)

#' 
#' **Calculate 95% CIs**
#' 
## -----------------------------------------------------------------------------
proj_lower <- apply(pb@t.star, 2, quantile, 0.025)
proj_upper <- apply(pb@t.star, 2, quantile, 0.975)
proj_est <- projected(mod_covs)[2,]

#' 
#' **Make plot**
#' 
## -----------------------------------------------------------------------------
plot(1:9, proj_est, ylim=c(0.2,0.5), type='o', pch=19,
     xlab="time", ylab="psi")
segments(1:9, proj_lower, 1:9, proj_upper)

#' 
#' ### Build a map of occupancy
#' 
#' Use matching map data in `Switzerland` dataset.
#' 
## -----------------------------------------------------------------------------
data(Switzerland)
head(Switzerland)

#' 
#' **Convert dataset to raster format**
#' 
#' (technically need to specify projection, but we won't worry about it)
#' 
## -----------------------------------------------------------------------------
library(terra)
sw_raster <- rast(Switzerland, type="xyz")
names(sw_raster) <- c("ele", "forest", "water")
plot(sw_raster)

#' 
## -----------------------------------------------------------------------------
pr_map <- predict(mod_covs, type="psi", newdata=sw_raster)
plot(pr_map)

#' 
#' ### Exercise
#' 
#' Use unmarked to fit the following models:
#' 
#' Model 1:
#' 
#' * initial occupancy: elevation + forest
#' * colonization: forest
#' * extinction: elevation
#' * detection: date

mod_exercise1 <-  colext(psiformula = ~scale(ele) + scale(forest), 
                         epsilonformula = ~scale(ele), 
                         gammaformula = ~scale(forest), 
                         pformula = ~scale(date), 
                         data = umf)

#' 
#' Model 2:
#' 
#' * initial occupancy: elevation + elevation^2 + forest
#' * colonization: forest
#' * extinction: elevation
#' * detection: date

mod_exercise2 <-  colext(psiformula = ~scale(ele) + scale(I(ele^2)) + scale(forest), 
                         epsilonformula = ~scale(ele), 
                         gammaformula = ~scale(forest), 
                         pformula = ~scale(date), 
                         data = umf)

#' 
#' If using `unmarked`, also make a map of initial occupancy using the `Switzerland` dataset for both models
#' 
#' Hint: you can use the `I()` function in formulas to run a calculation, such as
#' 
#' `~I(a^2)`
#' 
## -----------------------------------------------------------------------------


plotEffects(mod_exercise2, type="psi", "ele")



#' 
#' ## Model selection
#' 
#' Use the `fitList` and `modSel` functions.
#' 
## -----------------------------------------------------------------------------
fl <- fitList(mod_covs=mod_covs, mod_exercise1=mod_exercise1, mod_exercise2=mod_exercise2)
modSel(fl)

#' 
#' # An alternative: ubms
#' 
#' * `unmarked`-style interface
#' * Fit in Bayesian framework with `Stan`
#' 
## -----------------------------------------------------------------------------
library(ubms)

mod_stan <- stan_colext(psiformula = ~scale(ele), epsilonformula = ~1, 
                        gammaformula = ~scale(forest), pformula = ~scale(date), umf,
                        chains = 3, iter = 300, refresh=0)

#' 
## -----------------------------------------------------------------------------
traceplot(mod_stan)

#' 
## -----------------------------------------------------------------------------
mod_stan

#' 
#' **Compare `unmarked` and `ubms` output**
#' 
## -----------------------------------------------------------------------------
cbind(unmarked=coef(mod_covs), ubms=coef(mod_stan))

#' 
#' Some interesting helper functions:
#' 
## -----------------------------------------------------------------------------
plot_posteriors(mod_stan) #not skew or a line (prior bad set)

#' 
## -----------------------------------------------------------------------------
plot(mod_stan) #wreht

#' 
## -----------------------------------------------------------------------------
gof_test <- gof(mod_stan)
gof_test
plot(gof_test)

#' 
## -----------------------------------------------------------------------------
plot_effects(mod_stan, 'state')

#' 
