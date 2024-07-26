
# Exercise 1a, b and c (Dail-Madsen module)


# TASK: take the data simulated just before Section 2.5.2 and the BUGS model in Section 2.5.2
# and then make the following changes to that model (NOTE: you do NOT change anything in 
# the data set, rather use exactly that data set):

# (Ex. 1a): fit a linear habitat gradient in lambda. For this, just invent
#           some random habitat covariate, such as by 'hab <- rnorm(50)'
#
# (Ex. 1b): add a year effect in detection probability
# (Ex. 1c): Make the year effects random




# To solve these exercises, we have copied all the relevant base code from the code file associated with the Chapter 2 and then we start working with that as a template


# Execute function as before to get one data set to work with
set.seed(2017, kind = "L'Ecuyer")
str(data <- simDM0(nsites = 50, nsurveys = 3, nyears = 5, lambda = 4,
    phi = 0.8, gamma = 1.5, p = 0.7))

 ----------------------------------

# Bundle data set
str(bdata <- list(C = data$y, nsites = dim(data$y)[1], nsurveys = dim(data$y)[3],
    nyears = dim(data$y)[2]))
# List of 4
# $ C       : int [1:50, 1:5, 1:3] 3 4 4 1 2 5 3 5 7 5 ...
# $ nsites  : int 50
# $ nsurveys: int 3
# $ nyears  : int 5

# Specify model in BUGS language
cat(file = "DM1.txt","
model {
  # Priors
  lambda ~ dunif(0, 100)   # Initial site-specific abundance
  phi ~ dunif(0, 1)        # Apparent survival (omega in paper/unmarked)
  gamma ~ dunif(0, 5)      # Per-capita recruitment rate
  p ~ dunif(0, 1)          # Detection probability

  # Likelihood
  for(i in 1:nsites){
    # State process: initial condition
    N[i,1] ~ dpois(lambda)
    # State process: transition model
    for(t in 1:(nyears-1)){
      S[i,t+1] ~ dbin(phi, N[i,t])     # Survival process
      # R[i,t+1] ~ dpois(gamma)        # 'absolute' recruitment = 'constant'
      R[i,t+1] ~ dpois(N[i,t] * gamma) # per-capita recruitment = 'autoreg'
      N[i,t+1] <- S[i,t+1] + R[i,t+1]
    }
    # Observation process
    for(t in 1:nyears){
      for(j in 1:nsurveys){
        C[i,t,j] ~ dbin(p, N[i,t])
      }
    }
  }
}
")

# Initial values
Nst <- apply(data$y, c(1,2), max) + 2
Nst[, 2:5] <- NA                   # cols 2:5 of N are deterministic, N <- S + R.
R1 <- apply(data$y, c(1,2), max)   # Observed max. counts + 1 as inits
R1[,1] <- NA
inits <- function(){list( lambda = runif(1, 6, 16), phi = runif(1),
    gamma = runif(1), p = runif(1), N = Nst, R = R1 + 1 )}

# Parameters monitored
params <- c("lambda", "phi", "gamma", "p")

# MCMC settings
na <- 1000 ; ni <- 25000 ; nt <- 4 ; nb <- 5000 ; nc <- 3

# Call JAGS (ART 2 min), check convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "DM1.txt", n.adapt = na, n.chains = nc,
    n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

