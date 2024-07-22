
##########################################################################################################
###
### Frequentist and Bayesian analyses of a binomial proportion 
###      with a linear model with one covariate (= logistic regression)
###
##########################################################################################################


### Some code examples taken from the books of Royle & Dorazio (2008) and Ntzoufras (2009)


# Written first 11 February 2014
# Later changes: 21 Sep 2014, 25 Dec 2014, 4 Jan 2015, 3 Sep 2017
#   22 Nov 2020, 22 Aug 2022, 28 January 2024




### Simulate some data for simple binary logistic regression
############################################################

data.fn <- function(n = 100, alpha = 0, beta = 1){
#   n <- n                      # Sample size
   x <- runif(n, -3, 3)
   alpha <- alpha                    # Regression coefficients: intercept and slope     
   beta <- beta
   logit.p <- alpha + beta * x   # Linear predictor (on logit-link scale)
   p <- plogis(logit.p)          # Success probability
   y <- rbinom(n, 1, p)          # observed binary data (0/1) as n Bernoulli trials
   plot(x, y, main = "Observed (circles) and expected response (red line)",
      xlab = "Some covariate x", pch = 16, frame = FALSE)
   lines(sort(x), p[order(x)], col = "red", lwd = 3)
   legend('topleft', paste("True values:\nIntercept:", alpha, "\nSlope:", beta),
   bty = 'n', cex = 1.2)
   return(list(n=n, alpha=alpha, beta=beta, p=p, x=x, y=y))
}

data <- data.fn()

attach(data)


# Define NLL
negLogLike <- function(beta, y, x) {   # Define likelihood function
    beta0 <- beta[1]
    beta1 <- beta[2]
    psi <- plogis(beta0 + beta1*x)
#    likelihood <- psi^y * (1-psi)^(1-y) # same as:
   likelihood <- dbinom(y, 1, psi)
    return(-sum(log(likelihood)))
}




### (0) Could do infinite number of combinations of params 'by hand'
####################################################################
# Can look at (negative) log-likelihood for two parameter sets
negLogLike(c(0,0), y, x)
negLogLike(c(-3,2), y, x) # Lower is better!
negLogLike(c(10,0), y, x) # Lower is better!




### (1) Brute force: do a grid search for MLEs
##############################################
(mat <- as.matrix(expand.grid(seq(-10,10,0.1), seq(-10,10,0.1))))  ### Can vary resolution
# mat <- as.matrix(expand.grid(seq(-10,10,0.01), seq(-10,10,0.01)))  # This takes a while ....
str(mat)     # Try 40,000 combos of possible values for intercept and slope

nll <- array(NA, dim = nrow(mat))
system.time(
  for (i in 1:nrow(mat)){
    nll[i] <- negLogLike(mat[i,], y, x)
  }
)
which(nll == min(nll))           # Number of minimum value
mat[which(nll == min(nll)),]     # Values of params that make nll minimal




### (2) Maximum likelihood (numerical optimization of function)
###############################################################

# Minimize function for observed data and return MLE
inits <- c(beta0=0, beta1=0)
(out <- optim(inits, negLogLike, y=y, x=x, hessian=TRUE) )
mles <- out$par     # MLEs are pretty close to truth

# Make a table with estimates, SEs, and 95% CI
mle.table <- data.frame(Est=mles,
                        SE = sqrt(diag(solve(out$hessian))))
mle.table$lower <- mle.table$Est - 1.96*mle.table$SE
mle.table$upper <- mle.table$Est + 1.96*mle.table$SE
mle.table




### (3) Maximum likelihood (using canned function glm())
########################################################

# Estimate parameter on link scale
summary(fm <- glm(y ~ x, family = binomial, data = data))
pred <- predict(fm, type = "response", se = TRUE)
lines(sort(data$x), pred$fit[order(data$x)], col = "blue", lwd = 3, lty = 2)





# Plots of likelihood landscape with estimates of minimum
#########################################################
# Using colors in raster package
library(raster)
r <- rasterFromXYZ(data.frame(x = mat[,1], y = mat[,2], z = nll))
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

plot(r, col = mapPalette(100), main = "Negative log-likelihood \n(blue - glm MLEs, red - truth, point - grid search MLEs)", 
	xlab = "Intercept", ylab = "Slope", asp = 1)
(rrr <- range(nll[nll != 'Inf']))
contour(r, add = TRUE, levels = round(seq(rrr[1], rrr[2], 20)))
#contour(r, add = TRUE, levels = seq(60, 700, 10))

points(mat[which(nll == min(nll)),1], mat[which(nll == min(nll)),2], pch = 16)  # MLEs from grid search

# Compare with 'exact' maximum likelihood solution from glm()
fm <- glm(y ~ x, family = binomial, data = data)
abline(v = coef(fm)[1], col = "blue", lwd = 2)
abline(h = coef(fm)[2], col = "blue", lwd = 2)

# Add truth (values used to simulate the data)
abline(v = data$alpha, col = "red", lty = "dashed", lwd = 2)
abline(h = data$beta, col = "red", lty = "dashed", lwd = 2)


# can someone come up with an embellished version of this graph ?
# perhaps as a 3d graph ?





### Random walk MCMC for logistic regression (Metropolis algorithm)
##################################################################

# Code example partly taken from p. 58 of Ntzoufras (2009)

### Define function that does RW-MCMC
mcmc.fn <- function(y = y, x = x,            # The data
niters = 2500,                               # Number of MCMC iterations
mu.beta = c(0, 0),                           # Mean and sd of flat normal priors on two params 
s.beta = c(100, 100),
prop.s = c(0.2, 0.2),                        # SD of Normal proposal (determines step length of random walk)
beta0 = c(0, 0)){                            # Initial values

# Function that does random-walk MCMC for a binomial proportion (beta = logit intercept and slope)
# Code taken from p. 58 in the book by Ntzoufras (2009)

# Initalize calculations
  start.time = Sys.time()                    # Set timer
  beta.post <- array(NA, dim = c(niters,2))  # Set up matrix for posterior draws of beta
  acc.counter <- 0                           # Initialize counters for acceptance
  current.beta <- beta0                      # Initial value for theta

# Start MCMC algorithm
precision <- 700                             # Numerical stabilisation

for (t in 1:niters){                         # Repeat niters times

   if(t %% 500 ==0 )                        # report progress
   cat("iter", t, "\n")

   prop.beta <- rnorm(2, current.beta, prop.s)  # Propose new values for beta to compare against current values

   # calculate the linear predictors
   current.eta <- current.beta[1] + current.beta[2]*x # lp for current param values
   prop.eta <- prop.beta[1] + prop.beta[2]*x          # lp for proposed param values

   # if the linear predictor > precision then set it equal to the precision (to avoid p=0, p=1)
   current.eta[current.eta > precision] <- precision
   current.eta[current.eta < -precision] <- -precision
   prop.eta[prop.eta > precision] <- precision
   prop.eta[prop.eta < -precision] <- -precision

   # Compare likelihood times prior (which is proportional to posterior density) between the new (proposed) and the old (current values)
   loga <- (sum(y * prop.eta - log(1 + exp(prop.eta))) 
   - sum(y * current.eta - log(1 + exp(current.eta)))          
   + sum(dnorm(prop.beta, mu.beta, s.beta, log=TRUE)) 
   - sum(dnorm(current.beta, mu.beta, s.beta, log=TRUE)) )

   u <- log(runif(1))
   if (u < loga){ current.beta <- prop.beta    # If new (proposed) beta leads to a higher product of likelihood x prior,
                                               # then take new beta as the new current.beta (with p=1)
                                               #
                                               # If it leads to a smaller product of like x prior, then take the new beta with probability a
                                               # and keep the old value of beta with probability (1-a)
   acc.counter <- acc.counter + 1              #  Counts the number of acceptances
   }	        
#   browser()                                  # If unhashed, allows to inspect values of logu and loga at each iteration
   beta.post[t,] <- current.beta

   plot(beta.post[1:t,1], beta.post[1:t,2], type = "b",   # Cool animation
      xlab = "Intercept (alpha)", ylab = "Slope (beta)",
      main = "MCMC trajectory (red: truth, blue: sample means)")

abline(v = alpha, col = "red", lwd = 1)
abline(h = data$beta, col = "red", lwd = 1)
abline(h = mean(beta.post[1:t,2]), col = "blue", lwd = 1)
abline(v = mean(beta.post[1:t,1]), col = "blue", lwd = 1)

points(beta.post[1,1], beta.post[1,2], add = TRUE)
}

devAskNewPage(ask = TRUE)

acc.prob <- acc.counter/niters		# Acceptance probability
cat("Acceptance probability:", round(acc.prob, 2), "\n")
end.time = Sys.time()			# Stop time
elapsed.time = round(difftime(end.time, start.time, units='secs'), 2)  # Compute elapsed time
cat(paste(paste('Posterior samples drawn in ', elapsed.time, sep=''), ' seconds\n\n', sep='')) # Output run time

par(mfrow = c(3,2), mar = c(3,2,4,1))          # Plots
plot(beta.post[,1], type = "l",                # Plot 1: time-series plot of logit(alpha) 
   ylab = "logit(Intercept))", main = "Intercept, logit(alpha): Time-series")
abline(h = alpha, col = "red")		       # True value in data generation
abline(h = mean(beta.post[,1]), col = "blue")  # Posterior mean

plot(beta.post[,2], type = "l",                # Plot 2: time-series plot of logit(beta)
   ylab = "logit(Intercept))", main = "Slope, logit(beta): Time-series ")
abline(h = data$beta, col = "red")		       # True value in data generation
abline(h = mean(beta.post[,2]), col = "blue")  # Posterior mean

hist(beta.post[,1], breaks = 100, col = "grey",# Plot 3: Histogram of posterior of logit(alpha) with smoothed line
main = "Posterior (with truth red and smooth blue)", freq = FALSE)	
smooth <- density(beta.post[,1], adjust = 2)
abline(v = alpha, col = "red", lwd = 3)        # True value in data generation
lines(smooth$x, smooth$y, type = 'l', lwd = 2, col = "blue")

hist(beta.post[,2], breaks = 100, col = "grey",# Plot 4: Histogram of posterior of logit(beta) with smoothed line
   main = "Posterior (with truth red and smooth blue)", freq = FALSE)
smooth <- density(beta.post[,2], adjust = 2)
abline(v = data$beta, col = "red", lwd = 3)         # True value in data generation
lines(smooth$x, smooth$y, type = 'l', lwd = 2, col = "blue")

plot(acf(beta.post[,1], plot = FALSE), main = "Autocorrelation function", lwd = 3) # Plot 5: Autocorrelation function plot
plot(acf(beta.post[,2], plot = FALSE), main = "Autocorrelation function", lwd = 3) # Plot 6: Autocorrelation function plot

par(mfrow = c(1,1))

# Return MCMC samples
return(beta.post = beta.post)

} # end of function definition



### Execute function
####################

# NOTE: RStudio users may have to use real R to experience the full power of R !

#detach(data)

data <- data.fn()
attach(data)
out <- mcmc.fn(x=data$x, y=data$y, beta0 = c(2, 2), niters = 1500)
#detach(data)



# ********************************************************************** 
#
#   NOTE: Have to close graphics window after every run of function !!!
#
# ********************************************************************** 

#### Different MCMC settings (some of them crazy !)
# Effects of variance of proposal distribution (step length in parameter space exploration)
##########################################################################################
out <- mcmc.fn(y = data$y, x = data$x, niters = 1000, mu.beta = c(0,0), s.beta = c(100, 100),
   prop.s = c(0.2, 0.2), beta0 = c(0,0))       # Small variance in proposal distribution
# close window....

#out <- mcmc.fn(y = data$y, x = data$x, niters = 1000, mu.beta = c(0,0), s.beta = c(100, 100),
#   prop.s = c(0.4, 0.4), beta0 = c(0,0))       # Slightly larger variance in proposal distribution
# close window....

out <- mcmc.fn(y = data$y, x = data$x, niters = 1000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(1, 1), beta0 = c(0,0))          # Larger variance in proposal distribution
# close window....

#out <- mcmc.fn(y = data$y, x = data$x, niters = 1000, mu.beta = c(0,0), s.beta = c(100, 100), 
#   prop.s = c(2, 2), beta0 = c(0,0))          # Still larger variance in proposal distribution
# close window....

out <- mcmc.fn(y = data$y, x = data$x, niters = 1000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(5, 5), beta0 = c(0,0))           # Huge variance in proposal distribution
# close window....




# Check power of search
#######################
### Inits far out: this is 'beta0'
out <- mcmc.fn(y = data$y, x = data$x, niters = 1000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(0.2, 0.2), beta0 = c(10,10))


### Inits extremely far out, short chains
out <- mcmc.fn(y = data$y, x = data$x, niters = 1000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(0.2, 0.2), beta0 = c(100, 100))


### Inits extremely far out, longer chains
out <- mcmc.fn(y = data$y, x = data$x, niters = 3000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(0.2, 0.2), beta0 = c(100, 100))


### Inits extremely far out, chains NOT longer, but greater step length
out <- mcmc.fn(y = data$y, x = data$x, niters = 1000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(0.6, 0.6), beta0 = c(100, 100))
   
   
out <- mcmc.fn(y = data$y, x = data$x, niters = 3000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(0.2, 0.2), beta0 = c(-100, -100))




### Inits extremely far out, longer chains                       ##### ************ Cool !
out <- mcmc.fn(y = data$y, x = data$x, niters = 2000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(0.2, 0.2), beta0 = c(0, 100))


out <- mcmc.fn(y = data$y, x = data$x, niters = 2000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(0.2, 0.2), beta0 = c(100, 0))                        ##### ************ Cool U-shape !


out <- mcmc.fn(y = data$y, x = data$x, niters = 2000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(0.2, 0.2), beta0 = c(-100, 100))         # **** Straaaaaiiiiiight line !



# The longer the chains, the less important becomes initial
out <- mcmc.fn(y = data$y, x = data$x, niters = 10000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(0.2, 0.2), beta0 = c(-100, 100))         # **** Posterior means moves in more and more




# ******* Craaaaazy !
system.time(
out <- mcmc.fn(y = data$y, x = data$x, niters = 30000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(0.2, 0.2), beta0 = c(500, 500))
   )
# ..... and yet, it does find the solution eventually !!!!!
   


# ******* Craaaaazy !         (*** takes about 45 mins to finish  ***)
system.time(
out <- mcmc.fn(y = data$y, x = data$x, niters = 25000, mu.beta = c(0,0), s.beta = c(100, 100), 
   prop.s = c(0.2, 0.2), beta0 = c(0, 1000))
   )
# ..... and yet, it does find the solution eventually !!!!!








###### Execute all together (for Bayes and non-Bayes estimates)
##############################################################


#detach(data)

data <- data.fn()
attach(data)
out <- mcmc.fn(x=data$x, y=data$y, niters = 1000)


# Plot data and predicted response from MLE and Bayesian posterior inference
par(mfrow = c(1,1), mar = c(5,4,2,2))
# Data
plot(data$x, data$y, main = "Observed (circles) and expected response (red line)",
   xlab = "Covariate x")
lines(sort(data$x), data$p[order(data$x)], col = "red")

# MLEs (with function glm)
fm <- glm(y ~ x, family = binomial, data = data)
pred <- predict(fm, type = "response", se = TRUE)
lines(sort(data$x), pred$fit[order(data$x)], col = "blue")

# Bayesian posterior estimates
pm <- apply(out, 2, mean)          # Posterior means of intercept and slope
pred.pm <- plogis(pm[1] + pm[2] * data$x)
lines(sort(data$x), pred.pm[order(data$x)], col = "black")
legend(-3, 0.9, c('True values in data generation','MLEs', 'Bayesian posterior mean'), 
   col=c("red", "blue", "black"), lty = 1, cex=0.8, border = FALSE)

legend(0.4, 0.2, "MCMC is just like magic .....", cex = 0.8, col = "red", border = F)


### Could also add MLE's using own function minimisation code from above




# Plots of likelihood landscape with estimates of minimum (frequentist AND Bayesian)
####################################################################################
# Using colors in raster package

par(mfrow = c(1,1), mar = c(5,4,3,2), cex.main = 1)

library(raster)
r <- rasterFromXYZ(data.frame(x = mat[,1], y = mat[,2], z = nll))
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r, col = mapPalette(100), main = "Negative log-likelihood (red - truth, black - MLE (grid search),
   blue - MLE (glm), green - Post. mean)", xlab = "Intercept parameter", ylab = "Slope parameter")
(rrr <- range(nll[nll != 'Inf']))
contour(r, add = TRUE, levels = round(seq(rrr[1], rrr[2], 20)))
#contour(r, add = TRUE, levels = seq(60, 700, 10))

points(mat[which(nll == min(nll)),1], mat[which(nll == min(nll)),2], pch = 16)  # MLEs from grid search

# Compare with 'exact' maximum likelihood solution from glm()
fm <- glm(y ~ x, family = binomial, data = data)
abline(v = coef(fm)[1], col = "blue", lwd = 2)
abline(h = coef(fm)[2], col = "blue", lwd = 2)


# Add truth (values used to simulate the data)
abline(v = alpha, col = "red", lty = "dashed", lwd = 2)
abline(h = beta, col = "red", lty = "dashed", lwd = 2)

# Add Bayesian posterior means
abline(v = pm[1], col = "green", lty = "dashed", lwd = 2)
abline(h = pm[2], col = "green", lty = "dashed", lwd = 2)









###############################################
###############################################
###############################################
###############################################
###############################################
###############################################
########
######## HAVE YET TO DO THIS HERE BELOW .....
########
###############################################
###############################################
###############################################
###############################################
###############################################
###############################################





### Bayesian analysis using WinBUGS
###################################

# Write text file with model description
sink("model.txt")
cat("
model {
# Prior
 p ~ dbeta(shape1, shape2)
  

# Likelihood
y ~ dbin(p, N) 
}
",fill=TRUE)
sink()


# Bundle data
r = 20
n = 50
shape1 = 1
shape2 = 1
win.data <- list(y = r, N = n, shape1 = shape1, shape2 = shape2)


# Inits function
inits <- function(){ list(p=runif(1))}


# Parameters to estimate
params <- c("p")


# MCMC settings
nc = 3  ;  ni=2500  ;  nb=500  ;  nt=2


# Start Gibbs sampler
out.1.1 <- bugs(data = win.data, inits = inits, parameters = params, model = "model.txt", 
n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, debug = FALSE, bugs.directory = bugs.dir)
print(out.1.1, dig = 3)


# List of posterior samples
out.1.1$sims.list$p
mean(out.1.1$sims.list$p)
sd(out.1.1$sims.list$p)
quantile(out.1.1$sims.list$p, probs = c(0.025, 0.975))


# Show histo and smoothed estimate
hist(out.1.1$sims.list$p, breaks = 100, col = "grey", ylab = "", xlab = "Detection probability", cex = 1.5,
main = "Histogram of posterior samples", freq = FALSE)
smooth <- density(out.1.1$sims.list$p, adjust = 2)
lines(smooth$x, smooth$y, type = 'l', lwd = 3, col = "blue")


# Plot posterior distribution
smooth <- density(out.1.1$sims.list$p, adjust = 2)
plot(smooth$x, smooth$y, type = 'l', lwd = 3, col = "blue", xlab = "Parameter value", ylab = "", xlim = c(0,1),
main = "Prior and posterior in tadpole example", las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
post.mean <- smooth$x[which(smooth$y == max(smooth$y))]
# abline(v = post.mean, col = "blue", lwd = 2)

# Add scaled likelihood function values
expand.factor <- max(smooth$y) / max(lik)
lines(psi, expand.factor * lik, lwd = 3, col = "red")
# abline(v = 0.4, col = "red", lwd = 2, lty = "dashed")

# Add prior to plot
x <- seq(0,1,0.01)
lines(x, dbeta(x, 1,1), lwd = 3, col = "green")
legend(0.6, 5, c('Likelihood function (scaled)', 'Posterior dist.', 
   'Prior dist.'), col=c("red", "blue", "green"), lty = 5, lwd = 5, cex = 0.8)









































###### ARCHIVE FROM TADPOLE EXAMPLE




### (1) Maximum likelihood (brute force approach)
#################################################

# = "just try out what's gives the highest function value")
# Code example taken from p. 41 of Royle & Dorazio (2008)


# Define the data
r <- 20
N <- 50

# Do calculations
eps <- 0.0001	# Step length
psi <- seq(eps, 1-eps, by = eps) 		# Values of psi to try out
logLike = r*log(psi) + (N-r)*log(1-psi)		# Value of log-likelihood for each trial value
lik <- exp(logLike)				# Value of likelihood
MLE <- psi[logLike == max(logLike)]		# Determine MLE

# Return results
cat("Maximum likelihood estimate of psi: ", MLE, "\n\n")
plot(psi, logLike, ylab = "loglikelihood", xlab = "Value of parameter psi", las = 1, type = "l", lwd = 3)
abline(v = MLE, lwd = 3, col = "blue")




### (2) Maximum likelihood (numerical optimization of function)
###############################################################

# Code example taken from p. 42 of Royle & Dorazio (2008)


# Define the data
r <- 20
N <- 50

# Define negative log-likelihood function
nll <- function(psi) -dbinom(r, size = N, prob = psi, log = TRUE)

# Minimize function for observed data and return MLE
fit <- optim(par = 0.5, fn = nll, method = "BFGS")
fit
cat("Maximum likelihood estimate of psi: ", fit$par, "\n\n")




### (3) Maximum likelihood (using canned function glm())
########################################################

# Define the data
r <- 20
N <- 50

# Estimate parameter on link scale
fm <- glm(cbind(20,30) ~ 1, family = binomial)
summary(fm)
predict(fm, type = "response", se = TRUE)











### Random walk MCMC for binomial proportion
############################################

# Code example taken from p. 48 of Ntzoufras (2009)


### Define function that does RW-MCMC
mcmc.fn <- function(y = 399, N = 845, niters = 2500, mu.theta = 0, s.theta = 100, 
prop.s = 0.1, current.theta = 0){

# Function that does random-walk MCMC for a binomial proportion (theta = logit(p))
# Code taken from p. 48 in the book by Ntzoufras (2009)

# Initalize calculations
  start.time = Sys.time()		     # Set timer
  y <- y ; N <- N			     # Observed data
  niters <- niters			     # Number of MCMC iterations		
  mu.theta <- mu.theta; s.theta <- s.theta   # "Flat normal" prior
  prop.s <- prop.s			     # Proposal parameter (determines step length of random walk)
  theta <- numeric(niters)		     # Set up vector for posterior draws of theta
  acc.counter <- 0			     # Initialize counter for acceptance
  current.theta <- current.theta	     # Initial value for theta

# Start MCMC algorithm
for (t in 1:niters){			# Repeat niters times
   prop.theta <- rnorm(1, current.theta, prop.s)		# Propose a value prop.theta
   loga <- (( prop.theta * y - N * log(1 + exp(prop.theta)))	# Compare likelihood times prior (which is proportional to posterior density)
   - (current.theta * y - N * log(1 + exp(current.theta)))	# between the new (proposed) and the old (current value)
   + dnorm(prop.theta, mu.theta, s.theta, log = TRUE)
   - dnorm(current.theta, mu.theta, s.theta, log = TRUE) )

   u <- log(runif(1))
   if (u < loga){ current.theta <- prop.theta		# If new (proposed) theta leads to a higher product of likelihood times prior,
   							# then take it as new current.theta with prob a
							# Otherwise keep the old value of theta
							
   acc.counter <- acc.counter + 1		. 	#  Counts the number of acceptances
	        }
#   browser()						# If unhashed, allows to inspect values of logu and loga at each iteration
   theta[t] <- current.theta
}
p <- plogis(theta)			# Compute p from logit(p)
acc.prob <- acc.counter/niters		# Acceptance probability
cat("Acceptance probability:", round(acc.prob, 2), "\n")
end.time = Sys.time()			# Stop time
elapsed.time = round(difftime(end.time, start.time, units='secs'), 2)  # Compute elapsed time
cat(paste(paste('Posterior samples drawn in ', elapsed.time, sep=''), ' seconds\n\n', sep='')) # Output run time

par(mfrow = c(2,2))			# Plots of theta=logit(p) and of p
plot(theta, type = "l", ylab = "theta (=logit(p))")		# Plot 1: time-series plot of theta = logit(p)
plot(p, type = "l", ylim = c(0,1))				# Plot 2: time-series plot of p
abline(h = y/N, col = "red")		# Add maximum likelihood estimate
abline(h = mean(p), col = "blue")	# Add posterior mean
hist(p, breaks = 100, col = "grey", main = "", freq = FALSE)	# plot 3: Histogram of posterior with smoothed line
smooth <- density(p, adjust = 2)
lines(smooth$x, smooth$y, type = 'l', lwd = 2, col = "blue")
plot(acf(p, plot = FALSE), main = "", lwd = 3)			# Plot 4: Autocorrelation function plot
} # end of function



### Execute function
####################
mcmc.fn(y = 20, N = 50, niters = 25000, mu.theta = 0, s.theta = 100, 
prop.s = 0.1, current.theta = 0)

mcmc.fn(y = 20, N = 50, niters = 2500, mu.theta = 0, s.theta = 100, 
prop.s = 1, current.theta = 0)

mcmc.fn(y = 20, N = 50, niters = 2500, mu.theta = 0, s.theta = 100, 
prop.s = 2, current.theta = 0)

mcmc.fn(y = 20, N = 50, niters = 2500, mu.theta = 0, s.theta = 100, 
prop.s = 100, current.theta = 0)

mcmc.fn(y = 20, N = 50, niters = 2500, mu.theta = 0, s.theta = 100, 
prop.s = 0.1, current.theta = 10)






### Bayesian analysis using WinBUGS
###################################

# Write text file with model description
sink("model.txt")
cat("
model {
# Prior
 p ~ dbeta(shape1, shape2)
  

# Likelihood
y ~ dbin(p, N) 
}
",fill=TRUE)
sink()


# Bundle data
r = 20
n = 50
shape1 = 1
shape2 = 1
win.data <- list(y = r, N = n, shape1 = shape1, shape2 = shape2)


# Inits function
inits <- function(){ list(p=runif(1))}


# Parameters to estimate
params <- c("p")


# MCMC settings
nc = 3  ;  ni=2500  ;  nb=500  ;  nt=2


# Start Gibbs sampler
out.1.1 <- bugs(data = win.data, inits = inits, parameters = params, model = "model.txt", 
n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, debug = FALSE, bugs.directory = bugs.dir)
print(out.1.1, dig = 3)


# List of posterior samples
out.1.1$sims.list$p
mean(out.1.1$sims.list$p)
sd(out.1.1$sims.list$p)
quantile(out.1.1$sims.list$p, probs = c(0.025, 0.975))


# Show histo and smoothed estimate
hist(out.1.1$sims.list$p, breaks = 100, col = "grey", ylab = "", xlab = "Detection probability", cex = 1.5,
main = "Histogram of posterior samples", freq = FALSE)
smooth <- density(out.1.1$sims.list$p, adjust = 2)
lines(smooth$x, smooth$y, type = 'l', lwd = 3, col = "blue")


# Plot posterior distribution
smooth <- density(out.1.1$sims.list$p, adjust = 2)
plot(smooth$x, smooth$y, type = 'l', lwd = 3, col = "blue", xlab = "Parameter value", ylab = "", xlim = c(0,1),
main = "Prior and posterior in tadpole example", las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
post.mean <- smooth$x[which(smooth$y == max(smooth$y))]
# abline(v = post.mean, col = "blue", lwd = 2)

# Add scaled likelihood function values
expand.factor <- max(smooth$y) / max(lik)
lines(psi, expand.factor * lik, lwd = 3, col = "red")
# abline(v = 0.4, col = "red", lwd = 2, lty = "dashed")

# Add prior to plot
x <- seq(0,1,0.01)
lines(x, dbeta(x, 1,1), lwd = 3, col = "green")
legend(0.6, 5, c('Likelihood function (scaled)', 'Posterior dist.', 
   'Prior dist.'), col=c("red", "blue", "green"), lty = 5, lwd = 5, cex = 0.8)







### Some informative prior examples
a <- shape1
b <- shape2
m <- a / (a+b)
s <- sqrt(a*b/((a+b)^2 *(a+b+1)))

plot(x, dbeta(x, 10, 5), lwd = 3, col = "green")

plot(x, dbeta(x, 40, 60), lwd = 3, col = "green")



# Example with informative prior and same location

# Bundle data
shape1 = 4
shape2 = 6
win.data <- list(y = r, N = n, shape1 = shape1, shape2 = shape2)

# Start Gibbs sampler
out.4.6 <- bugs(data = win.data, inits = inits, parameters = params, model = "model.txt", 
n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, debug = FALSE, bugs.directory = bugs.dir)
print(out.4.6, dig = 3)

# Plot posterior distribution
smooth <- density(out.4.6$sims.list$p, adjust = 2)
plot(smooth$x, smooth$y, type = 'l', lwd = 3, col = "blue", xlab = "Parameter value", ylab = "", xlim = c(0,1),
main = "Prior and posterior in tadpole example", las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
post.mean <- smooth$x[which(smooth$y == max(smooth$y))]
# abline(v = post.mean, col = "blue", lwd = 2)

# Add scaled likelihood function values
expand.factor <- max(smooth$y) / max(lik)
lines(psi, expand.factor * lik, lwd = 3, col = "red")
# abline(v = 0.4, col = "red", lwd = 2, lty = "dashed")

# Add prior to plot
x <- seq(0,1,0.01)
lines(x, dbeta(x, shape1, shape2), lwd = 3, col = "green")
legend(0.6, 5, c('Likelihood function (scaled)', 'Posterior dist.', 'Prior dist.'), 
col=c("red", "blue", "green"), lty = 5, lwd = 5, cex = 0.8)
### Not much effect of informative prior





# Bundle data
shape1 = 40
shape2 = 60
win.data <- list(y = r, N = n, shape1 = shape1, shape2 = shape2)

# Start Gibbs sampler
out.40.60 <- bugs(data = win.data, inits = inits, parameters = params, model = "model.txt", 
n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, debug = FALSE, bugs.directory = bugs.dir)
print(out.40.60, dig = 3)

# Plot posterior distribution
smooth <- density(out.40.60$sims.list$p, adjust = 2)
plot(smooth$x, smooth$y, type = 'l', lwd = 3, col = "blue", xlab = "Parameter value", ylab = "", xlim = c(0,1),
main = "Prior and posterior in tadpole example", las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
post.mean <- smooth$x[which(smooth$y == max(smooth$y))]
# abline(v = post.mean, col = "blue", lwd = 2)

# Add scaled likelihood function values
expand.factor <- max(smooth$y) / max(lik)
lines(psi, expand.factor * lik, lwd = 3, col = "red")
# abline(v = 0.4, col = "red", lwd = 2, lty = "dashed")

# Add prior to plot
x <- seq(0,1,0.01)
lines(x, dbeta(x, shape1, shape2), lwd = 3, col = "green")
legend(0.6, 5, c('Likelihood function (scaled)', 'Posterior dist.', 'Prior dist.'), 
col=c("red", "blue", "green"), lty = 5, lwd = 5, cex = 0.8)
### Informative prior has effect as would an additional sample of data




# Bundle data
shape1 = 60
shape2 = 40
win.data <- list(y = r, N = n, shape1 = shape1, shape2 = shape2)

# Start Gibbs sampler
out.60.40 <- bugs(data = win.data, inits = inits, parameters = params, model = "model.txt", 
n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, debug = FALSE, bugs.directory = bugs.dir)
print(out.60.40, dig = 3)

# Plot posterior distribution
smooth <- density(out.60.40$sims.list$p, adjust = 2)
plot(smooth$x, smooth$y, type = 'l', lwd = 3, col = "blue", xlab = "Parameter value", ylab = "", xlim = c(0,1),
main = "Prior and posterior in tadpole example", las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
post.mean <- smooth$x[which(smooth$y == max(smooth$y))]
# abline(v = post.mean, col = "blue", lwd = 2)

# Add scaled likelihood function values
expand.factor <- max(smooth$y) / max(lik)
lines(psi, expand.factor * lik, lwd = 3, col = "red")
# abline(v = 0.4, col = "red", lwd = 2, lty = "dashed")

# Add prior to plot
x <- seq(0,1,0.01)
lines(x, dbeta(x, shape1, shape2), lwd = 3, col = "green")
legend(0.6, 8, c('Likelihood function (scaled)', 'Posterior dist.', 'Prior dist.'), 
col=c("red", "blue", "green"), lty = 5, lwd = 5, cex = 0.8)









##### All four in one plot
##########################

par(mfrow = c(2,2))

# Uninformative priors
######################
shape1 = 1
shape2 = 1

# Plot posterior distribution
smooth <- density(out.1.1$sims.list$p, adjust = 2)
plot(smooth$x, smooth$y, type = 'l', lwd = 3, col = "blue", xlab = "", ylab = "Density", xlim = c(0,1),
main = "Uninformative prior", las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
post.mean <- smooth$x[which(smooth$y == max(smooth$y))]

# Add scaled likelihood function values
expand.factor <- max(smooth$y) / max(lik)
lines(psi, expand.factor * lik, lwd = 3, col = "red")

# Add prior to plot
x <- seq(0,1,0.01)
lines(x, dbeta(x, 1,1), lwd = 3, col = "green")
legend(0.6, 5, c('Likelihood function', 'Posterior dist.', 'Prior dist.'), 
col=c("red", "blue", "green"), lty = 5, lwd = 3, cex = 0.8)


# Informative priors A
######################
shape1 = 4
shape2 = 6

# Plot posterior distribution
smooth <- density(out.4.6$sims.list$p, adjust = 2)
plot(smooth$x, smooth$y, type = 'l', lwd = 3, col = "blue", xlab = "", ylab = "", xlim = c(0,1),
main = "Informative prior A", las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
post.mean <- smooth$x[which(smooth$y == max(smooth$y))]

# Add scaled likelihood function values
expand.factor <- max(smooth$y) / max(lik)
lines(psi, expand.factor * lik, lwd = 3, col = "red")

# Add prior to plot
x <- seq(0,1,0.01)
lines(x, dbeta(x, shape1, shape2), lwd = 3, col = "green")
#legend(0.6, 5, c('Likelihood function (scaled)', 'Posterior dist.', 'Prior dist.'), 
# col=c("red", "blue", "green"), lty = 5, lwd = 5, cex = 0.8)



# Informative priors B
######################
# Bundle data
shape1 = 40
shape2 = 60

# Plot posterior distribution
smooth <- density(out.40.60$sims.list$p, adjust = 2)
plot(smooth$x, smooth$y, type = 'l', lwd = 3, col = "blue", xlab = "Parameter value", ylab = "Density", xlim = c(0,1),
main = "Informative prior B", las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
post.mean <- smooth$x[which(smooth$y == max(smooth$y))]

# Add scaled likelihood function values
expand.factor <- max(smooth$y) / max(lik)
lines(psi, expand.factor * lik, lwd = 3, col = "red")

# Add prior to plot
x <- seq(0,1,0.01)
lines(x, dbeta(x, shape1, shape2), lwd = 3, col = "green")
#legend(0.6, 5, c('Likelihood function (scaled)', 'Posterior dist.', 'Prior dist.'), 
# col=c("red", "blue", "green"), lty = 5, lwd = 5, cex = 0.8)



# Informative priors C
######################

# Bundle data
shape1 = 60
shape2 = 40

# Plot posterior distribution
smooth <- density(out.60.40$sims.list$p, adjust = 2)
plot(smooth$x, smooth$y, type = 'l', lwd = 3, col = "blue", xlab = "Parameter value", ylab = "", xlim = c(0,1),
main = "Informative prior C", las = 1, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
post.mean <- smooth$x[which(smooth$y == max(smooth$y))]

# Add scaled likelihood function values
expand.factor <- max(smooth$y) / max(lik)
lines(psi, expand.factor * lik, lwd = 3, col = "red")

# Add prior to plot
x <- seq(0,1,0.01)
lines(x, dbeta(x, shape1, shape2), lwd = 3, col = "green")
#legend(0.68, 8, c('Likelihood function (scaled)', 'Posterior dist.', 'Prior dist.'), 
# col=c("red", "blue", "green"), lty = 5, lwd = 5, cex = 0.8)

































































#################### ARCHIVE ###################################



# Original, with counter of acceptances
y <- 399
N <- 845
iters <- 25000
mu.theta <- 0
s.theta <- 100
prop.s <- 0.35
theta <- numeric(iters)
acc.prob <- 0
current.theta <- 0
for (t in 1:iters){
   prop.theta <- rnorm(1, current.theta, prop.s)
   loga <- (( prop.theta * y - N * log(1 + exp(prop.theta)))
   - (current.theta * y - N * log(1 + exp(current.theta)))
   + dnorm(prop.theta, mu.theta, s.theta, log = TRUE)
   - dnorm(current.theta, mu.theta, s.theta, log = TRUE) )
   u <- runif(1)
   u <- log(u)
   if (u < loga) { current.theta <- prop.theta
                    acc.prob <- acc.prob + 1
   }
   theta[t] <- current.theta
}

p <- plogis(theta)
par(mfrow = c(2,1))
plot(theta, type = "l")
plot(p, type = "l", ylim = c(0,1))
abline(h = y/N, col = "red")
abline(h = mean(p), col = "blue")
acc.prob
