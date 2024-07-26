

# Here's a question: Can we put the same covariate into the state 
# and the detection model in a distance-sampling protocol ?


# Marc Kéry, 26 July 2024

# "Simulation is the answer !"

# [Perhaps we can't *prove* things by simulation, 
#  but this is the most accessible method for answering 
#  such a question for a non-statistician.]

# We first take a data simulation function from the AHMbook package, simHDS(),
#   and adapt it so that we can specify effects of the same continuous covariate
#   on both parts of the model (abundance/density and detection)

# Then, we will repeatedly simulate a data set with effects of that covariate
#   on both model parts, and analyse the data set with the hierarchical distance
#   sampling function distsamp() in unmarked.

library(unmarked)
library(AHMbook)
?simHDS


# First, we need to define a new variant of the simHDS() function.
# The standard function in the AHMbook package lets one fit effects
# of a continuous explanatory variable, habitat, on density, and of 
# another continuous covariate, wind, on the detection function scale parameter sigma.
# (That is, on the observation part of the model.)
# We will now adapt that function so that we can specify an effect of habitat
# on detection also.

# We can copy-paste the function by just typing in R the function name
# *without* the parentheses. That is:

simHDS

# In the following is the result of these changes that I make.

# -------------------- Start of R code for function definition -----------------------
simHDSX <- function(type = c("line", "point"), nsites = 100, mean.lambda = 2, 
  beta.lam.hab = 1, mean.sigma = 1, beta.sig.hab = -1, beta.sig.wind = -0.5,
  B = 3, discard0 = TRUE, show.plot = TRUE) {

  # Simulate values for covariates
  habitat <- rnorm(nsites)
  wind <- runif(nsites, -2, 2)

  # Simulate true state (abundance N)
  lambda <- exp(log(mean.lambda) + beta.lam.hab * habitat)
  N <- rpois(nsites, lambda)
  N.true <- N

  # Simulate detection model
  sigma <- exp(log(mean.sigma) + beta.sig.hab * habitat +  beta.sig.wind * wind)
  data <- NULL

  for (i in 1:nsites) {
    if (N[i] == 0) {
      data <- rbind(data, c(i, NA, NA, NA, NA))
    next
    }
    if (type == "line") {
      d <- runif(N[i], 0, B)
      p <- exp(-d * d/(2 * (sigma[i]^2)))
      y <- rbinom(N[i], 1, p)
      u <- v <- rep(NA, N[i])
      d <- d[y == 1]
      u <- u[y == 1]
      v <- v[y == 1]
      y <- y[y == 1]
    }
    if (type == "point") {
      u <- runif(N[i], 0, 2 * B)
      v <- runif(N[i], 0, 2 * B)
      d <- sqrt((u - B)^2 + (v - B)^2)
      N.true[i] <- sum(d <= B)
      p <- exp(-d * d/(2 * (sigma[i]^2)))
      pp <- ifelse(d <= B, 1, 0) * p
      y <- rbinom(N[i], 1, pp)
      u <- u[y == 1]
      v <- v[y == 1]
      d <- d[y == 1]
      y <- y[y == 1]
    }
    if (sum(y) > 0) 
      data <- rbind(data, cbind(rep(i, sum(y)), y, u, v, d))
    else data <- rbind(data, c(i, NA, NA, NA, NA))
    }
    colnames(data) <- c("site", "y", "u", "v", "d")
    if (discard0) 
      data <- data[!is.na(data[, 2]), ]
    if (show.plot) {
      if (type == "line") {
        op <- par(mfrow = c(1, 3))
        on.exit(par(op))
        tryPlot <- try({
          hist(data[, "d"], col = "lightblue", breaks = 20, 
          main = "Frequency of distances", xlab = "Distance")
          ttt <- table(data[, 1])
          n <- rep(0, nsites)
          n[as.numeric(rownames(ttt))] <- ttt
          plot(habitat, n, main = "Observed counts (n) vs. habitat", pch = 19)
          plot(wind, n, main = "Observed counts (n) vs. wind speed", pch = 19)
        }, silent = TRUE)
            if (inherits(tryPlot, "try-error")) 
                tryPlotError(tryPlot)
        }
        if (type == "point") {
            op <- par(mfrow = c(2, 2))
            on.exit(par(op))
            tryPlot <- try({
                plot(data[, "u"], data[, "v"], pch = 16, 
				  main = "Located individuals in point transects", 
                  xlim = c(0, 2 * B), ylim = c(0, 2 * B), col = data[, 
                    1], asp = 1)
                points(B, B, pch = "+", cex = 3, col = "black")
                plotrix::draw.circle(B, B, B)
                hist(data[, "d"], col = "lightblue", breaks = 20, 
                  main = "Frequency of distances", xlab = "Distance")
                ttt <- table(data[, 1])
                n <- rep(0, nsites)
                n[as.numeric(rownames(ttt))] <- ttt
                plot(habitat, n, main = "Observed counts (n) vs. habitat", pch = 19)
                plot(wind, n, main = "Observed counts (n) vs. wind speed", pch = 19)
            }, silent = TRUE)
            if (inherits(tryPlot, "try-error")) 
                tryPlotError(tryPlot)
        }
    }
    list(type = type, nsites = nsites, mean.lambda = mean.lambda, 
        beta.lam.hab = beta.lam.hab, mean.sigma = mean.sigma,
		beta.sig.hab = beta.sig.hab, beta.sig.wind = beta.sig.wind,  
        B = B, data = data, habitat = habitat, wind = wind, N = N, 
        N.true = N.true)
}
# -------------------- End of R code for function definition -----------------------

# Once we have defined the new function by executing the past ~100 lines of code,
# we can use the new function to simulate data.
# Here are two examples, one for each of the two typical types of DS protocols
# (i.e., line or point 'transect' sampling).

# Line transect protocol with explicit default argument values
str(dat <- simHDSX(type = "line", nsites = 100, mean.lambda = 2, 
  beta.lam.hab = 1, mean.sigma = 1, beta.sig.hab = -1, beta.sig.wind = -0.5,
  B = 3, discard0 = TRUE, show.plot = TRUE))


# Point transect protocol with explicit default argument values
str(dat <- simHDSX(type = "point", nsites = 100, mean.lambda = 2, 
  beta.lam.hab = 1, mean.sigma = 1, beta.sig.hab = -1, beta.sig.wind = -0.5,
  B = 3, discard0 = TRUE, show.plot = TRUE))



# The next step is that we develop the code necessary for the entire simulation by
# creating a single data set and then writing code for data formatting and 
# data analysis. Only once we have written this code, we can then (in the bottom part of this script)
# package all into a loop over many data sets and thus run the actual simulation study.

# Thus, we now write code to simulate one data set and then fitting the data-generating model
# to it using the distsamp() function in unmarked.


# Create a point transect data set
# -------------------------------
set.seed(1)
str(dat <- simHDSX(type = "point", nsites = 250, mean.lambda = 100, 
  beta.lam.hab = 1, mean.sigma = 1, beta.sig.hab = -1, beta.sig.wind = 0,
  B = 5, discard0 = TRUE, show.plot = TRUE))

# Note that here we can assume that mean.sigma and B are measured in units of 100 metres.
# In addition, mean.lambda refers to a *square* with side length 2 * B. That is, here,
# to a square with area = (2*dat$B)^2 = 100 area units.

# However, for the point transect, we have an area of pi * 5^2 = 78.53982.
# Thus, what we will estimate below as 'N' in a point-count DS model
# will be only 0.785 times that value.

# Extract distance sampling and habitat information
data <- dat$data
habitat <- dat$habitat

# Prepare the DS data
delta <- 1                          # Bin width in units of 100m
dclass <- data[,'d'] %/% delta + 1  # Convert distances to distance category data
dist.breaks <- seq(0, 100*dat$B, by = 100*delta)  # in units of metres
nD <- length(dist.breaks)-1         # Number of distance intervals
ytab <- table(c(data[,'site'], 1:dat$nsites), c(dclass, rep(NA, dat$nsites))) # Get response data
y <- matrix(c(ytab), nrow = dat$nsites, ncol = nD)
head(ytab)  ;  head(y)              # Check whether formatting still OK

# Put data into unmarked data frame
summary(umf <- unmarkedFrameDS(y = y, siteCovs = data.frame(hab = dat$habitat),
  dist.breaks = dist.breaks, survey = "point", unitsIn = "m"))

# Fit Null model
summary(fm0 <- distsamp(formula = ~ 1 ~ 1, data = umf, keyfun = "halfnorm",
               output="abund", unitsOut = "ha"))

# Fit data-generating model with habitat effect in both abundance and detection
summary(fm1 <- distsamp(formula = ~ hab ~ hab, data = umf, keyfun = "halfnorm",
               output="abund", unitsOut = "ha"))

# Call:
# distsamp(formula = ~hab ~ hab, data = umf, keyfun = "halfnorm", 
    # output = "abund", unitsOut = "ha")

# Abundance (log-scale):
            # Estimate     SE    z P(>|z|)
# (Intercept)    4.319 0.0467 92.5 0.0e+00
# hab            0.993 0.0487 20.4 2.1e-92

# Detection (log-scale):
            # Estimate     SE     z   P(>|z|)
# (Intercept)    4.606 0.0187 246.1  0.00e+00
# hab           -0.989 0.0292 -33.9 2.77e-251

# AIC: 2669.482 
# Number of sites: 250
# optim convergence code: 0
# optim iterations: 132 
# Bootstrap iterations: 0 

# Survey design: point-transect
# Detection function: halfnorm
# UnitsIn: m
# UnitsOut: ha 




# Now we are ready to package the above code for data simulation,
# data formatting and data analysis and use a for loop to
# repeat this a large number of times, e.g., 100, 1000 or more.
# It's most insightful when we vary the values of 
# the parameters and possibly also of the sample sizes
# (e.g., number of sites, number of distance classes) among
# the individual data sets. Then, we can see whether any of these
# has an effect on 'the results' (which might for instance be
# the degree of bias of the estimates or their precision).
# ----------------------------------------------------------------

# Choose number of simulated data sets
simrep <- 10000

# Set up R objects to hold the true values and the estimates and
# the indicator for convergence from function optim() (which,
# remember, is what unmarked uses internally to fit all models.
true.vals <- array(NA, dim = c(simrep, 4))
colnames(true.vals) <- names(coef(fm1))
estimates <- array(NA, dim = c(4, 4, simrep))
dimnames(estimates) <- list(names(coef(fm1)), c('Estimate', 'SE', 'z', 'P(>|z|)'), NULL)
convergence <- numeric(simrep)

# Launch simulation 
system.time(
for(k in 1:simrep){
  cat(paste('\n\n*** Iteration', k, '***\n\n')) # Counter

  # Randomly pick some values for the 4 parameters and save them
  # We pick them randomly over some range of values that appear interesting 
  # to us
  mlam <- runif(1, 50, 200)
  blam <- runif(1, -0.5, 0.5)
  msig <- runif(1, 0.5, 1.5)
  bsig <- runif(1, -0.5, 0.5)
  true.vals[k,] <- c(mlam, blam, msig, bsig)     # Save the true values
  
  # Simulate a data set using these values of the four parameters
  dat <- simHDSX(type = "point", nsites = 250, mean.lambda = mlam, 
    beta.lam.hab = blam, mean.sigma = msig, beta.sig.hab = bsig, beta.sig.wind = 0,
    B = 5, discard0 = TRUE, show.plot = FALSE)

  # Extract and format data needed for distsamp()
  # The bottom two lines contain some sausage making to ensure that we always get nD
  # number of distance intervals
  data <- dat$data
  habitat <- dat$habitat
  delta <- 1                               # Bin width in units of 100m
  dclass <- data[,'d'] %/% delta + 1       # Convert distances to distance category data
  dist.breaks <- seq(0, 100*dat$B, by = 100*delta)  # in units of metres
  nD <- length(dist.breaks)-1              # Number of distance intervals
  ytab <- table(c(data[,'site'], 1:(dat$nsites+nD)), c(dclass, rep(NA, dat$nsites), 1:nD)) # Get response data .... some sausage making here
  y <- matrix(c(ytab), nrow = dat$nsites+nD, ncol = nD)[1:dat$nsites,] # more sausage making

  # Put data into unmarked data frame
  umf <- unmarkedFrameDS(y = y, siteCovs = data.frame(hab = dat$habitat),
    dist.breaks = dist.breaks, survey = "point", unitsIn = "m")

  # Fit data-generating model with habitat effect in both abundance and detection
  fm <- distsamp(formula = ~ hab ~ hab, data = umf, keyfun = "halfnorm",
               output="abund", unitsOut = "ha")
  convergence[k] <- fm@opt$convergence     # Save optim convergence indicator

  # Save point estimates
  tmp <- summary(fm@estimates)
  estimates[,,k] <- as.matrix(rbind(tmp$state, tmp$det))
} 
)

# 10^5 take about 25 minutes on my laptop


# Remember abundance intercept is for 100 areal units in simulated data,
# but the intercept for abundance in the fitted model refers to 
# pi * dat$B^2 = 78.53982 areal units


# Visual summary of results: compare estimates with known truth
# The 1:1 line is shown in red. In an ideal world, all the estimates
# will cluster tightly around this 1:1 line (i.e., the estimates would
# close to their true values).
# The blue dashed line shows the line of best fit for the regression
# of the estimates on the truth. 
# Match of the red and the blue lines indicate unbiasedness of the 
# estimators.

par(mfrow = c(2, 2), mar = c(6,6,6,4), cex.axis = 1.5, cex.lab = 1.5, 
  cex.main = 1.5)
plot(true.vals[,1] * (pi * dat$B^2 *0.01), exp(estimates[1,1,]), xlab = 'True', ylab = 'Estimated', 
  main = 'Abundance intercept', pch = 19, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, lwd = 2, col = rgb(1,0,0,1))
abline(lm(exp(estimates[1,1,]) ~ I(true.vals[,1] * (pi * dat$B^2 *0.01))), lwd = 3, 
  lty = 3, col = rgb(0,0,1,1))

plot(true.vals[,2], estimates[2,1,], xlab = 'True', ylab = 'Estimated', 
  main = 'Abundance slope(hab)', pch = 19, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, lwd = 2, col = rgb(1,0,0,1))
abline(lm(estimates[2,1,] ~ true.vals[,2]), lwd = 3, lty = 3, col = rgb(0,0,1,0.4))

plot(100*true.vals[,3], exp(estimates[3,1,]), xlab = 'True', ylab = 'Estimated', ylim = c(40, 160),
  main = 'Det. function intercept', pch = 19, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, lwd = 2, col = rgb(1,0,0,1))
abline(lm(estimates[3,1,] ~ I(100*true.vals[,3])), lwd = 3, lty = 3, col = rgb(0,0,1,0.4))

plot(true.vals[,4], estimates[4,1,], xlab = 'True', ylab = 'Estimated', 
  main = 'Det. function slope(hab)', pch = 19, col = rgb(0,0,0,0.3), frame = FALSE)
abline(0, 1, lwd = 2, col = rgb(1,0,0,1))
abline(lm(estimates[4,1,] ~ true.vals[,4]), lwd = 3, lty = 3, col = rgb(0,0,1,0.4))

