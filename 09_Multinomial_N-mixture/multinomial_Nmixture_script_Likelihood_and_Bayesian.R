
#
# Part 1, demonstration on simulating multinomial data for capture-recapture type studies
#

set.seed(2015)                         # Initialize RNG
 
# Simulate covariate values and local population size for each point
x <- rnorm(100)
N <- rpois(100, lambda=exp(-1 + 1*x) ) # Intercept and slope equal to 1
table(N)                               # Summarize; it looks like a territorial birds 
#N
# 0  1  2  3  4  5  6
#72 17  6  1  2  1  1
 

 
# Define detection probabilities (p) for both observers
p1 <- 0.8
p2 <- 0.6
 
# Construct the multinomial cell probabilities (pi)
cellprobs <- c(p1*p2, p1*(1-p2), (1-p1)*p2, (1-p1)*(1-p2))

# Create a matrix to hold the data
y <- matrix(NA, nrow=100, ncol=4)
dimnames(y) <- list(1:100, c("11", "10", "01", "00"))

# Loop over sites and generate data with function rmultinom()
for(i in 1:100){
   y[i,] <- rmultinom(1, N[i], cellprobs)
}

# Remove 4th column ("not detected") and summarize results
y <- y[,-4]
apply(y, 2, sum)
#11 10 01
#23 17  6


# Package the data up into an unmarkedFrame
library(unmarked)
umf <- unmarkedFrameMPois(y = y,
    type = "double")

# Fit models: multinomPois order of formulas: detection, abundance
fit <- multinomPois(~ 1 ~ 1, umf)

# Create an obsCov to fit observer effect....

#
# Part 2: look at some real data
#

library(unmarked)
data(ovendata)
ovendata.list$data[11:20,]   # Look at a snippet of data set
#      [,1] [,2] [,3] [,4]
# [1,]    0    0    0    0
# [2,]    0    0    0    0
# [3,]    2    0    1    0
# [4,]    1    0    0    0
# [5,]    1    0    0    0
# [6,]    0    0    0    0
# [7,]    0    1    0    1
# [8,]    0    1    0    0
# [9,]    0    0    0    0
#[10,]    0    0    0    0
 
apply(ovendata.list$data,2,sum)  # Removals in occasion 1-4
# [1] 49 16  5  7


# Package the data up into an unmarkedFrame

ovenFrame <- unmarkedFrameMPois(y = ovendata.list$data,
    siteCovs = as.data.frame(scale(ovendata.list$covariates[,-1])),
    type = "removal")

# Fit models: multinomPois order of formulas: detection, abundance
fm0 <- multinomPois(~ 1 ~ 1, ovenFrame)
fm1 <- multinomPois(~ 1 ~ ufc, ovenFrame)
fm2 <- multinomPois(~ 1 ~ trba, ovenFrame)
fm3 <- multinomPois(~ 1 ~ ufc + trba, ovenFrame)
fm4 <- multinomPois(~ 1 ~ ufc + trba + ufc:trba, ovenFrame)
fm5 <- multinomPois(~ ufc ~ ufc + trba, ovenFrame)
fm6 <- multinomPois(~ ufc ~ ufc + trba + ufc:trba, ovenFrame)

# Rank models by AIC
ms <- fitList(
"lam(.)p(.)"                                = fm0,
"lam(ufc)p(.)"                              = fm1,
"lam(trba)p(.)"                             = fm2,
"lam(ufc+trba)p(.)"                         = fm3,
"lam(ufc+trba+ufc:trba)p(.)"                = fm4,
"lam(ufc+trba)p(ufc)"                       = fm5,
"lam(ufc+trba+ufc:trba)p(ufc)"              = fm6)

(ms1 <- modSel(ms))
#                             nPars    AIC delta AICwt cumltvWt
#lam(trba)p(.)                    3 324.77  0.00 0.284     0.28
#lam(ufc)p(.)                     3 325.73  0.96 0.176     0.46
#lam(ufc+trba)p(.)                4 326.14  1.37 0.143     0.60
#lam(.)p(.)                       2 326.28  1.51 0.134     0.74
#lam(ufc+trba)p(ufc)              5 326.63  1.86 0.112     0.85
#lam(ufc+trba+ufc:trba)p(.)       5 327.17  2.40 0.086     0.93
#lam(ufc+trba+ufc:trba)p(ufc)     6 327.72  2.95 0.065     1.00
 

 
coef(ms1)[,1:4]  # Only first 4 columns shown
#                             lambda(Int) lambda(trba) lambda(ufc) lambda(ufc:trba)
#lam(trba)p(.)                  0.1062626   -0.2202827          NA               NA
#lam(ufc)p(.)                   0.1134562           NA   0.1789146               NA
#lam(ufc+trba)p(.)              0.1023681   -0.1708609   0.1002939               NA
#lam(.)p(.)                     0.1296349           NA          NA               NA
#lam(ufc+trba)p(ufc)            0.1069856   -0.1704169   0.1333443               NA
#lam(ufc+trba+ufc:trba)p(.)     0.1547003   -0.1499694   0.1563209        0.1298563
#lam(ufc+trba+ufc:trba)p(ufc)   0.1568600   -0.1503608   0.1854883        0.1264758


# Table with everything you could possibly need
 output <- as(ms1, "data.frame")
 
#
#
# Use of gmultmix function to do the same thing. NOTE: gmultmix fits a class of open population models so it has an 
#   extra formula called "phiformula" -- we discuss this later, so please ignore it for now
#
#

ovenFrame <- unmarkedFrameGMM(ovendata.list$data,
    siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
    numPrimary=1, #it could be pseudo open
    type = "removal")


fm0 <- gmultmix(lambdaformula = ~1, 
                phiformula = ~1, #models variation from one primary period and another
                pformula = ~1,
                data=ovenFrame)


# Fit Poisson models
fm1 <- gmultmix(~ ufc, ~ 1, ~  1, data = ovenFrame)
fm2 <- gmultmix(~ trba, ~ 1, ~ 1, data = ovenFrame)
fm3 <- gmultmix(~ ufc + trba, ~ 1, ~ 1, data = ovenFrame)
fm4 <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ 1, data = ovenFrame)
# Maybe p also depends on understory foliage?
fm5 <- gmultmix(~ ufc + trba, ~ 1, ~ ufc, data = ovenFrame)
fm6 <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ ufc, data = ovenFrame)

# Fit analogous NegBin models
fm0nb <- gmultmix(~ 1, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm1nb <- gmultmix(~ ufc, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm2nb <- gmultmix(~ trba, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm3nb <- gmultmix(~ ufc + trba , ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm4nb <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
# maybe p also depends on understory foliage?
fm5nb <- gmultmix(~ ufc + trba, ~ 1, ~ ufc, mixture = "NB", data = ovenFrame)
fm6nb <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ ufc, mixture = "NB",
data = ovenFrame)


# Rank models by AIC
gms <- fitList(
"lam(.)p(.)"                                = fm0,
"lam(ufc)p(.)"                              = fm1,
"lam(trba)p(.)"                             = fm2,
"lam(ufc+trba)p(.)"                         = fm3,
"lam(ufc+trba+ufc:trba)p(.)"                = fm4,
"lam(ufc+trba)p(ufc)"                       = fm5,
"lam(ufc+trba+ufc:trba)p(ufc)"              = fm6,
"NB,lam(.)p(.)"                             = fm0nb,
"NB,lam(ufc)p(.)"                           = fm1nb,
"NB,lam(trba)p(.)"                          = fm2nb,
"NB,lam(ufc+trba)p(.)"                      = fm3nb,
"NB,lam(ufc+trba+ufc:trba)p(.)"             = fm4nb,
"NB,lam(ufc+trba)p(ufc)"                    = fm5nb,
"NB,lam(ufc+trba+ufc:trba)p(ufc)"           = fm6nb)

(gms1 <- modSel(gms))

# Table with everything you could possibly need
output <- as(gms1, "data.frame")

#NB and poisson are in ∆AIC<2
fm2nb
#the dispersion parameter is 4.57
exp(4.57)

#
# Model fit assessment using parboot.
#
# need to create a function that computes fit statistics from the fit object in question
#

fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2,na.rm=TRUE)
    chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
    }

gof <- parboot(fm1, fitstats, nsim=100, report=1)

set.seed(2015)
(gof <- parboot(fm2, fitstats, nsim = 1000, report = 1))



#
#
# Part 2: Custom multinomial models , e.g., Capture-recapture models
#
#
#

#
# see help file ?piFuns
#

#
# piFun for a capture-recapture protocol
#

crPiFun <- function(p) {
  p1 <- p[,1] #detection probability fro historic count 1
  p2 <- p[,2] #detection probability for the second history count 
  p3 <- p[,3]
   cbind("001" = (1 - p1) * (1 - p2) *      p3,
         "010" = (1 - p1) *      p2  * (1 - p3),
         "011" = (1 - p1) *      p2  *      p3,
         "100" =      p1  * (1 - p2) * (1 - p3),
         "101" =      p1  * (1 - p2) *      p3,
         "110" =      p1  *      p2  * (1 - p3),
         "111" =      p1  *      p2  *      p3)
}


# 
# Chandler's ALFL data
#

alfl <- read.csv(system.file("csv", "alfl.csv", package="unmarked"))
head(alfl, 5)

#         id survey interval1 interval2 interval3
#1 crick1_05      1         1         1         1
#2 crick1_05      3         1         0         1
#3   his1_05      1         0         1         1
#4   his1_05      1         1         1         1
#5   his1_05      2         0         1         1


alfl.covs <- read.csv(system.file("csv", "alflCovs.csv",package="unmarked"),
 row.names=1)
head(alfl.covs)
#          struct woody time.1 time.2 time.3 date.1 date.2 date.3
#crick1_05   5.45  0.30   8.68   8.73   5.72      6     25     34
#his1_05     4.75  0.05   9.43   7.40   7.58     20     32     54
#hisw1_05   14.70  0.35   8.25   6.70   7.62     20     32     47
#hisw2_05    5.05  0.30   7.77   6.23   7.17     20     32     47
#kenc1_05    4.15  0.10   9.57   9.55   5.73      8     27     36
#kenc2_05    9.75  0.40   9.10   9.12   9.12      8     27     36


#
# Create encounter histories and then subset to the first survey of the season
#

alfl$captureHistory <- paste(alfl$interval1, alfl$interval2, alfl$interval3, sep="")
alfl$captureHistory <- factor(alfl$captureHistory,
    levels=c("001", "010", "011", "100", "101", "110", "111"))
alfl$id <- factor(alfl$id, levels=rownames(alfl.covs))

alfl.v1 <- alfl[alfl$survey==1,]
alfl.H1 <- table(alfl.v1$id, alfl.v1$captureHistory)
head(alfl.H1, 5)

#            001 010 011 100 101 110 111
#  crick1_05   0   0   0   0   0   0   1
#  his1_05     0   0   1   0   0   0   1
#  hisw1_05    0   0   0   0   0   0   0
#  hisw2_05    0   0   0   0   0   0   1
#  kenc1_05    0   0   0   0   0   0   0

#
## grab survey 2 and 3 here just in case we use them (later)
#
alfl.v2 <- alfl[alfl$survey==2,]
alfl.H2 <- table(alfl.v2$id, alfl.v2$captureHistory)

alfl.v3 <- alfl[alfl$survey==3,]
alfl.H3 <- table(alfl.v3$id, alfl.v3$captureHistory)

#
# Create a covariate "which interval" was the individual detected

#
intervalMat <- matrix(c('1','2','3'), 50, 3, byrow=TRUE)
class(alfl.H1) <- "matrix"
o2y <- matrix(1, 3, 7)

# Create the unmarkedFrame

umf.cr1 <- unmarkedFrameMPois(y=alfl.H1,
   siteCovs=alfl.covs[,c("woody", "struct", "time.1")],
   obsCovs=list(interval=intervalMat), 
   obsToY=o2y, 
   piFun="crPiFun")
# summary(umf.cr1)

#4 candidate models
M0 <- multinomPois(~ 1 ~ 1, umf.cr1)
Mt <- multinomPois(~ interval - 1 ~ 1, umf.cr1)
Mx <- multinomPois(~ time.1 ~ 1, umf.cr1)
(M0.woody <- multinomPois(~ 1 ~ woody, umf.cr1)) 
#more woody vegetation increase mean abundance of the flycatcher

# Put the models into a fitlist
fl <- modSel(fitList(M0, Mt, Mx, M0.woody))

# Create a plot of the effect of woody vegetation
nd <- data.frame(woody=seq(0, 0.8, length=50))
E.abundance <- predict(M0.woody, type="state", newdata=nd, appendData=TRUE)
plot(Predicted ~ woody, E.abundance, type="l", ylim=c(0, 6), ylab="Alder flycatchers / plot", xlab="Woody vegetation cover", frame = F)
lines(lower ~ woody, E.abundance, col=gray(0.7))
lines(upper ~ woody, E.abundance, col=gray(0.7))

##
##
## Class exercise: ALFL data
##      
##

#Fit NB abundance models to the flycatcher data using the gmultmix 
#function and compare to Poisson by AIC (also using gmultmix)

str(alfl.H1)

umf.cr2 <- unmarkedFrameGMM(y=alfl.H1,
                            siteCovs=alfl.covs[,c("woody", "struct", "time.1")],
                            numPrimary = 1,
                            obsCovs=list(interval=intervalMat), 
                            obsToY=o2y,
                            piFun="crPiFun")


# Fit Poisson models
fm0 <- gmultmix(lambdaformula = ~1, 
                phiformula = ~1, #models variation from one primary period and another
                pformula = ~1,
                data=umf.cr2)
fmt <- gmultmix(~ 1, ~ 1, ~  interval, data = umf.cr2)
fmx <- gmultmix(~ 1, ~ 1, ~ time.1, data = umf.cr2)
fm0.woody <- gmultmix(~ woody, ~ 1, ~ 1, data = umf.cr2)

# Fit analogous NegBin models
fm0nb <- gmultmix(~ 1, ~ 1, ~ 1, mixture = "NB", data = umf.cr2)
fmtnb <- gmultmix(~ 1, ~ 1, ~  interval, mixture = "NB", data = umf.cr2)
fmxnb <- gmultmix(~ 1, ~ 1, ~ time.1, mixture = "NB", data = umf.cr2)
fm0.woodynb <- gmultmix(~ woody, ~ 1, ~ 1, mixture = "NB", data = umf.cr2)


# Rank models by AIC
gms <- fitList(
  "Plam(.)p(.)"                                = fm0,
  "Plam(.)p(interval)"                              = fmt,
  "Plam(.)p(time)"                             = fmx,
  "Plam(woody)p(.)"                         = fm0.woody,
  "NB,lam(.)p(.)"              = fm0nb,
  "NB,lam(.)p(interval)"                     = fmtnb,
  "NB,lam(.)p(time)"            = fmxnb,
  "NB,lam(woody)p(.)"                              = fm0.woodynb)

(gms2 <- modSel(gms))

# Table with everything you could possibly need
output2 <- as(gms2, "data.frame")

#NB and poisson are in ∆AIC<2
fm0.woodynb
#the dispersion parameter is 8.43
exp(8.43)

#
# Model fit assessment using parboot.
#
# need to create a function that computes fit statistics from the fit object in question
#

fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2,na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

#Do a goodness-of-fit analysis of the NB and Poisson models

gof2 <- parboot(fm0.woody, fitstats, nsim=100, report=1)


# Create a plot of the effect of woody vegetation
nd <- data.frame(woody=seq(0, 0.8, length=50))
E.abundance <- predict(fm0.woody, type="lambda", newdata=nd, appendData=TRUE)
plot(Predicted ~ woody, E.abundance, type="l", ylim=c(0, 6), ylab="Alder flycatchers / plot", xlab="Woody vegetation cover", frame = F)
lines(lower ~ woody, E.abundance, col=gray(0.7))
lines(upper ~ woody, E.abundance, col=gray(0.7))


####################################################################################
#
#
#  PART 3: Bayesian analysis
#
#
#
####################################################################################

###
### Bayesian analysis of 3-part multinomial model for the ovenbird data
###
###

library(unmarked)
data(ovendata)
ovenFrame <- unmarkedFrameMPois(y = ovendata.list$data,
    siteCovs = as.data.frame(scale(ovendata.list$covariates[,-1])),
    type = "removal")


# Harvest the data and bundle it up for sending to BUGS
# y = removal counts (4 successive time periods -- "mental" removal)

y <- as.matrix(getY(ovenFrame))
ncap <- apply(y, 1, sum)   # number of individuals removed per point (detected)

# Bundle the data
data <- list(y = y, M = nrow(y), n = ncap, X=as.matrix(siteCovs(ovenFrame)))
str(data)                  # Good practice to always inspect your BUGS data

# Write BUGS model
cat("
model {

# Prior distributions
p0 ~ dunif(0,1) #
alpha0 <- logit(p0) #
alpha1 ~ dnorm(0, 0.01)
beta0 ~ dnorm(0, 0.01)
beta1 ~ dnorm(0, 0.01)
beta2 ~ dnorm(0, 0.01)
beta3 ~ dnorm(0, 0.01)

for(i in 1:M){ # Loop over sites
   # Conditional multinomial cell probabilities (four removal periods of 3-5 minutes)
   pi[i,1] <- p[i] #removed in the third period
   pi[i,2] <- p[i]*(1-p[i]) #miss in the first, removed in the second
   pi[i,3] <- p[i]*(1-p[i])*(1-p[i]) #miss in the first two, removed in the third
   pi[i,4] <- p[i]*(1-p[i])*(1-p[i])*(1-p[i]) #miss in the first three, removed in the fourth
   pi0[i] <- 1 - (pi[i,1] + pi[i,2] + pi[i,3] + pi[i,4]) #miss in all periods
   pcap[i] <- 1 - pi0[i]
   for(j in 1:4){
      pic[i,j] <- pi[i,j] / pcap[i] #conditional pi probabilities
   }

   # logit-linear model for detection: understory cover effect
   logit(p[i]) <- alpha0 + alpha1 * X[i,1] #regression model for the periods as function of understory

   # Model specification, three parts:
   y[i,1:4] ~ dmulti(pic[i,1:4], n[i]) # component 1 uses the conditional #multinomial
                                       #    cell probabilities
   n[i] ~ dbin(pcap[i], N[i])          # component 2 is a model for the #binomial
                                       #    observed sample size 
   N[i] ~ dpois(lambda[i])             # component 3 is the process model (Poisson)

   # log-linear model for abundance: UFC + TRBA + UFC:TRBA #model for expected pop size
 log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1]
}
}
",fill=TRUE, file="model.txt")

# Initial values
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), beta0 = runif(1), N = ncap+2)
}

# Parameters monitored
parameters <- c("p0", "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3")

# MCMC settings
nc <- 3   ;   ni <- 6000   ;   nb <- 1000   ;   nt <- 1

# Run JAGS and print posterior summary
library(jagsUI)
out1<- jags(data, inits, parameters, "model.txt", n.thin=nt, 
n.chains=nc, n.burnin=nb, n.iter=ni)

print(out1, 3)



#
#
# Add a GoF analysis here using posterior predictive check. See book (not covered in lecture)
#
#

# Write BUGS model
cat("
model {

# Prior distributions
p0 ~ dunif(0,1)
alpha0 <- logit(p0)
alpha1 ~ dnorm(0, 0.01)
beta0 ~ dnorm(0, 0.01)
beta1 ~ dnorm(0, 0.01)
beta2 ~ dnorm(0, 0.01)
beta3 ~ dnorm(0, 0.01)

for(i in 1:M){ # Loop over sites
   # Conditional multinomial cell probabilities
   pi[i,1] <- p[i]
   pi[i,2] <- p[i]*(1-p[i])
   pi[i,3] <- p[i]*(1-p[i])*(1-p[i])
   pi[i,4] <- p[i]*(1-p[i])*(1-p[i])*(1-p[i])
   pi0[i] <- 1 - (pi[i,1] + pi[i,2] + pi[i,3] + pi[i,4])
   pcap[i] <- 1 - pi0[i]
   for(j in 1:4){
      pic[i,j] <- pi[i,j] / pcap[i]
   }

   # logit-linear model for detection: understory cover effect
   logit(p[i]) <- alpha0 + alpha1 * X[i,1]

   # Model specification, three parts:
   y[i,1:4] ~ dmulti(pic[i,1:4], n[i]) # component 1 uses the conditional
                                       #    cell probabilities
   n[i] ~ dbin(pcap[i], N[i])          # component 2 is a model for the
                                       #    observed sample size 
   N[i] ~ dpois(lambda[i])             # component 3 is the process model

   # log-linear model for abundance: UFC + TRBA + UFC:TRBA
 log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1]
}

for(i in 1:M){
  # Simulate a new data set
   n.pred[i] ~ dbin(pcap[i], N[i])
   y.pred[i, 1:4] ~ dmulti(pic[i,1:4], n[i]) #note this is condl on ncap[i]

   for(k in 1:4){
      e1[i,k] <- pic[i,k] * n[i]
      resid1[i,k] <- pow(pow(y[i,k], 0.5) - pow(e1[i,k], 0.5), 2)
      resid1.pred[i,k] <- pow(pow(y.pred[i,k], 0.5) - pow(e1[i,k], 0.5), 2)
   }
   e2[i] <- pcap[i] * lambda[i]
   resid2[i] <- pow(pow(n[i], 0.5) - pow(e2[i], 0.5), 2)
   resid2.pred[i] <- pow(pow(n.pred[i], 0.5) - pow(e2[i], 0.5), 2)
}
fit1.data <- sum(resid1[,])         # Fit statistic for data y
fit1.pred <- sum(resid1.pred[,])

fit2.data<- sum(resid2[])           # Fit statistic for data n
fit2.pred<- sum(resid2.pred[])


}
",fill=TRUE, file="model2.txt")

# Initial values
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), beta0 = runif(1), N = ncap+2)
}


# MCMC settings
nc <- 3   ;   ni <- 6000   ;   nb <- 1000   ;   nt <- 1

# Parameters monitored

parameters <- c("N", "p0", "beta0", "beta1", "beta2", "beta3",
   "fit1.data", "fit1.pred", "fit2.data", "fit2.pred")

out2 <- jags (data, inits, parameters, "model2.txt",n.thin=nt, n.chains=nc, 
n.burnin=nb, n.iter=ni, parallel=TRUE)

print(out2, digits=3)
       
mean(out2$sims.list$fit1.pred > out2$sims.list$fit1.data)
#[1] 0.714

mean(out2$sims.list$fit2.pred > out2$sims.list$fit2.data)
#[1] 0.4066667

  
 

###
### Poisson formulation of the model
###


# Specify model in BUGS language
cat("
model {

# Prior distributions
p0 ~ dunif(0,1)
alpha0 <- logit(p0)
alpha1 ~ dnorm(0, 0.01)
beta0 ~ dnorm(0, 0.01)
beta1 ~ dnorm(0, 0.01)
beta2 ~ dnorm(0, 0.01)
beta3 ~ dnorm(0, 0.01)

for(i in 1:M){
   # logit-linear model for detection: understory cover effect
   logit(p[i]) <- alpha0 + alpha1 * X[i,1]
   # log-linear model for abundance: UFC + TRBA + UFC:TRBA
   log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1]

   # Poisson parameter = multinomial cellprobs x expected abundance (y_ih~Poisson(pi_ih,lambda_i))
   pi[i,1] <- p[i] * lambda[i]
   pi[i,2] <- p[i] * (1-p[i]) * lambda[i]
   pi[i,3] <- p[i] * (1-p[i]) * (1-p[i]) * lambda[i]
   pi[i,4] <- p[i] * (1-p[i]) * (1-p[i]) * (1-p[i]) * lambda[i]

   for(j in 1:4){
      y[i,j] ~ dpois(pi[i,j])
   }
   # Generate predictions of N[i]
   N[i] ~ dpois(lambda[i])
}
}
",fill=TRUE,file="modelP.txt")

# Bundle up the data and inits
data <- list(y = y, M = nrow(y), X = as.matrix(siteCovs(ovenFrame)))
inits <- function(){
  list (p0 = runif(1), alpha1=runif(1), beta0=runif(1), beta1=runif(1), beta2=runif(1), 
beta3=runif(1))
}

# Define parameters to save and MCMC settings
parameters <- c("p0", "alpha1", "beta0", "beta1", "beta2", "beta3" )
nc <- 3   ;   ni <- 6000   ;   nb <- 1000   ;   nt <- 1

 

# Call WinBUGS from R and summarize marginal posteriors 
out4 <- jags(data, inits, parameters, "modelP.txt", n.thin=nt,
   n.chains = nc, n.burnin = nb, n.iter = ni,parallel=FALSE)
print(out4, 3)

#########################################################################
#define categorical for site (1 for sites 1:35, 2 for sites 36:70)


library(unmarked)
data(ovendata)

#Generate a new category 1-2 for the sites

cat <- c(rep(1,35), rep(2,35))

ovenFrame <- unmarkedFrameMPois(y = ovendata.list$data,
                                siteCovs = as.data.frame(scale(ovendata.list$covariates[,-1])),
                                type = "removal")

# Harvest the data and bundle it up for sending to BUGS
# y = removal counts (4 successive time periods -- "mental" removal)

y <- as.matrix(getY(ovenFrame))
ncap <- apply(y, 1, sum)   # number of individuals removed per point (detected)

# Bundle the data
data <- list(y = y, M = nrow(y), n = ncap, #Andy Royle put a new object here, cat = c(rep(1,35), rep(2,35))
             cat = cat, 
             X=as.matrix(siteCovs(ovenFrame)))
str(data)                  # Good practice to always inspect your BUGS data

# Write BUGS model
cat("
model {

# Prior distributions
p0 ~ dunif(0,1) #
alpha0 <- logit(p0) #
alpha1 ~ dnorm(0, 0.01)
beta0[1] ~ dnorm(0, 0.01)
beta0[2] ~ dnorm(0, 0.01)

beta1 ~ dnorm(0, 0.01)
beta2 ~ dnorm(0, 0.01)
beta3 ~ dnorm(0, 0.01)
beta0.diff <- beta0[1] - beta0[2]

for(i in 1:M){ # Loop over sites
   # Conditional multinomial cell probabilities (four removal periods of 3-5 minutes)
   pi[i,1] <- p[i] #removed in the third period
   pi[i,2] <- p[i]*(1-p[i]) #miss in the first, removed in the second
   pi[i,3] <- p[i]*(1-p[i])*(1-p[i]) #miss in the first two, removed in the third
   pi[i,4] <- p[i]*(1-p[i])*(1-p[i])*(1-p[i]) #miss in the first three, removed in the fourth
   pi0[i] <- 1 - (pi[i,1] + pi[i,2] + pi[i,3] + pi[i,4]) #miss in all periods
   pcap[i] <- 1 - pi0[i]
   for(j in 1:4){
      pic[i,j] <- pi[i,j] / pcap[i] #conditional pi probabilities
   }

   # logit-linear model for detection: understory cover effect
   logit(p[i]) <- alpha0 + alpha1 * X[i,1] #regression model for the periods as function of understory

   # Model specification, three parts:
   y[i,1:4] ~ dmulti(pic[i,1:4], n[i]) # component 1 uses the conditional #multinomial
                                       #    cell probabilities
   n[i] ~ dbin(pcap[i], N[i])          # component 2 is a model for the #binomial
                                       #    observed sample size 
   N[i] ~ dpois(lambda[i])             # component 3 is the process model (Poisson)

   # log-linear model for abundance: UFC + TRBA + UFC:TRBA #model for expected pop size
 log(lambda[i])<- beta0[cat[i]] + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1]
}
}
",fill=TRUE, file="model_SiteCat.txt")

# Initial values
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), N = ncap+2)
}

# Parameters monitored
parameters <- c("p0", "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3")

# MCMC settings
nc <- 3   ;   ni <- 6000   ;   nb <- 1000   ;   nt <- 1

# Run JAGS and print posterior summary
library(jagsUI)
out_SiteCat<- jags(data, inits, parameters, "model_SiteCat.txt", n.thin=nt, 
            n.chains=nc, n.burnin=nb, n.iter=ni)

print(out_SiteCat, 3)

#########################################################################
#
#
#  PART 4: Open populations (multi-year data)
#
#
########################################################################

##########################################################################
#
# Analysis of the ALFL data using open models 
# Next we create the unmarkedFrame from the ALFL data
#
#########################################################################

crPiFun <- function(p) {
  p1 <- p[,1]
  p2 <- p[,2]
  p3 <- p[,3]
   cbind("001" = (1 - p1) * (1 - p2) *      p3,
         "010" = (1 - p1) *      p2  * (1 - p3),
         "011" = (1 - p1) *      p2  *      p3,
         "100" =      p1  * (1 - p2) * (1 - p3),
         "101" =      p1  * (1 - p2) *      p3,
         "110" =      p1  *      p2  * (1 - p3),
         "111" =      p1  *      p2  *      p3)
}

alfl <- read.csv(system.file("csv", "alfl.csv", package = "unmarked"))
alfl.covs <- read.csv(system.file("csv", "alflCovs.csv", package = "unmarked"), row.names = 1)
head(alfl.covs)

alfl$captureHistory <- paste(alfl$interval1, alfl$interval2,alfl$interval3, sep = "")
alfl$captureHistory <- factor(alfl$captureHistory,
                              levels = c("001", "010", "011", "100", "101", "110", "111"))
alfl$id <- factor(alfl$id, levels = rownames(alfl.covs))
head(alfl, 5)

# Standardize some obsCovs
time <- as.matrix(alfl.covs[,c("time.1", "time.2", "time.3")] )
time <- time - median(time)

date <- as.matrix(alfl.covs[,c("date.1", "date.2", "date.3")] )
date <- date - median(date)

# Note: we analyzed survey == 1 in ch. 7 in AHM1, here use surveys for 3 days within the same breeding
#   season. For each day we assumed closure, but perhaps not across days
#
alfl.v1 <- alfl[alfl$survey == 1,]
alfl.H1 <- table(alfl.v1$id, alfl.v1$captureHistory)
alfl.v2 <- alfl[alfl$survey == 2,]
alfl.H2 <- table(alfl.v2$id, alfl.v2$captureHistory)
alfl.v3 <- alfl[alfl$survey == 3,]
alfl.H3 <- table(alfl.v3$id, alfl.v3$captureHistory)

str(alfl.v1)
str(alfl.H1)

# Arrange the data in wide format
Ywide <- cbind(alfl.H1, alfl.H2, alfl.H3)
str(Ywide)

# unmarkedFrame for a static model just using 1 primary sample
intervalMat <- matrix(c("1", "2", "3"), 50, 3, byrow = TRUE)
class(alfl.H1) <- "matrix"
o2y <- matrix(1, 3, 7)  # This is the correct obsToY matrix for capture-recapture. Trust me.

# SESSION variables -- Occasion is index to primary... should be YearlySiteCov but using this to illustrate
occ <- matrix(NA, nrow = 50, ncol = 9)
occ[, 1:3]<- 1
occ[, 4:6]<- 2
occ[, 7:9]<- 3

#
#
# First we ignore temporal variability and fit some closed models to the stacked data
#
#

# Stack the data
Ystacked <- rbind(Ywide[,1:7], Ywide[,8:14], Ywide[,15:21])
sc <- alfl.covs[,c("woody", "struct")]

# Stack the siteCovs
sc <- rbind(sc, sc, sc) # Repeat the site covs matrix
sc <- cbind(sc, date = c(date), time = c(time)) # YearlySiteCovs become siteCovs

# Specify numPrimary = 1 for the stacked data
summary(stacked.umf <- unmarkedFrameGMM(y = Ystacked, siteCovs = sc, numPrimary = 1,
                                        obsToY = o2y, piFun = "crPiFun") )

#
# Use gmultmix here. Remember order of formulas: lambda, phi, p. phi = ~1 for stacked data
# For stacked data, covariates across primary sessions could be on p or lambda and it's not usually clear where to put them.
# We demonstrate that here using models with time of day on lambda or p or both
#
m0 <- gmultmix( ~1, ~1, ~1, data = stacked.umf, mixture = "P")
m1 <- gmultmix( ~1, ~1, ~time, data = stacked.umf, mixture = "P")
m1b <- gmultmix( ~time, ~1, ~1, data = stacked.umf, mixture = "P")
m1c <- gmultmix( ~time, ~1, ~time, data = stacked.umf, mixture = "P")

# Models with date on lambda or p or both
m2 <- gmultmix( ~1, ~1, ~date, data = stacked.umf, mixture = "P")
m2b <- gmultmix( ~date, ~1, ~1, data = stacked.umf, mixture = "P")
m2c <- gmultmix( ~date, ~1, ~date, data = stacked.umf, mixture = "P")

# Models with both time and date .... all run in no time
m3 <- gmultmix( ~1, ~1, ~time + date, data = stacked.umf, mixture = "P")

# Add woody and struct to the best model
m4 <- gmultmix( ~woody + struct + date, ~1, ~date, data = stacked.umf, mixture = "P")

fl<- fitList(m0, m1, m1b, m1c, m2, m2b, m2c, m3)
modSel(fl)

# What does it mean? It is not clear.

##################################################################################
#
# Next we fit slightly open models (TE) where N(t) is random for each replicate survey
#
##################################################################################


summary(alfl.umf <- unmarkedFrameGMM(y = Ywide,
                                     siteCovs = alfl.covs[,c("woody", "struct")], 
                                     obsCovs = list(occ = occ),
                                     numPrimary = 3, #key to tell the sesions!!!
                                     yearlySiteCovs = list(time = time, date = date),
                                     obsToY = o2y, piFun = "crPiFun") )

# 
# Fit a suite of models with detection covariates on p or phi (availability) or both
# ORDER OF FORMULAS:  lambdaformula (abundance), phiformula (availability), pformula (detection), 
# Here we think of things like time of day as affecting availability most likely. Given availability, detection should not
# depend on time of day. But this is debateable!
#

m0 <- gmultmix( ~1, ~1, ~1, data = alfl.umf, mixture = "P")
m1 <- gmultmix( ~1, ~time, ~1, data = alfl.umf, mixture = "P")
m1b <- gmultmix( ~1, ~time, ~time, data = alfl.umf, mixture = "P")
m1c <- gmultmix( ~1, ~1, ~time, data = alfl.umf, mixture = "P")
m2 <- gmultmix( ~1, ~date, ~1, data = alfl.umf, mixture = "P")
m2b <- gmultmix( ~1, ~date, ~date, data = alfl.umf, mixture = "P")
m2c <- gmultmix( ~1, ~1, ~date, data = alfl.umf, mixture = "P")
m3 <- gmultmix( ~1, ~time + date, ~1, data = alfl.umf, mixture = "P")
m3b <- gmultmix( ~1, ~time + date, ~time, data = alfl.umf, mixture = "P")
m3c <- gmultmix( ~1, ~date, ~time, data = alfl.umf, mixture = "P")
m3d <- gmultmix( ~1, ~time, ~date, data = alfl.umf, mixture = "P")
m3e <- gmultmix( ~1, ~time, ~time + date, data = alfl.umf, mixture = "P")
m3f <- gmultmix( ~1, ~1, ~time + date, data = alfl.umf, mixture = "P")
m3g <- gmultmix( ~1, ~time + date, ~date, data = alfl.umf, mixture = "P")
m3h <- gmultmix( ~1, ~time + date, ~time + date, data = alfl.umf, mixture = "P")
m3i <- gmultmix( ~1, ~date, ~time + date, data = alfl.umf, mixture = "P")


# We organize these model fits into a fitList and then take a look at the model selection table:
(fl <- modSel(fitList(m0, m1, m1b, m1c, m2, m2b, m2c, m3, m3b, m3c, m3d, m3e, m3f, m3g, m3h, m3i)) )

#best: 
m3g
#time and date affect availability, date affect detection probability (but not a big effect)

# Now I'm going to add habitat to the abundance model
m4 <- gmultmix( ~woody + struct, ~time + date, ~date, data = alfl.umf, mixture = "P")
m4
#positive effect of woody vegetation, slight effect of understory structure

# Carry out a parametric bootstrap goodness-of-fit test (ART 7 min)
(pb <- parboot(m3g, fitstats, nsim = 30, report = 1, ncores = 6))
 


#
#
# Now we could try the new D-M model!  But note: we cannot yet handle custom piFuns so we have to convert the data to removal....
#   (IN PROGRESS)
#

# We have to make removal data or 2-session data from these data b/c multimixOpen does not handle more general cases yet. To make
# removal data we have to combine the encounter histories to pool those where 1st capture was in period 1, those where 1st capture
# was in period 2, and those where 1st capture was in period 3. Like this:

y<- Ywide
y.rem<- matrix(NA,nrow=50, ncol= 9)
y.rem[1:50, 1]<-  rowSums(y[1:50,4:7])     # Session 1 all encounter histories where first capture was in sample period 1
y.rem[1:50, 4] <- rowSums( y[1:50,11:14])  # Session 2 all encounter histories where first capture was in sample period 2
y.rem[1:50, 7]<-  rowSums( y[1:50,18:21] ) 
y.rem[1:50, 2]<-  rowSums( y[1:50, 2:3] )
y.rem[1:50, 5]<-  rowSums( y[1:50, 9:10]) 
y.rem[1:50, 8]<-  rowSums( y[1:50, 16:17])
y.rem[1:50, 3]<-    ( y[1:50, 1] )  
y.rem[1:50, 6]<-    ( y[1:50, 8] )
y.rem[1:50, 9]<-    ( y[1:50, 15] )


# Looks right:
apply(y.rem,2,sum)

##   [1] 40  8  2 24  5  2 11  6  0

# Compare and contrast MMO vs GMM 

umf1 <- unmarkedFrameMMO(y = y.rem, numPrimary=3, 
                         siteCovs=alfl.covs[,c("woody", "struct", "time.1")], 
                         type="removal")

umf2<- unmarkedFrameGMM(y = y.rem, numPrimary=3, 
                        siteCovs=alfl.covs[,c("woody", "struct", "time.1")], 
                        type="removal")   # Same as before

M0.1 <- multmixOpen(~ woody + struct , ~1, ~1, ~ 1, umf1, K = 100)
M0.2 <- gmultmix(~woody + struct, ~1, ~ 1, umf2, K = 100)

# Temporary emigration model for the removal data of the flycatcher data
M0.2

#Full dynamic model
M0.1

#intercept of abundance is lower in full dynamic

#We can model migration dynamics!!!!!

#####################################################################
#
#
# PART 5: Bayesian analysis of Open population models for multi-session data: Using ALFL data.
#
# Previously we used unmarked. We input and structure the data in the same manner. 
#
#####################################################################


alfl <- read.csv(system.file("csv", "alfl.csv", package = "unmarked"))
alfl.covs <- read.csv(system.file("csv", "alflCovs.csv", package = "unmarked"), row.names = 1)
head(alfl.covs)

alfl$captureHistory <- paste(alfl$interval1, alfl$interval2,
alfl$interval3, sep = "")
alfl$captureHistory <- factor(alfl$captureHistory,
levels = c("001", "010", "011", "100", "101", "110", "111"))
alfl$id <- factor(alfl$id, levels = rownames(alfl.covs))
head(alfl, 5)

# Standardize some obsCovs
time <- as.matrix(alfl.covs[,c("time.1", "time.2", "time.3")] )
time <- time - median(time)

date <- as.matrix(alfl.covs[,c("date.1", "date.2", "date.3")] )
date <- date - median(date)

# Note: we analyzed survey == 1 in ch. 7 in AHM1, here use all surveys
alfl.v1 <- alfl[alfl$survey == 1,]
alfl.H1 <- table(alfl.v1$id, alfl.v1$captureHistory)
alfl.v2 <- alfl[alfl$survey == 2,]
alfl.H2 <- table(alfl.v2$id, alfl.v2$captureHistory)
alfl.v3 <- alfl[alfl$survey == 3,]
alfl.H3 <- table(alfl.v3$id, alfl.v3$captureHistory)
 
Ywide <- cbind(alfl.H1, alfl.H2, alfl.H3)

# unmarkedFrame for a static model just using 1 primary sample
intervalMat <- matrix(c("1", "2", "3"), 50, 3, byrow = TRUE)
occ <- matrix(NA, nrow = 50, ncol = 9)
occ <- col(occ)


crPiFun <- function(p) {
  p1 <- p[,1]
  p2 <- p[,2]
  p3 <- p[,3]
   cbind("001" = (1 - p1) * (1 - p2) *      p3,
         "010" = (1 - p1) *      p2  * (1 - p3),
         "011" = (1 - p1) *      p2  *      p3,
         "100" =      p1  * (1 - p2) * (1 - p3),
         "101" =      p1  * (1 - p2) *      p3,
         "110" =      p1  *      p2  * (1 - p3),
         "111" =      p1  *      p2  *      p3)
}


#
# Now we set things up to run the temporary emigration model in BUGS 
# 

#
# We use a 3-d data structure which is very convenient in BUGS model specifications
#
y3d <- array(NA, dim = c(nrow(Ywide), 7, 3) ) # Create 3D array and fill it
y3d[,,1] <- Ywide[,1:7]
y3d[,,2] <- Ywide[,8:14]
y3d[,,3] <- Ywide[,15:21]
nsessions <- 3 # Number of primary sessions
nsites <- nrow(Ywide) # Number of sites
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and sessions

 
# Bundle the data
str(bdata <- list(y3d = y3d, 
                  nsites = nsites, 
                  nsessions = nsessions,
                  nobs = nobs, woody = alfl.covs[,"woody"], 
                  struct = alfl.covs[,"struct"],
                  date = date, time = time, pi = pi))
 
# Specify model in BUGS language
cat(file="CR_TE.txt", "
model {
# Priors
# Abundance parameters
beta0 ~ dnorm(0, 0.01)
beta1 ~ dnorm(0, 0.01)
beta2 ~ dnorm(0, 0.01)
# Availability parameters
gamma0 ~ dnorm(0, 0.01)
gamma1 ~ dnorm(0, 0.01)
gamma2 ~ dnorm(0, 0.01)
# Detection parameters
alpha0 ~ dnorm(0, 0.01)
alpha1 ~ dnorm(0, 0.01)

# Likelihood
for (i in 1:nsites) {
  # Linear model for lambda
  log(lambda[i]) <- beta0 + beta1*woody[i] + beta2*struct[i]

  for(k in 1:nsessions){
      # Linear models for phi and p
      logit(phi[i,k]) <- gamma0 + gamma1*date[i,k] + gamma2*time[i,k]
      logit(p[i,k]) <- alpha0 + alpha1*date[i,k]
      # Define multinomial cell probabilities
        cp[i,k,1] <- (1-p[i,k])*(1-p[i,k])*p[i,k]
        cp[i,k,2] <- (1-p[i,k])*p[i,k]*(1-p[i,k])
        cp[i,k,3] <- (1-p[i,k])*p[i,k]*p[i,k]
        cp[i,k,4] <- p[i,k]*(1-p[i,k])*(1-p[i,k])
        cp[i,k,5] <- p[i,k]*(1-p[i,k])*p[i,k]
        cp[i,k,6] <- p[i,k]*p[i,k]*(1-p[i,k])
        cp[i,k,7] <- p[i,k]*p[i,k]*p[i,k]
        cp[i,k,8] <- 1-sum(cp[i,k,1:7])
        cellprobs.cond[i,k,1] <- cp[i,k,1]/sum(cp[i,k,1:7])
        cellprobs.cond[i,k,2] <- cp[i,k,2]/sum(cp[i,k,1:7])
        cellprobs.cond[i,k,3] <- cp[i,k,3]/sum(cp[i,k,1:7])
        cellprobs.cond[i,k,4] <- cp[i,k,4]/sum(cp[i,k,1:7])
        cellprobs.cond[i,k,5] <- cp[i,k,5]/sum(cp[i,k,1:7])
        cellprobs.cond[i,k,6] <- cp[i,k,6]/sum(cp[i,k,1:7])
        cellprobs.cond[i,k,7] <- cp[i,k,7]/sum(cp[i,k,1:7])

# Conditional 4-part version of the model
  pdet[i,k] <- sum(cp[i,k, 1:7])
  pmarg[i,k] <- pdet[i,k]*phi[i,k]
#Part 4: Multinomial
  y3d[i,1:7,k] ~ dmulti(cellprobs.cond[i,k,1:7], nobs[i,k])
#Part 3: Number of detected individuals
  nobs[i,k] ~ dbin(pmarg[i,k], M[i]) #n individuals observed conditional to the super population M with marginal detection
#Part 2: Number of available individuals
  Navail[i,k] ~ dbin(phi[i,k], M[i])
} # end k loop

M[i] ~ dpois(lambda[i]) #Part 1: Abundance model
} # end i loop

# Derived quantities
for(k in 1:nsessions){
  Davail[k] <- mean(phi[,k])*exp(beta0)/area
}
Mtotal <- sum(M[])
area <- (pi*50*50)/10000 # 50 m point counts, D in units per ha
Dtotal <- exp(beta0)/area
}
")


# Initial values
Navail.st <- apply(y3d, c(1, 3), sum)
Mst <- apply(Navail.st, 1, max, na.rm = TRUE) + 2
inits <- function() list(M = Mst, sigma = 100)

# Parameters monitored
params <- c("beta0" , "beta1", "beta2", "alpha0", "alpha1", "gamma0",
"gamma1", "gamma2", "Mtotal", "Davail", "Dtotal")

# MCMC settings
na <- 1000 ; ni <- 10000 ; nb <- 5000 ; nt <- 5 ; nc <- 3

# Call JAGS (ART 0.2 min), gauge convergence, summarize posteriors
out7 <- jags(bdata, inits, params, "CR_TE.txt", n.adapt = na, n.iter = ni,
n.burnin = nb, n.thin = nt, n.chains = nc, parallel = TRUE)
 
print(out7, 3)

