#' ---
#' title: An R Markdown document converted from "occupancy-abundance-unmarked.ipynb"
#' output: html_document
#' ---
#' 
#' # Occupancy-abundance models in `unmarked`
#' 
#' ### Ken Kellner
#' 
#' ## Introduction
#' 
#' This model is in development and so is not available in `unmarked` from CRAN.
#' 
## -----------------------------------------------------------------------------
# remotes::install_github("rbchan/unmarked", ref="occuRNMulti")
library(unmarked)

#' 
#' We will use the function `occuRNMulti`*
#' 
#' \* *Name will almost certainly change*
#' 
#' ![](occuRNMulti.png)
#' 
#' ## Model support
#' 
#' Supports up to 3 species, and most possible combinations:
#' 
#' sp1 --> sp2
#' 
#' sp1 --> sp2 --> sp3
#' 
#' sp1 --> sp3 <-- sp2
#' 
#' sp2 <-- sp1 --> sp3
#' 
#' sp1 --> sp2 --> sp3 & sp1 --> sp3
#' 
#' ## Set up the input data
#' 
#' Mesocarnivores again.
#' 
## -----------------------------------------------------------------------------
data(MesoCarnivores)
lapply(MesoCarnivores, head) # look at raw data

#' 
#' ### Formatting the data
#' 
#' We will create an `unmarkedFrameOccuMulti` exactly as we did in the multispecies occupancy section.
#' 
#' ![](umf_om_construction.png)
#' 
#' ### y list
#' 
#' We'll run a two-species model
#' 
## -----------------------------------------------------------------------------
ylist <- MesoCarnivores[c(2:3)]
lapply(ylist, head) # look at first few rows

#' 
#' ### Site covariates
#' 
#' The last element of the `MesoCarnivores` list is a data frame of site covariates:
#' 
## -----------------------------------------------------------------------------
site_covs <- MesoCarnivores$sitecovs
head(site_covs)

#' 
#' ### Construct the unmarked frame
#' 
## -----------------------------------------------------------------------------
umf <- unmarkedFrameOccuMulti(y=ylist, siteCovs=site_covs)
head(umf)

#' 
## -----------------------------------------------------------------------------
summary(umf)

#' 
#' ## Set up the formulas
#' 
#' ![](occuRNMulti.png)
#' 
#' ### `stateformulas`
#' 
#' A named list with one element per species.
#' 
#' The elements should be in order of dominance, starting with the top-level dominant species.
#' 
#' **For dominant species:**
#' 
#' The list element contains a single formula for the abundance model.
#' 
#' ```r
#' coyote = ~1
#' ```
#' 
#' $$\mathrm{log}(\lambda_{coy}) = \beta_{0,coy}$$
#' 
#' 
#' **For subordinate species:**
#' 
#' The list element is itself a list of multiple formulas.
#' 
#' * One formula for the part of the linear predictor not related to the dominant species
#' * One formula representing the part of the linear predictor corresponding to each dominant species
#' 
#' ```r
#' redfox = list(~1, coyote = ~1)
#' ```
#' 
#' $$\mathrm{log}(\lambda_{fox}) = \beta_{0,fox} + (\gamma_{0,coy->fox}) \cdot N_{coy}$$
#' 
#' With covariates:
#' 
#' ```r
#' redfox = list(~x, coyote = ~z)
#' ```
#' 
#' $$\mathrm{log}(\lambda_{fox}) = \beta_{0,fox} + \beta_{x,fox} + (\gamma_{0,coy->fox} + \gamma_{z,coy->fox} \cdot z) \cdot N_{coy}$$
#' 
## -----------------------------------------------------------------------------
sf <- list(coyote = ~1,
           redfox = list(~1, coyote = ~1))

#' 
#' ### `detformulas`
#' 
#' This one is easier, there is just one formula per species, so there should be 2 total.
#' 
#' Intercept-only models for both species:
#' 
## -----------------------------------------------------------------------------
df <- list(coyote = ~1, redfox = ~1)
df

#' 
#' ### `modelOccupancy`
#' 
#' A character vector of species names for which you want to model occupancy rather than abundance (must be subordinate species)
#' 
#' ```r
#' modelOccupancy = 'redfox'
#' ```
#' 
#' ## Run `occuRNMulti`
#' 
#' We now have all the pieces we need (`unmarkedFrameOccuMulti`, `stateformulas`, `detformulas`) needed to run a model.
#' 
## -----------------------------------------------------------------------------
mod_null <- occuRNMulti(detformulas=df, stateformulas=sf, data=umf, modelOccupancy="redfox")
summary(mod_null)

#' 
#' ## Model with covariates
#' 
#' * Effects of disturbance, housing density, longitude, and latitude on both species
#' * Effect of housing density on coyote -> red fox interaction
#' * Effect of trail on both detection models
#' 
## -----------------------------------------------------------------------------
both <- ~scale(Dist_5km) + scale(Longitude) + scale(Latitude) + scale(HDens_5km)

sf <- list(
  coyote = both,
  redfox = list(both, coyote = ~scale(HDens_5km))
)

#' 
## -----------------------------------------------------------------------------
df <- list(coyote = ~Trail, redfox = ~Trail)

#' 
## -----------------------------------------------------------------------------
mod_covs <- occuRNMulti(detformulas=df, stateformulas=sf, data=umf, modelOccupancy="redfox", threads=3)
summary(mod_covs)

#' 
#' ### Comparison with `NIMBLE`
#' 
#' Read in saved results from the same model fit with `NIMBLE`:
#' 
## -----------------------------------------------------------------------------
(nim_sum <- as.data.frame(readRDS('nimble_summary.Rds')))

#' 
## -----------------------------------------------------------------------------
nim_sum$type <- "nimble"

unm_sum <- data.frame(est = coef(mod_covs),
                      low = coef(mod_covs) - 1.96 * SE(mod_covs),
                      up = coef(mod_covs) + 1.96 * SE(mod_covs),
                      type = "unmarked")

all_sum <- rbind(nim_sum, unm_sum)
all_sum$par <- rownames(nim_sum)

#' 
## -----------------------------------------------------------------------------
library(ggplot2)
options(repr.plot.width=15, repr.plot.height=10)
pos <- position_dodge(0.5)

ggplot(data = all_sum, aes(x=par, col=type)) +
  geom_errorbar(aes(ymin=low, ymax=up), width=0.2, position=pos) +
  geom_point(aes(y=est), position=pos) +
  theme_bw(base_size=16) +
  labs(col=element_blank(), x=element_blank(), y="Estimate and 95% CI") +
  theme(panel.grid = element_blank(), 
        legend.position.inside=c(0.8,0.2), legend.text=element_text(size=20))

#' 
#' ## Plotting effect of `Hdens_5km`
#' 
#' Make a sequence of values and plug into a `newdata` data frame.
#' 
## -----------------------------------------------------------------------------
hdens_rng <- range(siteCovs(umf)$HDens_5km)
hdens_seq <- seq(hdens_rng[1], hdens_rng[2], length.out=100)

nd <- data.frame(HDens_5km = hdens_seq, 
                 Dist_5km = mean(siteCovs(umf)$Dist_5km),
                 Longitude = mean(siteCovs(umf)$Longitude),
                 Latitude = mean(siteCovs(umf)$Latitude),
                 Trail = 0)

head(nd)

#' 
#' Run `predict`. It generates a list of dataframes, one per species.
#' 
## -----------------------------------------------------------------------------
pr <- predict(mod_covs, type = "state", newdata = nd)
lapply(pr, tail)

#' 
#' Note that the predictions will match the type of model for each species.
#' 
#' So `coyote` represents abundances, and `redfox` represents occupancy, since we set `modelOccupancy = "redfox"`.
#' 
#' ### Create plot
#' 
#' First combine the outputs into a single data frame for `ggplot2`.
#' 
## -----------------------------------------------------------------------------
plot_data <- do.call(rbind, pr)
plot_data$species <- rep(c("coyote", "redfox"), each=100)
plot_data$Hdens <- hdens_seq
head(plot_data)

#' 
## -----------------------------------------------------------------------------
library(ggplot2)

ggplot(data = plot_data, aes(x = Hdens)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=species), alpha=0.2) +
  geom_line(aes(y = Predicted, col=species)) +
  facet_wrap("species", scales='free') +
  theme_bw(base_size=20) +
  theme(panel.grid = element_blank())

#' 
#' The huge error bars for red fox are likely due to the high uncertainty in the interaction parameter estimates.
#' 
#' ### Fix abundance of a species
#' 
#' We might want to get predictions for a specific abundance value of a dominant species.
#' 
#' For example here, with coyote absent (abundance = 0)
#' 
#' We just need to add the species abundance to the `newdata` data frame:
#' 
## -----------------------------------------------------------------------------
nd$coyote <- 0
pr <- predict(mod_covs, type = "state", newdata = nd)
lapply(pr, tail)

#' 
#' Note that while still high, the SEs around the occupancy estimates for red fox are much smaller than before. That's because we are "zeroing out" the interaction terms and their associated uncertainty.
#' 
## -----------------------------------------------------------------------------
plot_data <- do.call(rbind, pr)
plot_data$species <- rep(c("coyote", "redfox"), each=100)
plot_data$Hdens <- hdens_seq
ggplot(data = plot_data, aes(x = Hdens)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=species), alpha=0.2) +
  geom_line(aes(y = Predicted, col=species)) +
  facet_wrap("species", scales='free') +
  theme_bw(base_size=20) +
  theme(panel.grid = element_blank())

#' 
#' Overall red fox predicted occupancy is lower when coyote is absent.
#' 
#' That makes sense given the positive interaction intercept in the model.
#' 
#' If you want to try `occuRNMulti` get in touch or install yourself using the code I provided at the beginning.
#' 
#' There are definitely missing features and probably bugs!
#' 
