##########################################################
### Setup
##########################################################
rm(list = ls())

library(dplyr)
library(parallel)
library(MatchIt)
library(gdata)
library(survival)
library(marginaleffects)
library(lme4)
library(lmtest) #coeftest
library(sandwich) #vcovCL

source("~/Matching/Matching-Applied/Scripts/match_functions.R")
source("~/Matching/Matching-Applied/Scripts/simul.R")
source("~/Matching/Matching-Applied/Scripts/mod_definitions.R")

dir.create("../Results/", showWarnings = F)

args = commandArgs(trailingOnly=TRUE)
map_idx = args[1]
#map_idx = 181

############################
# SPECIFICATIONS #
############################
# One variable - linear
seeds = 1:100

### simulation specs
mod=2
mod_specs <- get_mod_specs(mod)

### matching specs
dat="data"
predictor <- "treat"
confounder <- "var1"
f1 <- paste0(predictor, " ~ ", paste0(confounder, collapse = "+"))
confounder <- paste0("var", 1:mod_specs$pk)
f2 <- paste0(predictor, " ~ ", paste0(confounder, collapse = "+"))
match.formulas <- c(f2)
estimand = c("ATC", "ATT", "ATE")
methods=c("null","nearest", "quick", "optimal", "genetic", "full")
distances=c("glm", "gam", "randomforest", "cbps", "nnet", "bart", "robust_mahalanobis", "scaled_euclidean")
ratio=c(1)
#caliper
#replace

### estimation specs
outcome <- "y"
f1 <- paste0(outcome, " ~ ", predictor)
f2 <- paste0(outcome, " ~ ", predictor, "+", paste0(confounder, collapse = "+"))
stat.formulas <- c(f1, f2)
family = "gaussian"
covar.names <- paste0("c(", paste0("\"", confounder, "\"", collapse = ",") ,")")

### expand all and set seed
dlist <- append(mod_specs,
                list(
                  dat=dat, f=match.formulas, estimand=estimand, method=methods, distance=distances, ratio = ratio,
                  stat.formula = stat.formulas, family=family, covar.names = covar.names, stringsAsFactors = F))
map_full <- do.call(expand.grid, dlist)
### filters
# non ATE methods: nearest, optimal, genetic
map_full <- map_full %>% filter(!(estimand == "ATE" & method %in% c("nearest", "optimal", "genetic")))

map <- map_full[map_idx,]

out <- mclapply(seeds, function(seed) {
  set.seed(seed)
  
  ##########################################################
  ### Generate simulation data
  ##########################################################
  
  simul <- SimulateRegressionTreatment(n=map$n, pk = map$pk, treat_p = map$treat_p, treat_beta=map$treat_beta, nu_conf = map$nu_conf, nu_xy = map$nu_xy, complexity = map$complexity)
  data <- data.frame(y=simul$ydata[,1], treat=simul$treat, simul$xdata)
  
  ############################
  # RUN #
  ############################
  
  res <- main(data=data, map=map)
  res <- with(map, data.frame(map_idx=map_idx, seed=seed,
                              res))
  return(res)
#})
}, mc.cores = 4)

out <- do.call(rbind.data.frame, out)

############################
# SAVE #
############################

map_savepath <- paste0("../Results/res_mod",mod,"_map.csv")
if(!file.exists(map_savepath)) {
  write.csv(map_full, map_savepath, row.names=F)
}

if (map$method == "null") {
  write.csv(out, paste0("../Results/res_mod",mod,"_j", map_idx, "_null.csv"), row.names = F)
} else {
  write.csv(out, paste0("../Results/res_mod",mod,"_j", map_idx, ".csv"), row.names = F)
}

#################
###### END ######
#################
combinations <- do.call(expand.grid, mod_specs)
