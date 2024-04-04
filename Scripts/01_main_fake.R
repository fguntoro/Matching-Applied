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
library(stringr)

source("~/Matching/Matching-Applied/Scripts/match_functions.R")
source("~/Matching/Matching-Applied/Scripts/simul.R")
source("~/Matching/Matching-Applied/Scripts/mod_definitions.R")

setwd("/rds/general/user/fg520/home/Matching/Matching-Applied/Scripts/")

dir.create("../Results/", showWarnings = F)
dir.create("../Results/unmatched/", showWarnings = F)
dir.create("../Results/matched/", showWarnings = F)

args = commandArgs(trailingOnly=TRUE)
map_idx = args[1]
#map_idx = 12

############################
# SPECIFICATIONS #
############################
# One variable - linear
seeds = 1:100
exp=2
map_full <- get_exp_specs(exp)
map <- map_full[map_idx,]

out <- mclapply(seeds, function(seed) {
  set.seed(seed)
  
  ##########################################################
  ### Generate simulation data
  ##########################################################
  
  simul <- SimulateRegressionTreatment(n=map$n, pk = map$pk, treat_p = map$treat_p, treat_beta=map$treat_beta, nu_conf = map$nu_conf, nu_xy = map$nu_xy, complexity = map$complexity, ev_xy=map$ev_xy)
  
  ############################
  # RUN #
  ############################
  
  res <- main(data=simul, map=map)
  res <- with(map, data.frame(map_idx=map_idx, seed=seed,
                              res))
  return(res)
#})
}, mc.cores = 24)

out <- do.call(rbind.data.frame, out)
out <- cbind(exp=exp, out)

if(map$method == "null") {
  start_var <- "time"
  end_var <- "adjusted_pval"
} else {
  start_var <- "ESS.total.perc"
  end_var <- "marginal_pval"
}

out.summary.mean <- out %>%
  summarise(across(start_var:end_var, \(x) mean(x, na.rm=T))) %>%
  rename_with(~str_c(., "_mean"), everything())

out.summary.sd <- out %>%
  summarise(across(start_var:end_var, \(x) sd(x, na.rm=T))) %>%
  rename_with(~str_c(., "_sd"), everything())

out.summary <- cbind(exp=exp, map_idx=map_idx, n_seeds=length(seeds), map, out.summary.mean, out.summary.sd)

############################
# SAVE #
############################

map_savepath <- paste0("../Results/res_exp",exp,"_map.csv")
write.csv(map_full, map_savepath, row.names=T)

if (map$method == "null") {
  write.csv(out, paste0("~/../ephemeral/Matching/Results/res_exp",exp,"_j", map_idx, "_null_iter.csv"), row.names = F)
  write.csv(out.summary, paste0("../Results/res_exp",exp,"_j", map_idx, "_null_summary.csv"), row.names = F)
} else if (pk == 1) {
  write.csv(out, paste0("~/../ephemeral/Matching/Results/res_exp",exp,"_j", map_idx, "_pk1_iter.csv"), row.names = F)
  write.csv(out.summary, paste0("../Results/res_exp",exp,"_j", map_idx, "_pk1_summary.csv"), row.names = F)
} else {
  write.csv(out, paste0("~/../ephemeral/Matching/Results/res_exp",exp,"_j", map_idx, "_iter.csv"), row.names = F)
  write.csv(out.summary, paste0("../Results/res_exp",exp,"_j", map_idx, "_summary.csv"), row.names = F)
}

#################
###### END ######
#################