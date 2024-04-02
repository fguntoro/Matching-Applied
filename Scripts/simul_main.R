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
library("lmtest") #coeftest
library("sandwich") #vcovCL

source("~/Matching/Matching-Applied/Scripts/match_functions.R")
args = commandArgs(trailingOnly=TRUE)


##########################################################
### Matching on simul dataset
##########################################################
data <- tibble(
    # x is our confounder
    x = runif(1000, -1, 4),
    # it affects the probabiity of receiving the treatment
    # in a NON-LINEAR way (step function)
    prob_d = ifelse(x > 0.5 & x < 2.5, 0.1, 0.9),
    d = rbinom(1000, 1, prob_d),
    noise = rnorm(1000, sd = 0.1),
    # for simplification, the treatment effect is homogeneous
    treat_effect = 1,
    # x also effects the outcome in a non-linear way
    y = sin(x) + d*treat_effect + noise
    # y = 2 * x + d*treat_effect + noise
  ) %>% 
  mutate(d_factor = factor(d,
                           levels=c(0,1), labels=c("Untreated",
                                                   "Treated")))


##########################################################
### MatchIt
##########################################################

dir.create("../Results/", showWarnings = F)

seed = 1
set.seed(seed)

map_idx = args[1]

covar_names = c("d", "x")

############################
# SPECIFICATIONS #
############################
outcome <- "y"
predictor <- "d"
confounder <- "x"

f1 <- paste0(predictor, " ~ ", paste0(confounder, collapse = "+"))
match.formulas <- c(f1)


f1 <- paste0(outcome, " ~ ", predictor)
f2 <- paste0(outcome, " ~ ", predictor, "+", paste0(confounder, collapse = "+"))
stat.formulas <- c(f1, f2)
family = "gaussian"

dat="data"
methods=c("nearest", "quick", "optimal", "genetic", "full")
distances=c("glm", "gam", "randomforest", "cbps", "nnet", "bart", "robust_mahalanobis", "scaled_euclidean")
# distances=c("glm", "gam", "elasticnet", "randomforest", "cbps", "nnet", "bart", "robust_mahalanobis", "scaled_euclidean")
# ratio=c(1,3,5)

covar.names <- paste0("c(", paste0("\"", covar_names, "\"", collapse = ",") ,")")

############################
# RUN #
############################
map_full <- expand.grid(seed=seed, f=match.formulas, dat=dat, method=methods, distance=distances, stringsAsFactors = F, ratio = 1,
                   stat.formula = stat.formulas, family=family, covar.names = covar.names)
#map <- map[map_idx,]

res_df <- data.frame()
for (map_idx in 1:nrow(map_full)) {
  map <- map_full[map_idx,]
  res <- main(seed=map$seed, formula.str=map$f, data_str=map$dat, method=map$method, distance=map$distance, ratio = map$ratio,
              stat.formula = map$stat.formula, family=map$family, covar.names=map$covar.names)
  res_df <- rbind(res_df, res)
}

res <- main(seed=map$seed, formula.str=map$f, data_str=map$dat, method=map$method, distance=map$distance, ratio = map$ratio,
            stat.formula = map$stat.formula, family=map$family, covar.names=map$covar.names)

stats_out <- StatsFun(stat.formula=map$stat.formula, covar.names=map$covar.names, data=data, family = map$family, weights = rep(1,nrow(data)))

res <- cbind(data.frame(jobid = map_idx), res)

############################
# SAVE #
############################

saveRDS(res, paste0("../Results/res_mod2_unmatched_j", map_idx, ".rds"))
write.csv(res, paste0("../Results/res_mod2_unmatched_j", map_idx, ".csv"), row.names = F)