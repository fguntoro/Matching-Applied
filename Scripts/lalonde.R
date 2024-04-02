rm(list = ls())

##########################################################
### Matching on lalonde dataset
##########################################################
data <- MatchIt::lalonde

##########################################################
### MatchIt
##########################################################
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

dir.create("../Results/", showWarnings = F)

seed = 1
set.seed(seed)

map_idx = args[1]

covar_names <- colnames(data)[c(-1, -9)]

############################
# SPECIFICATIONS #
############################
f1 <- paste0("treat ~ ", paste0(c("age", "educ"), collapse = "+"))
f2 <- paste0("treat ~ ", paste0(covar_names, collapse = "+"))
match.formulas <- c(f1, f2)

outcome <- "re78"
predictor <- "treat"
f1 <- paste0(outcome, " ~ ", predictor)
f2 <- paste0(outcome, " ~ ", predictor, paste0(c("age", "educ"), collapse = "+"))
f3 <- paste0(outcome, " ~ ", predictor, paste0(covar_names, collapse = "+"))
stat.formulas <- c(f1, f2, f3)
family = "gaussian"

dat="data"
methods=c("nearest", "quick", "optimal", "genetic", "full")
distances=c("glm", "gam", "elasticnet", "randomforest", "cbps", "nnet", "bart", "robust_mahalanobis", "scaled_euclidean")
# ratio=c(1,3,5)

covar.names <- paste0("c(", paste0("\"", covar_names, "\"", collapse = ",") ,")")

############################
# RUN #
############################
map <- expand.grid(seed=seed, f=match.formulas, dat=dat, method=methods, distance=distances, stringsAsFactors = F, ratio = 1,
                   stat.formula = stat.formulas, family=family, covar.names = covar.names)
map <- map[map_idx,]

res <- main(seed=map$seed, formula.str=map$f, data_str=map$dat, method=map$method, distance=map$distance, ratio = map$ratio,
            stat.formula = map$stat.formula, family=map$family, covar.names=map$covar.names)

res <- cbind(data.frame(jobid = map_idx), res)

############################
# SAVE #
############################

saveRDS(res, paste0("../Results/res_mod2_unmatched_j", map_idx, ".rds"))
write.csv(res, paste0("../Results/res_mod2_unmatched_j", map_idx, ".csv"), row.names = F)