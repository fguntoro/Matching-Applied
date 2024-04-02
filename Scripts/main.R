rm(list=ls())

library(dplyr)
library(parallel)
library(MatchIt)
library(gdata)
library(survival)
library(marginaleffects)

source("~/Matching/Matching-Applied/Scripts/match_functions.R")
args = commandArgs(trailingOnly=TRUE)

dir.create("../Results/", showWarnings = F)

seed = 1
set.seed(seed)

data <- readRDS("~/CVD/data/ukb_mod2_unmatched.rds")
#data <- data[sample(nrow(data), 10000),]

map_idx = args[1]

f1 <- paste0("result ~ ", paste0(c("sexMale", "age"), collapse = "+"))
dat="data"
methods=c("nearest", "quick") # "optimal", "genetic", "full"
distances=c("glm", "gam", "elasticnet", "randomforest", "cbps", "nnet", "bart", "robust_mahalanobis", "scaled_euclidean")
# ratio=c(1,3,5)

map <- expand.grid(seed=seed, f=f1, dat=dat, method=methods, distance=distances, stringsAsFactors = F, ratio = 1,
                   outcome="case", predictor="result", covar_names='c("sexMale", "age")')
map <- map[map_idx,]

# test with matching map all base
# peakRAM({
#   timestart <- Sys.time()

#   cl <- makeCluster(3, outfile="")
#   clusterEvalQ(cl, {
#     library(MatchIt)
#     library(gdata)
#     library(survival)
#     library(marginaleffects)
#     library(dplyr)
#   })
#   clusterExport(cl, c("data", "MatchingFun", "StatsFun", "StatsGet"))
#   results <- clusterMap(cl, main, seed, formula_str=map$f, data_str=map$dat, method=map$method, distance=map$distance, ratio = map$ratio,
#                         outcome=map$outcome, predictor=map$predictor, covar_names=map$covar_names)
#   stopCluster(cl)
#   
#   res <- do.call(bind_rows, results)
#   
#   timeend <- Sys.time()
#   timetaken <- timeend - timestart
#   print(timetaken)
# }
# )
# 
res <- main(seed=map$seed, formula_str=map$f, data_str=map$dat, method=map$method, distance=map$distance, ratio = map$ratio,
     outcome=map$outcome, predictor=map$predictor, covar_names=map$covar_names)
 
# match_out <- MatchingFun(formula_str=map$f, data_str=map$dat, method=map$method, distance=map$distance, ratio = map$ratio)
# stats_out <- StatsFun(outcome=map$outcome, predictor=map$predictor, covar_names=map$covar_names, data=match_out$m.dat, family = "binomial", weights = match_out$m.dat$weights)
res <- cbind(data.frame(jobid = map_idx), res)

saveRDS(res, paste0("../Results/res_mod2_unmatched_j", map_idx, ".rds"))
write.csv(res, paste0("../Results/res_mod2_unmatched_j", map_idx, ".csv"), row.names = F)