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

args = commandArgs(trailingOnly=TRUE)
map_idx = args[1]
#map_idx = 1

data_name= "cancer"
foldername = paste0(data_name)
dir.create(paste0("../Results/", foldername), showWarnings = F)
dir.create(paste0("~/../ephemeral/Matching/Results/", foldername), showWarnings = F)

############################
# SPECIFICATIONS #
############################
seeds = 1:2
exp=2

if (data_name == "lalonde") {
  data = MatchIt::lalonde
  data = model.matrix(~ 0 + ., data)#[,-1]
  data = as.data.frame(data)
  map_full <- get_exp_specs(exp, real_dat = data)
  map_full <- map_full[c(-1:-8)]
  map_full$f <- "treat ~ age + educ + raceblack + racehispan + racewhite + married + nodegree + re74 + re75"
  map_full$stat.formula <- "re78 ~ treat + age + educ + racehispan + racewhite + married + nodegree + re74 + re75"
} else if (data_name == "cancer") {
  # Load data
  data = readRDS("~/Cancer/Data/data_unmatched.rds")
  col_dict <- read.csv("~/Cancer/column_dictionary.csv")
  col_names <- subset(col_dict, Include == 1 & !(Type %in% c("Metabolic","Biochemistry","CancerRegistry", "DeathRegistry")))$CodingName
  col_names <- col_names[which(col_names != "date_recr")]
  
  # Filter by vaccination date
  vax_date <- as.Date(c("2020-12-01"))
  data <- subset(data, is.na(specdate) | specdate <= vax_date)
  
  # Filter by number of year censored pre-covid 2020
  n_censor_year = 5
  data$deceased_out <- if_else(data$deceased == 1 & data$death_by_any_cancer == 1, 1,
                               if_else(data$deceased == 1 & data$death_by_any_cancer == 0, NA, 0))
  data$deceased_censored <- if_else(data[,"date_diagnosis_latest_any_cancer_precovid"] >= (as.Date("2020-01-01", format="%Y-%m-%d") - 365.25 * n_censor_year), NA, data$deceased_out)
  data$deceased_censored <- factor(data$deceased_censored, levels = c(0,1))
  data <- subset(data, !is.na(deceased_censored))
  
  # Keep relevant columns
  covar_names <- col_names[col_names %in% colnames(data)]
  data <- data %>%
    dplyr::select(all_of(c("deceased_censored", "variable", covar_names)))

  data = model.matrix(~ 0 + ., data, contrasts.arg = lapply(data[, sapply(data, is.factor), drop = FALSE],
                                                          contrasts, contrasts = FALSE))
  data = as.data.frame(data)
  covar_names_sub <- colnames(data)[c(-1:-5)]
  
  map_full <- get_exp_specs(exp)
  map_full <- map_full[c(-1:-8)]
  map_full$f <- paste0("variable1 ~ ", paste0(covar_names_sub, collapse = "+"))
  # "treat ~ age + educ + raceblack + racehispan + racewhite + married + nodegree + re74 + re75"
  map_full$stat.formula <- paste0("deceased_censored1 ~ variable1 +", paste0(covar_names_sub, collapse = "+"))
  map_full$family <- "quasibinomial"
}

map <- map_full[map_idx,]
out <- lapply(seeds, function(seed) {
  set.seed(seed)
  
  ############################
  # RUN #
  ############################
  
  res <- main(data=data, map=map, full_covar = T)
  res <- with(map, data.frame(map_idx=map_idx, seed=seed,
                              res))
  return(res)
  })
#}, mc.cores = 24)

out <- do.call(rbind.data.frame, out)
out <- cbind(exp=exp, out)

if(map$method == "null") {
  start_var <- "time.sec"
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

map_savepath <- paste0("../Results/", foldername, "/res_exp",exp,"_map.csv")
write.csv(map_full, map_savepath, row.names=T)

if (map$method == "null") {
  write.csv(out, paste0("~/../ephemeral/Matching/Results/", foldername,"/res_j", map_idx, "_null_iter.csv"), row.names = F)
  write.csv(out.summary, paste0("../Results/", foldername, "/res_j", map_idx, "_null_summary.csv"), row.names = F)
} else if (!is.null(map$pk)) {
  if (map$pk == 1) {
    write.csv(out, paste0("~/../ephemeral/Matching/Results/", foldername, "/res_j", map_idx, "_pk1_iter.csv"), row.names = F)
    write.csv(out.summary, paste0("../Results/", foldername, "/res_j", map_idx, "_pk1_summary.csv"), row.names = F)
  } else {
    write.csv(out, paste0("~/../ephemeral/Matching/Results/", foldername, "/res_j", map_idx, "_norm_iter.csv"), row.names = F)
    write.csv(out.summary, paste0("../Results/", foldername, "/res_j", map_idx, "_norm_summary.csv"), row.names = F)
  }
} else {
  write.csv(out, paste0("~/../ephemeral/Matching/Results/", foldername, "/res_j", map_idx, "_norm_iter.csv"), row.names = F)
  write.csv(out.summary, paste0("../Results/", foldername, "/res_j", map_idx, "_norm_summary.csv"), row.names = F)
}

#################
###### END ######
#################