library(dplyr)

get_exp_specs <- function(exp = 1) {
  ### get default
  dlist_default = get_default_spec()
  
  sub_exps <- c("n", "pk", "nu_xy", "nu_conf", "treat_p", "treat_beta", "ev_xy", "complexity")
  
  map_full <- data.frame()
  
  if (exp == 1) {
    ### Experiment 1: Simulation specs
    for(sub_exp in sub_exps) {
      dlist <- dlist_default
      if (sub_exp == "n") {
        dlist$n = c(1000, 10000, 100000)
      }
      if (sub_exp == "pk") {
        dlist$pk = c(1, 10, 100)
      }
      if (sub_exp == "nu_xy") {
        dlist$pk = c(100)
        dlist$nu_xy = c(0.01, 0.02, 0.05, 0.1, 1)
      }
      if (sub_exp == "nu_conf") {
        dlist$pk = c(100)
        dlist$nu_conf = c(0.01, 0.02, 0.05, 0.1, 1)
      }
      if (sub_exp == "treat_p") {
        dlist$treat_p = c(0.2, 0.5, 0.8)
      }
      if (sub_exp == "treat_beta") {
        dlist$treat_beta = c(0.1, 0.5, 1)
      }
      if (sub_exp == "ev_xy") {
        dlist$ev_xy = c(0.5, 0.7, 0.99)
      }
      if (sub_exp == "complexity") {
        dlist$complexity = c("linear", "sq", "sin")
      }
      map_tmp <- do.call(expand.grid, dlist)
      map_full <- rbind(map_full, map_tmp)
      
    }
    

  }
  
  ### filters
  
  map_full <- map_full %>%
    # non ATE methods: nearest, optimal, genetic
    filter(!(estimand == "ATE" & method %in% c("nearest", "optimal", "genetic"))) %>%
    # optimal, full above 10000 n too much memory, genetic too long
    filter(!(n>10000 & method %in% c("optimal", "full", "genetic")))
  map_full <- map_full %>% distinct()
  
  return(map_full)
}

get_default_spec <- function() {
  # simulation specs
  n = 1000
  pk = 1
  nu_xy = 1
  nu_conf = 1
  treat_p = 0.5
  treat_beta = 1
  ev_xy=0.99
  complexity = "linear"
  
  ### matching specs
  dat="data"
  treat <- "treat"
  #confounder <- "var1"
  #f1 <- paste0(treat, " ~ ", paste0(confounder, collapse = "+"))
  # confounder <- paste0("var", 1:pk)
  # f2 <- paste0(treat, " ~ ", paste0(confounder, collapse = "+"))
  f3 <- "simul_selected"
  match.formulas <- c(f3)
  
  estimand = c("ATT")
  methods=c("null","nearest", "quick", "optimal", "genetic", "full")
  distances=c("glm", "gam", "randomforest", "cbps", "nnet", "bart", "robust_mahalanobis", "scaled_euclidean")
  ratio=1
  caliper="null"
  replace="F"
  discard="none"
  
  ### estimation specs
  outcome <- "y"
  # f1 <- paste0(outcome, " ~ ", treat, "+", paste0(confounder, collapse = "+"))
  f2 <- "simul_selected"
  stat.formulas <- c(f2)
  family = "gaussian"
  
  ### return list
  dlist <- list(n=n, pk = pk,nu_xy = nu_xy,nu_conf = nu_conf,treat_p = treat_p,treat_beta = treat_beta,ev_xy=ev_xy, complexity = complexity,
                dat=dat, f=match.formulas, estimand=estimand, method=methods, distance=distances, ratio = ratio,caliper=caliper, replace=replace,
                    stat.formula = stat.formulas, family=family, stringsAsFactors = F)
  return(dlist)
}

get_mod_specs <- function(mod = 1) {
  ### simulation specs
  
  n = c(1000)
  pk = 1
  nu_xy = 1
  nu_conf = c(1)
  treat_p = c(0.2, 0.5, 0.8)
  treat_beta = c(0.1, 0.5, 1)
  complexity = "linear"
  
  if(mod == 1) {
    type="onevar_linear"
    
  } else if (mod == 2) {
    type="multivar_linear"
    pk = c(100)
    nu_xy = c(0.1, 1)
    nu_conf = c(0.1, 1)
    
  } else if (mod == 3) {
    type="onevar_nonlinear"
    complexity="sin"
    
  } else if (mod == 4) {
    type= "multivar_nonlinear"
    pk = c(100)
    nu_xy = c(0.1, 1)
    nu_conf = c(0.1, 1)
    treat_p = c(0.2, 0.5, 0.8)
    complexity = "sin"
  }
  
  mod <- list(mod=mod, type=type, n=n, pk=pk, nu_xy=nu_xy, nu_conf=nu_conf, treat_p=treat_p, treat_beta=treat_beta, complexity=complexity)
  
  return(mod)
}

