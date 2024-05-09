main <- function(data, map, full_covar = F) {
  estimand <- map$estimand
  method <- map$method
  distance <- map$distance
  ratio <- map$ratio
  family <- map$family
  treat_beta = map$treat_beta
  caliper <- map$caliper
  replace <- map$replace
  discard <- map$discard
  
  if(class(data) == "simulation_regression") {
    if(map$f == "simul_selected"){
      true_conf <- data$theta_conf * data$theta
      confounder <- paste0("var", which(true_conf==1))
      formula.str <- paste0("treat", " ~ ", paste0(confounder, collapse = "+"))
    } else {
      formula.str <- map$f
    }
    
    if(map$stat.formula == "simul_selected") {
      confounder <- paste0("var", which(data$theta_conf==1))
      stat.formula <- paste0("y", " ~ treat + ", paste0(confounder, collapse = "+"))
    } else {
      stat.formula <- map$stat.formula
    }
    
    data <- with(data, data.frame(y=ydata[,1], treat=treat, xdata))
    
  } else {
    data <- data
    formula.str <- map$f
    stat.formula <- map$stat.formula
    
    #TODO
  }
  
  match.fm <- as.formula(formula.str)
  stat.fm <- as.formula(stat.formula)
  
  match_out <- MatchingFun(formula=match.fm, data=data, estimand=estimand, method=method, distance=distance, ratio=ratio, caliper, replace, discard, full_covar)
  stats_out <- StatsFun(formula=stat.fm, data=match_out$m.dat, family = family, weights = match_out$m.dat$weights, method = method, treat_beta=treat_beta)
  
  res <- data.frame(as.list(c(map, match_out$res, stats_out)))
  
  return(res)
}

rescale <- function(x){(x-min(x))/(max(x)-min(x))}

MatchingFun <- function(formula, data, estimand, method, distance, ratio, caliper, replace, discard, full_covar){
  print(paste0("Running ", method, " with distance ", distance, " and ratio ", ratio))
  
  covar.names <- labels(terms(formula))
  
  timeStart <- Sys.time()
  if (distance == "gbm") {
    distance.options = list(n.cores=8)
  } else if (distance == "nnet") {
    distance.options = list(size=50)
  } else {
    distance.options = list()
  }
  
  if (ratio == 0) {
    ratio = NULL
  } else {
    ratio = ratio
  }
  
  if (caliper == "null") {
    caliper = NULL
  } else {
    caliper = as.numeric(caliper)
  }
  
  if (replace == "F") {
    replace=F
  } else if (replace == "T") {
    replace=T
  }
  
  suppressWarnings(
  if (method == "genetic") {
    m.out <- matchit(formula, data = data, estimand = estimand, caliper=caliper, replace=replace, discard=discard,
                     method = method, distance = distance, distance.options = distance.options, pop.size = 100, ratio=ratio)
  } else if (method == "optimal") {
    m.out <- matchit(formula, data = data, estimand = estimand, caliper=caliper, replace=replace, discard=discard,
                     method = method, distance = distance, distance.options = distance.options, tol = 1e-5, ratio=ratio)
  } else if (method == "null") {
    m.out <- matchit(formula, data = data, estimand = estimand, caliper=caliper, replace=replace, discard=discard,
                     method = NULL)
  } else {
    m.out <- matchit(formula, data = data, estimand = estimand, caliper=caliper, replace=replace, discard=discard,
                     method = method, distance = distance, distance.options = distance.options, ratio = ratio)
  }
  )
  
  #rm(data)
  
  timeEnd <- Sys.time()
  time <- difftime(timeEnd, timeStart, units="secs")
  print(time)
  
  m.dat <- if(replace == T & method == "nearest" & distance %in% c("scaled_euclidean", "robust_mahalanobis")) get_matches(m.out) else match.data(m.out)
  
  summ <- summary(m.out, addlvariables=covar.names, standardize=T, improvement=T, pair.dist=F)

  if(method == "null") {
    res <- c(unmatrix(summ$nn, byrow = T))
    
    if(length(covar.names) > 1) {
      # Get SMD, var ratio, and KS
      summ.covars <- cbind(summ$sum.all[covar.names,c(3,4,6)])
      summ.covars[,1] <- abs(summ.covars[,1])
      colnames(summ.covars) <- c("SMD", "Var.ratio", "KS")
      
      if(!full_covar) {
        summ.covars <- as.matrix(data.frame(mean = apply(summ.covars,2,mean, na.rm = T),
                                            median = apply(summ.covars,2,median,na.rm = T),
                                            max = apply(summ.covars,2,max)))
      }
      
      summ.covars <- unmatrix(summ.covars, byrow = T)
      
    } else {
      # Get SMD, var ratio, and KS
      summ.covars <- c(summ$sum.all[covar.names,c(3,4,6)])
      summ.covars[1] <- abs(summ.covars[1])
      names(summ.covars) <- c("SMD", "Var.ratio", "KS")
    }
    
  } else {
    # Get sample sizes including ESS, ESS/total%
    row.names(summ$nn)[c(1,3)] <- c("All.ESS", "Matched.ESS")
    ESS.total.perc <- round(sum(summ$nn[3,]) / sum(summ$nn[1,]) * 100, 2)
    res <- c(unmatrix(summ$nn, byrow = T), ESS.total.perc = ESS.total.perc)
    
    if(length(covar.names) > 1) {
      # Get SMD, var ratio, and KS
      summ.covars <- cbind(summ$sum.matched[covar.names,c(3,4,6)], summ$reduction[covar.names, c(1,2,4)])
      summ.covars[,1] <- abs(summ.covars[,1])
      colnames(summ.covars) <- c("SMD", "Var.ratio", "KS", "SMD.PBR", "Var.ratio.log.PBR", "KS.PBR")
      
      if(!full_covar) {
        summ.covars <- as.matrix(data.frame(mean = apply(summ.covars,2,mean, na.rm = T),
                                            median = apply(summ.covars,2,median,na.rm = T),
                                            max = apply(summ.covars,2,max)))
      }

      summ.covars <- unmatrix(summ.covars, byrow = T)
      
    } else {
      # Get SMD, var ratio, and KS
      summ.covars <- c(summ$sum.matched[covar.names,c(3,4,6)], summ$reduction[covar.names, c(1,2,4)])
      summ.covars[1] <- abs(summ.covars[1])
      names(summ.covars) <- c("SMD", "Var.ratio", "KS", "SMD.PBR", "Var.ratio.log.PBR", "KS.PBR")
    }
    
  }
  res <- c(res, time.sec=time, summ.covars)
  res <- data.frame(as.list(res))
  res <- res %>% dplyr::select(-contains("distance.", ignore.case=F))
  
  out <- list(m.dat=m.dat, res=res)
  
  return(out)
}


StatsFun <- function(types=c("unadjusted", "adjusted", "conditional", "mixed", "marginal"), formula, data, family="quasibinomial", weights, method, treat_beta) {
  
  if (family == "gaussian") {
    types=c("unadjusted", "adjusted", "marginal")
    comparison = "differenceavg"
  } else {
    types=c("unadjusted", "adjusted", "marginal")
    comparison = "lnoravg"
  }
  
  if (method == "null") {
    types=c("unadjusted", "adjusted")
  }
  
  outcome <- as.character(formula[[2]])
  covar.names <- labels(terms(formula))
  expo <- covar.names[1]
  covar.names <- covar.names[-1]
  
  # dat_outcome <- data[,outcome]
  # dat_treat <- factor(data[,treat], levels = c(0,1))
  #tab <- table(dat_outcome, dat_treat, dnn=c(outcome, treat))
  #print(tab)
  #res <- c(sum(tab), tab[1,1],tab[2,1], tab[1,2], tab[2,2])
  #names(res) <- c("total_n", do.call(paste0, expand.grid(outcome, c(0,1), "_", treat, c(0,1))))
  res <- c()
  
  # if (family != "gaussian" & 0 %in% tab) {
  #   stats_empty <- rep(NA, length(types)*4)
  #   names(stats_empty) <- do.call(paste0, expand.grid(types,"_", c("estimate", "confint.lower", "confint.upper", "pval")) %>% arrange_all)
  #   res <- c(res, stats_empty)
  # 
  # } else {
    t.test.res <- t.test(data[which(data[expo] == 1), outcome], data[which(data[expo] == 0), outcome])
    # t.test.res <- t.test(data$y[which(data$treat == 1)], data$y[which(data$treat == 0)])
    t.test.res <- c(t.test.res$estimate, t.test.res$p.value)
    names(t.test.res) <- c("treat1.mean", "treat2.mean", "t.test.pval")
    
    for (type in types) {
      
      if (type == "unadjusted") {
        # unadjusted
        fm <- formula(paste0(outcome, "~", expo))
        
        #if (family == "binomial") {
        #  fit <- glm(fm, data = data, family="quasibinomial", weights = weights)
        #} else {
          fit <- glm(fm, data = data, family=family, weights = weights)
        #}
        
      } else if (type == "adjusted") {
        # adjusted
        fm <- formula(paste0(outcome, "~", expo, "+", paste0(covar.names, collapse = "+")))
        fit <- glm(fm, data = data, family=family, weights = weights)
        
      } else if (type == "conditional") {
        # conditional logistic regression
        fm <- formula(paste0(outcome, "~", expo, "+", paste0(covar.names, collapse = "+"), " + strata(subclass)" ))
        fit <- clogit(fm, method="approximate", data = data, weights = weights)
        
      } else if (type == "mixed") {
        # mixed logistic regression
        fm <- formula(paste0(outcome, "~", expo, "+ (1 | subclass) +", paste0(covar.names, collapse = "+")))
        if (family == "gaussian") {
          fit <- suppressWarnings(lmer(fm, data = data, weights = weights))
        } else {
          fit <- suppressWarnings(glmer(fm, data = data, weights = weights, family=family))
        }
        
      } else if (type == "marginal") {
        fm <- formula(paste0(outcome, "~", expo, "+", paste0(covar.names, collapse = "+")))
        fit <- glm(fm, data = data, family=family, weights = weights)
        # marginal effect
        fit <- avg_comparisons(fit, variables = expo,
                               vcov = ~subclass,
                               newdata = subset(data, get(expo) == 1), comparison = comparison,
                               wts = "weights")
      }
      res <- c(res, StatsGet(fit, type, family, treat_beta))
    }
  # }
  
  res <- data.frame(as.list(t.test.res), as.list(res))
}


StatsGet <- function(fit, type, family, treat_beta = NULL) {
  
  if (type %in% c("marginal")) {
    summ <- summary(fit)
    estimate <- summ$estimate
    std.error <- summ$std.error
    lower.bound <- summ$conf.low
    upper.bound <- summ$conf.high
    pval = summ$p.value
    mse <- NA

  } else if (type %in% c("mixed")) {
    summ <- summary(fit)
    estimate <- summ$coefficients[2,1]
    std.error <- summ$coefficients[2,2]
    
    degrees.freedom = fit$df.null #TODO fix this error return df, or just get confint below
    alpha = 0.05
    t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
    margin.error <- t.score * std.error
    lower.bound <- estimate - margin.error
    upper.bound <- estimate + margin.error
    
    #fit.confint <- confint(fit, method="Wald")
    #lower.bound <- fit.confint[4,1]
    #upper.bound <- fit.confint[4,2]
    mse <- mean(residuals(fit, type="deviance")^2)
    
    # to get pval
    fm <- update.formula(formula(fit), . ~ . -treat)
    
    fit0 <- update(fit, formula=drop.terms(terms(fit), 1, keep.response=TRUE))
    fit0 <- update(fit,. ~. -get(expo))
    tmp <- data %>%
      dplyr::select(-subclass) %>%
      apply(2, rescale) %>%
      as.data.frame()
    
    data <- fit@frame
    colnames(data)[ncol(data)] <- "weights"
    if (family == "gaussian") {
      fit0 <- suppressWarnings(lmer(fm, data = data, weights = weights))
    } else {
      fit0 <- suppressWarnings(glmer(fm, data = data, weights = weights, family=family))
    }
    pval <- anova(fit0, fit, test="Chisq")[2,5]
    # pval = suppressWarnings(drop1(fit, test="Chisq"))[2,4]
    
  } else {
    
    if (type %in% c("unadjusted", "adjusted")) {
      idx = 2
    } else if (type %in% c("conditional")) {
      idx = 1
    }
    summ <- summary(fit)
    estimate <- summary(fit)$coefficients[idx, 1]
    std.error <- summary(fit)$coefficients[idx, 2]
    degrees.freedom = fit$df.null
    alpha = 0.05
    t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
    margin.error <- t.score * std.error
    lower.bound <- estimate - margin.error
    upper.bound <- estimate + margin.error
    
    #fit.confint <- suppressMessages(confint.default(fit))
    #lower.bound <- fit.confint[idx,1]
    #upper.bound <- fit.confint[idx,2]
    pval = summary(fit)$coefficients[idx, 4]
    mse <- mean(residuals(fit, type="deviance")^2)
  }
  
  
  if(is.null(treat_beta)) {
    res <- data.frame(
      estimate = estimate,
      std.error = std.error,
      confint.lower = lower.bound,
      confint.upper = upper.bound,
      mse = mse
    )
  } else {
    res <- data.frame(
      estimate = estimate,
      std.error = std.error,
      confint.lower = lower.bound,
      confint.upper = upper.bound,
      squared.error = (estimate - treat_beta)^2,
      confint.coverage = as.integer(lower.bound < treat_beta & treat_beta < upper.bound),
      mse = mse
    )
  }
  

  
  if (family != "gaussian") {
    res <- exp(res)
  }

  res <- cbind(res, pval=pval)
  
  colnames(res) <- paste0(type, "_", colnames(res))
  
  return(res)
}
