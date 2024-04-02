main <- function(data, map) {
  formula.str <- map$f
  estimand <- map$estimand
  method <- map$method
  distance <- map$distance
  ratio <- map$ratio
  stat.formula <- map$stat.formula
  family <- map$family
  covar.names <- eval(parse(text = map$covar.names))
  treat_beta = map$treat_beta
  
  match_out <- MatchingFun(formula.str=formula.str, data=data, estimand=estimand, method=method, distance=distance, ratio=ratio, covar.names=covar.names)
  stats_out <- StatsFun(stat.formula=stat.formula, covar.names=covar.names, data=match_out$m.dat, family = family, weights = match_out$m.dat$weights, method = method, treat_beta=treat_beta)
  
  res <- data.frame(as.list(c(map, match_out$res, stats_out)))
  
  return(res)
}

rescale <- function(x){(x-min(x))/(max(x)-min(x))}

MatchingFun <- function(formula.str, data, estimand, method, distance, ratio, covar.names){
  print(paste0("Running ", method, " with distance ", distance, " and ratio ", ratio))
  
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
  
  suppressWarnings(
  if (method == "genetic") {
    m.out <- matchit(as.formula(formula.str), data = data, estimand = estimand,
                     method = method, distance = distance, distance.options = distance.options, pop.size = 100, ratio=ratio)
  } else if (method == "optimal") {
    m.out <- matchit(as.formula(formula.str), data = data, estimand = estimand,
                     method = method, distance = distance, distance.options = distance.options, tol = 1e-5, ratio=ratio)
  } else if (method == "null") {
    m.out <- matchit(as.formula(formula.str), data = data, estimand = estimand,
                     method = NULL)
  } else {
    m.out <- matchit(as.formula(formula.str), data = data, estimand = estimand,
                     method = method, distance = distance, distance.options = distance.options, ratio = ratio)
  }
  )
  
  #rm(data)
  
  timeEnd <- Sys.time()
  time <- difftime(timeEnd, timeStart, units="secs")
  print(time)
  
  m.dat <- match.data(m.out)
  
  summ <- summary(m.out, addlvariables=covar.names, standardize=T, improvement=T, pair.dist=F)

  # Get sample sizes including ESS, ESS/total%
  row.names(summ$nn)[c(1,3)] <- c("All.ESS", "Matched.ESS")
  ESS.total.perc <- round(sum(summ$nn[3,]) / sum(summ$nn[1,]) * 100, 2)
  res <- c(unmatrix(summ$nn, byrow = T), ESS.total.perc = ESS.total.perc)
           
  # if(method == "null") {
  #   res <- c(unmatrix(summ$nn, byrow = T),
  #            rep(NA, 7*length(covar.names)),
  #            time=time)
  # } else {
  #   res <- c(unmatrix(summ$nn, byrow = T),
  #            unmatrix(summ$sum.matched, byrow = T),
  #            time=time)
  # }
  
  # Get SMD, var ratio, and KS
  summ.covars <- cbind(summ$sum.matched[covar.names,c(3,4,6)], summ$reduction[covar.names, c(1,2,4)])
  summ.covars[,3] <- abs(summ.covars[,1])
  colnames(summ.covars) <- c("SMD", "Var.ratio", "KS", "SMD.PBR", "Var.ratio.log.PBR", "KS.PBR")
  
  if(length(covar.names) > 1) {
    summ.covars <- as.matrix(data.frame(mean = colMeans(summ.covars, na.rm = T),
                                        median = apply(summ.covars,2,median,na.rm = T),
                                        max = apply(summ.covars,2,max)))
    summ.covars <- unmatrix(summ.covars, byrow = T)
  }
  
  res <- c(res, summ.covars, time.sec=time)
  res <- data.frame(as.list(res))
  res <- res %>% dplyr::select(-contains("distance.", ignore.case=F))
  
  out <- list(m.dat=m.dat, res=res)
  
  return(out)
}


StatsFun <- function(types=c("unadjusted", "adjusted", "conditional", "mixed", "marginal"), stat.formula, covar.names, data, family="quasibinomial", weights, method, treat_beta) {
  
  if (family == "gaussian") {
    types=c("unadjusted", "adjusted", "mixed", "marginal")
    comparison = "differenceavg"
  } else {
    comparison = "lnoravg"
  }
  
  if (method == "null") {
    types=c("unadjusted", "adjusted")
  }
  
  stat.formula.split <- strsplit(stat.formula, "~")[[1]]
  outcome <- gsub(" ", "", stat.formula.split[1])
  
  stat.formula.predictors <- strsplit(stat.formula.split[2], "\\+")[[1]]
  predictor <- gsub(" ", "", stat.formula.predictors[1])
  
  dat_outcome <- data[,outcome]
  dat_predictor <- factor(data[,predictor], levels = c(0,1))
  
  #tab <- table(dat_outcome, dat_predictor, dnn=c(outcome, predictor))
  #print(tab)
  #res <- c(sum(tab), tab[1,1],tab[2,1], tab[1,2], tab[2,2])
  #names(res) <- c("total_n", do.call(paste0, expand.grid(outcome, c(0,1), "_", predictor, c(0,1))))
  res <- c()
  
  # if (family != "gaussian" & 0 %in% tab) {
  #   stats_empty <- rep(NA, length(types)*4)
  #   names(stats_empty) <- do.call(paste0, expand.grid(types,"_", c("estimate", "confint.lower", "confint.upper", "pval")) %>% arrange_all)
  #   res <- c(res, stats_empty)
  # 
  # } else {
    t.test.res <- t.test(data$y[which(data$treat == 1)], data$y[which(data$treat == 0)])
    t.test.res <- c(t.test.res$estimate, t.test.res$p.value)
    names(t.test.res) <- c("treat1.mean", "treat2.mean", "t.test.pval")
    
    for (type in types) {
      
      if (type == "unadjusted") {
        # unadjusted
        fm <- formula(paste0(outcome, "~", predictor))
        fit <- glm(fm, data = data, family=family, weights = weights)
        
      } else if (type == "adjusted") {
        # adjusted
        fm <- formula(paste0(outcome, "~", predictor, "+", paste0(covar.names, collapse = "+")))
        fit <- glm(fm, data = data, family=family, weights = weights)
        
      } else if (type == "conditional") {
        # conditional logistic regression
        fm <- formula(paste0(outcome, "~", predictor, "+", paste0(covar.names, collapse = "+"), " + strata(subclass)" ))
        fit <- clogit(fm, method="approximate", data = data, weights = weights)
        
      } else if (type == "mixed") {
        # conditional logistic regression
        fm <- formula(paste0(outcome, "~", predictor, "+ (1 | subclass) +", paste0(covar.names, collapse = "+")))
        if (family == "gaussian") {
          fit <- suppressWarnings(lmer(fm, data = data, weights = weights))
        } else {
          fit <- suppressWarnings(glmer(fm, data = data, weights = weights, family=family))
        }
        
      } else if (type == "marginal") {
        fm <- formula(paste0(outcome, "~", predictor, "+", paste0(covar.names, collapse = "+")))
        fit <- glm(fm, data = data, family=family, weights = weights)
        # marginal effect
        fit <- avg_comparisons(fit, variables = predictor,
                               vcov = ~subclass,
                               newdata = subset(data, get(predictor) == 1), comparison = comparison,
                               wts = "weights")
      }
      res <- c(res, StatsGet(fit, type, family, treat_beta))
    }
  # }
  
  res <- data.frame(as.list(t.test.res), as.list(res))
}


StatsGet <- function(fit, type, family, treat_beta) {
  
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
    fit.confint <- confint(fit, method="Wald")
    lower.bound <- fit.confint[4,1]
    upper.bound <- fit.confint[4,2]
    mse <- mean(residuals(fit, type="pearson")^2)
    
    # to get pval
    fm <- update.formula(formula(fit), . ~ . -treat)
    data <- fit@frame
    colnames(data)[ncol(data)] <- "weights"
    if (family == "gaussian") {
      fit0 <- suppressWarnings(lmer(fm, data = data, weights = weights))
    } else {
      fit0 <- suppressWarnings(glmer(fm, data = data, weights = weights, family=family))
    }
    pval <- anova(fit0, fit)[2,8]
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
    fit.confint <- suppressMessages(confint(fit, method="Wald"))
    lower.bound <- fit.confint[idx,1]
    upper.bound <- fit.confint[idx,2]
    pval = summary(fit)$coefficients[idx, 4]
    mse <- mean(residuals(fit, type="pearson")^2)
  }
  
  # sample.se <- summary(fit)$coef[2,2]
  # degrees.freedom = fit$df.null
  # alpha = 0.05
  # t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
  # margin.error <- t.score * sample.se
  # lower.bound <- estimate - margin.error
  # upper.bound <- estimate + margin.error
  
  res <- data.frame(
    estimate = estimate,
    std.error = std.error,
    confint.lower = lower.bound,
    confint.upper = upper.bound,
    squared.error = (estimate - treat_beta)^2,
    confint.coverage = as.integer(lower.bound < treat_beta & treat_beta < upper.bound),
    mse = mse
  )
  
  if (family != "gaussian") {
    res <- exp(res)
  }

  res <- cbind(res, pval=pval)
  
  colnames(res) <- paste0(type, "_", colnames(res))
  
  return(res)
}
