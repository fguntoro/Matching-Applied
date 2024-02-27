main <- function(seed, formula_str, data_str, method, distance, ratio, outcome, predictor, covar_names) {
  set.seed(seed)
  
  match_out <- MatchingFun(formula_str, data_str, method, distance, ratio)
  stats_out <- StatsFun(outcome=outcome, predictor=predictor, covar_names=covar_names, data=match_out$m.dat, family = "binomial", weights = match_out$m.dat$weights)
  
  res <- data.frame(as.list(c(seed=seed,match_out$res, stats_out)))
  
  return(res)
}

MatchingFun <- function(formula_str, data_str, method, distance, ratio){
  print(paste0("Running ", method, " with distance ", distance, " and ratio ", ratio))
  
  data <- get(data_str)
  
  timeStart <- Sys.time()
  if (distance == "gbm") {
    distance.options = list(n.cores=1)
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
  
  if (method == "genetic") {
    m.out <- matchit(as.formula(formula_str), data = data,
                     method = method, distance = distance, distance.options = distance.options, pop.size = 100, ratio=ratio)
  } else if (method == "optimal") {
    m.out <- matchit(as.formula(formula_str), data = data,
                     method = method, distance = distance, distance.options = distance.options, tol = 1e-5, ratio=ratio)
  } else {
    m.out <- matchit(as.formula(formula_str), data = data,
                     method = method, distance = distance, distance.options = distance.options, ratio = ratio)
  }
  
  timeEnd <- Sys.time()
  time <- difftime(timeEnd, timeStart, units="secs")
  print(time)
  
  m.dat <- match.data(m.out)
  
  call.dat <- c(f=formula_str,
                data=data_str,
                method=method,
                distance=distance)
  
  summ <- summary(m.out)
  
  res <- c(unmatrix(summ$nn),
           unmatrix(summ$sum.matched),
           time=time)
  res <- data.frame(as.list(call.dat), as.list(res))
  
  out <- list(m.dat=m.dat, res=res)
  
  return(out)
}


StatsFun <- function(types=c("unadjusted", "adjusted", "conditional", "marginal"), outcome, predictor, covar_names, data, family="binomial", weights) {
  
  call.dat <- c(outcome=outcome,
                predictor=predictor,
                covar_names=covar_names)
  
  covar_names <- eval(parse(text = covar_names))
  
  tab <- table(data[,outcome], data[,predictor])
  print(tab)
  res <- c(sum(tab), tab[1,1],tab[2,1], tab[1,2], tab[2,2])
  names(res) <- c("total_n", do.call(paste0, expand.grid(outcome, c(0,1), "_", predictor, c(0,1))))
  
  if (0 %in% tab) {
    stats_empty <- rep(NA, length(types)*4)
    names(stats_empty) <- do.call(paste0, expand.grid(types,"_", c("estimate", "confint.lower", "confint.upper", "pval")) %>% arrange_all)
    res <- c(res, stats_empty)
    
  } else {
    
    for (type in types) {
      
      if (type == "unadjusted") {
        # unadjusted logistic regression
        fm <- formula(paste0(outcome, "~", predictor))
        fit <- glm(fm, data = data, family=family, weights = weights)
        
      } else if (type == "adjusted") {
        # adjusted logistic regression
        fm <- formula(paste0(outcome, "~", predictor, "+", paste0(covar_names, collapse = "+")))
        fit <- glm(fm, data = data, family=family, weights = weights)
        
      } else if (type == "conditional") {
        # conditional logistic regression
        fm <- formula(paste0(outcome, "~", predictor, "+", paste0(covar_names, collapse = "+"), " + strata(subclass)" ))
        fit <- clogit(fm, method="approximate", data = data, weights = weights)
        
      } else if (type == "marginal") {
        fm <- formula(paste0(outcome, "~", predictor, "+", paste0(covar_names, collapse = "+")))
        fit <- glm(fm, data = data, family=family, weights = weights)
        # marginal logistic regression
        fit <- avg_comparisons(fit, variables = predictor,
                               vcov = ~subclass,
                               newdata = subset(data, get(predictor) == 1), comparison = "lnoravg",
                               wts = "weights")
      }
      res <- c(res, StatsGet(fit, type))
    }
  }
  
  res <- data.frame(as.list(call.dat), as.list(res))
}


StatsGet <- function(fit, type) {
  
  if (type %in% c("marginal")) {
    summ <- summary(fit)
    sample.estimate <- summ$estimate
    lower.bound <- summ$conf.low
    upper.bound <- summ$conf.high
    pval = summ$p.value
    
  } else {
    
    if (type %in% c("unadjusted", "adjusted")) {
      idx = 2
    } else if (type %in% c("conditional")) {
      idx = 1
    }
    sample.estimate <- coef(fit)[idx]
    lower.bound <- suppressMessages(confint(fit)[idx,1])
    upper.bound <- suppressMessages(confint(fit)[idx,2])
    pval = summary(fit)$coefficients[idx, "Pr(>|z|)"]
  }
  
  # sample.se <- summary(fit)$coef[2,2]
  # degrees.freedom = fit$df.null
  # alpha = 0.05
  # t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
  # margin.error <- t.score * sample.se
  # lower.bound <- sample.estimate - margin.error
  # upper.bound <- sample.estimate + margin.error
  
  res <- data.frame(
    estimate = exp(sample.estimate),
    confint.lower = exp(lower.bound),
    confint.upper = exp(upper.bound),
    pval = pval
  )
  colnames(res) <- paste0(type, "_", colnames(res))
  
  return(res)
}

f=c(f1, f2)
f1 <- paste0("result ~ ", paste0(vars_selected, collapse = "+"))
dat="data1"
methods=c("nearest", "optimal", "full", "quick")
distances=c("glm", "gam", "gbm", "elasticnet", "randomforest", "cbps", "nnet", "bart", "robust_mahalanobis", "scaled_euclidean")
ratio=c(1,3,5)

data <- dataUnmatched
idx <- sample(nrow(data), 500)
data1 <- data[idx,]
seed=1

library(dplyr)
library(parallel)

# test with matching map all base
{
  timestart <- Sys.time()
  
  map <- expand.grid(seed=seed, f=f1, dat=dat, method=methods, distance=distances, stringsAsFactors = F, ratio = 1,
                     outcome="case", predictor="result", covar_names='c("sexMale", "age")')
  
  cl <- makeCluster(5, outfile="")
  clusterEvalQ(cl, {
    library(MatchIt)
    library(gdata)
    library(survival)
    library(marginaleffects)
    library(dplyr)
  })
  clusterExport(cl, c("data1", "MatchingFun", "StatsFun","StatsGet"))
  results <- clusterMap(cl, main, seed, formula_str=map$f, data_str=map$dat, method=map$method, distance=map$distance, ratio = map$ratio,
                        outcome=map$outcome, predictor=map$predictor, covar_names=map$covar_names)
  stopCluster(cl)
  
  res <- do.call(bind_rows, results)
  
  timeend <- Sys.time()
  timetaken <- timeend - timestart
  print(timetaken)
}
