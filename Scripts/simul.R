####
# Simulation 1
# simple linear confounder
####
library(dplyr)

simul1 <- function(n) {
  # Data-generating code
  df <- tibble(
    # x is our confounder
    x = runif(n, -1, 4),
    # it affects the probability of receiving the treatment
    # in a NON-LINEAR way (step function)
    # prob_d = ifelse(x > 0.5 & x < 2.5, 0.1, 0.9),
    prob_d = 0.5,
    d = rbinom(n, 1, prob_d),
    noise = rnorm(n, sd = 0.5),
    # for simplification, the treatment effect is homogeneous
    treat_effect = 1,
    # x also effects the outcome in a non-linear way
    # y = sin(x) + d*treat_effect + noise
    y = 2 * x + d*treat_effect + noise
  ) %>% 
    mutate(d_factor = factor(d,
                             levels=c(0,1), labels=c("Untreated",
                                                     "Treated")))
  return(df)
}

####
# Simulation 2
# non-linear confounder
####

simul2 <- function(n) {
  # Data-generating code
  df <- tibble(
    # x is our confounder
    x = runif(n, -1, 4),
    # it affects the probabiity of receiving the treatment
    # in a NON-LINEAR way (step function)
    prob_d = 0.5,
    #prob_d = ifelse(x > 0.5 & x < 2.5, 0.1, 0.9),
    d = rbinom(n, 1, prob_d),
    noise = rnorm(n, sd = 0.1),
    # for simplification, the treatment effect is homogeneous
    treat_effect = 1,
    # x also effects the outcome in a non-linear way
    y = sin(x) + d*treat_effect + noise
    # y = 2 * x + d*treat_effect + noise
  ) %>% 
    mutate(d_factor = factor(d,
                             levels=c(0,1), labels=c("Untreated",
                                                     "Treated")))
  return(df)
}

####
# Simulation 3
# Additive
####
# just pluses, multivariate

simul3 <- function(n) {
  # Data-generating code
  df <- tibble(
    # x is our confounder
    x1 = runif(n, -1, 4),
    x2 = runif(n, -1, 4),
    # it affects the probabiity of receiving the treatment
    # in a NON-LINEAR way (step function)
    prob_d = 0.5,
    d = rbinom(n, 1, prob_d),
    noise = rnorm(n, sd = 0.1),
    # for simplification, the treatment effect is homogeneous
    treat_effect = 1,
    # x also effects the outcome in a non-linear way
    y = 1 * x1 + 2 * x2 + d*treat_effect + noise
    # y = 2 * x + d*treat_effect + noise
  ) %>% 
    mutate(d_factor = factor(d,
                             levels=c(0,1), labels=c("Untreated",
                                                     "Treated")))
  return(df)
}

df <- simul3(1000)


####
# Simulation 4
# Interaction
####

# essentially need to check whether matching works here

# function to generate data, fit correct model, and extract p-values of
# interaction 
f1 <- function(n = 1000){
  # generate data
  g <- gl(n = 2, k = n/2)
  x <- rnorm(n, mean = 10, sd = 2)
  y <- 1.2 + 0.8*(g == "2") + 0.6*x + -0.5*x*(g == "2") + rnorm(n, sd = 1.1)
  d <- data.frame(y, g, x)
  # fit correct model
  m1 <- lm(y ~ g + x + g:x, data = d)
  sm1 <- summary(m1)
  sm2 <- coeftest(m1, vcov. = vcovHC(m1))
  # get p-values using usual SE and robust ES
  c(usual = sm1$coefficients[4,4], 
    robust = sm2[4,4])
}

# Some code for plotting and replicate function
# ggplot(d,
#        aes(x, y,
#            color = g)) +
#   geom_point() + 
#   labs(color = "Treatment status")
# 
# # run the function 1000 times
# r_out <- replicate(n = 1000, expr = f1())
# 
# # get proportion of times we correctly reject Null of no interaction at p < 0.05
# # (ie, estimate power)
# apply(r_out, 1, function(x)mean(x < 0.05))

####
# Simulation 4
# additional confounding (e.g. time-variants)
####

# change probability of assignment
# step-wise probability

### MSE/RMSE of beta
### proportion of pvalue

###################################################
# NEW SimulateRegression FG
# add treatment and confounding
###################################################
# n = 100
# pk = 10
# xdata = NULL
# family = "gaussian" 
# q = 1
# theta = NULL
# nu_xy = 0.2
# beta_abs = c(0.1, 1)
# beta_sign = c(-1, 1)
# continuous = TRUE
# ev_xy = 0.7
# nu_conf = 0.5
# complexity = NULL
# nu_nl = 0.2

attach(getNamespace("fake"))

SimulateRegressionTreatment <- function (n = 100, pk = 10, xdata = NULL, family = "gaussian", 
          q = 1, theta = NULL, nu_xy = 0.2, beta_abs = c(0.1, 1), 
          beta_sign = c(-1, 1), continuous = TRUE, ev_xy = 0.7, nu_conf = 0.2, treat_p = 0.5, treat_beta = 1, complexity = "linear", nu_nl = 0.2) 
{
  if (is.null(xdata)) {
    if (is.null(pk) | is.null(n)) {
      stop("Argument 'xdata' must be provided if 'pk' and 'n' are not provided.")
    }
  }
  if (length(ev_xy) != q) {
    ev_xy <- rep(ev_xy[1], q)
  }
  if (length(nu_xy) != q) {
    nu_xy <- rep(nu_xy[1], q)
  }
  if (!is.null(xdata)) {
    n <- nrow(xdata)
    p <- ncol(xdata)
  } else {
    p <- sum(pk)
    xsimul <- SimulateGraphical(n = n, pk = pk, theta = NULL, 
                                implementation = HugeAdjacency, topology = "random", 
                                nu_within = 0, nu_between = 0, nu_mat = NULL, v_within = 0, 
                                v_between = 0, v_sign = c(-1, 1), continuous = TRUE, 
                                pd_strategy = "diagonally_dominant", ev_xx = NULL, 
                                scale_ev = TRUE, u_list = c(1e-10, 1), tol = .Machine$double.eps^0.25, 
                                scale = TRUE, output_matrices = FALSE)
    xdata <- xsimul$data
  }
  if (!is.null(theta)) {
    if (q == 1) {
      if (is.vector(theta)) {
        theta <- cbind(theta)
      }
    }
    if (ncol(theta) != q) {
      stop("Arguments 'theta' and 'q' are not compatible. Please provide a matrix 'theta' with 'q' columns.")
    }
    if (nrow(theta) != p) {
      stop("Please provide a matrix 'theta' with as many columns as predictors.")
    }
    theta <- ifelse(theta != 0, yes = 1, no = 0)
  }
  if (is.null(theta)) {
    theta <- SamplePredictors(pk = p, q = q, nu = nu_xy, 
                              orthogonal = FALSE)
  }

  beta <- theta
  if (continuous) {
    beta <- beta * matrix(stats::runif(n = nrow(beta) * 
                                         ncol(beta), min = min(beta_abs), max = max(beta_abs)), 
                          nrow = nrow(beta), ncol = ncol(beta))
  } else {
    beta <- beta * matrix(base::sample(beta_abs, size = nrow(beta) * 
                                         ncol(beta), replace = TRUE), nrow = nrow(beta), 
                          ncol = ncol(beta))
  }
  beta <- beta * matrix(base::sample(beta_sign, size = nrow(beta) * 
                                       ncol(beta), replace = TRUE), nrow = nrow(beta), ncol = ncol(beta))
  
  if (nu_conf != 0) {
    # Confounding + treatment
    theta_conf <- SamplePredictors(pk = p, q = q, nu = nu_conf, 
                                   orthogonal = FALSE)
    theta_conf_prob <- theta_conf * matrix(stats::runif(n= nrow(theta_conf) * ncol(theta_conf), min = 0, max = 1))
    treat_prob <- rescale(xdata %*% theta_conf_prob)
    treat_prob <- treat_prob^log(treat_p, mean(treat_prob)) # transform to treat_p
  } else {
    #generate warning to ignore nu_conf
    theta_conf <- matrix(0,ncol = q, nrow=p)
    treat_prob = treat_p
  }
  
  treat <- rbinom(n, 1, treat_prob)
  
  ydata <- matrix(NA, ncol = q, nrow = nrow(xdata))
  if (family == "gaussian") {
    for (j in 1:q) {
      ypred <- xdata %*% beta[, j]
      
      if(complexity != "linear") {
        theta_nl <- SamplePredictors(pk = p, q = q, nu = nu_nl, 
                                     orthogonal = F)
        if(complexity == "sq") {
          ypred[,theta_nl] = ypred[,theta_nl] ^ 2
        } else if (complexity == "sin") {
          ypred[,theta_nl] = sin(ypred[,theta_nl] * pi)
        } else {
          stop('Argument "complexity" only supports "sq" or "sin"')
        }
      } else {
        theta_nl <- matrix(0,ncol = q, nrow=p)
      }

      ypred = ypred + treat_beta * treat
      sigma <- sqrt((1 - ev_xy[j])/ev_xy[j] * stats::var(ypred))
      ydata[, j] <- stats::rnorm(n = n, mean = ypred, 
                                 sd = sigma)
    }
  }
  if (family == "binomial") {
    for (j in 1:q) {
      crude_log_odds <- xdata %*% beta[, j] + treat_beta * treat
      s_max <- max(abs(crude_log_odds))/log(0.51/0.49)
      s_min <- min(abs(crude_log_odds))/log(0.99/0.01)
      argmax_scaling_factor <- stats::optimise(f = TuneCStatisticLogit, 
                                               crude_log_odds = crude_log_odds, auc = ev_xy[j], 
                                               lower = 1/s_max, upper = 1/s_min)
      scaling_factor <- argmax_scaling_factor$minimum
      beta[, j] <- beta[, j] * scaling_factor
      log_odds <- crude_log_odds * scaling_factor
      proba <- 1/(1 + exp(-log_odds))
      ydata[, j] <- stats::rbinom(n = n, size = 1, prob = proba)
    }
  }
  rownames(ydata) <- rownames(xdata)
  colnames(ydata) <- paste0("outcome", 1:q)
  rownames(beta) <- rownames(theta) <- colnames(xdata)
  colnames(beta) <- colnames(theta) <- colnames(ydata)
  out <- list(xdata = xdata, ydata = ydata, beta = beta, theta = theta, treat= treat, treat_beta=treat_beta, theta_conf = theta_conf, theta_nl = theta_nl)
  class(out) <- "simulation_regression"
  return(out)
}
