rm(list = ls())

library(fake)

set.seed(12345)

n=1000
ev_xy = 1

simul <- SimulateRegression(n=n, pk= 2, family = "gaussian", nu_xy = 1, ev_xy = ev_xy)

Heatmap(cor(simul$xdata),
        legend_range = c(-1, 1),
        col = c("navy", "white", "darkred")
)

plot(simul$xdata[,1], simul$ydata)

fit <- lm(simul$ydata ~ simul$xdata)
summary(fit)

Concordance(observed = simul$ydata, predicted = fit$fitted.values)

data.frame(simul$beta, fit$coefficients[-1])

j = 1
ypred <- simul$xdata %*% simul$beta[, j]
sigma <- sqrt((1 - ev_xy[j])/ev_xy[j] * stats::var(ypred))
ydata <- stats::rnorm(n = n, mean = ypred, 
                           sd = sigma)

Concordance(observed = simul$ydata, predicted = ydata)

ev_xy = 0.8
simul <- SimulateRegression(n=n, pk= 1, family = "gaussian", nu_xy = 1, ev_xy = ev_xy)
df <- data.frame(simul$xdata)

# binary treatment, could think about multinomial treatment
df$treat <- rbinom(n,1,0.5)
d_treat <- 1.5
j = 1
ypred <- simul$xdata %*% simul$beta[, j] + d_treat * df$treat
ypred2 <- simul$xdata[,1] * simul$beta[1, j] + simul$xdata[,2] * simul$beta[2, j] + d_treat * df$treat

sigma <- sqrt((1 - ev_xy[j])/ev_xy[j] * stats::var(ypred))
df$ydata <- stats::rnorm(n = n, mean = ypred, 
                      sd = sigma)


ggplot(df,
       aes(var1, ydata,
           color = as.factor(treat))) +
  geom_point() + 
  labs(color = "Treatment status")



fit <- lm(ydata ~ ., data = df)
summary(fit)
data.frame(c(simul$beta, d_treat), fit$coefficients[-1])

fit <- lm(ydata ~ treat)
summary(fit)

m.out <- matchit(treat ~ ., data = data.frame(simul$xdata), method = "quick", distance = "mahalanobis")
summary(m.out)
m.dat <- match.data(m.out)

fit <- lm(ydata ~ treat, weights = m.dat$weights)
summary(fit)

fit <- lm(ydata ~ simul$xdata + treat, weights = m.dat$weights)
summary(fit)
data.frame(c(simul$beta, d_treat), fit$coefficients[-1])

fit <- lmer(ydata ~ treat + (1 | subclass), data = m.dat, weights = m.dat$weights)
summary(fit)


#####
df <- simul1(1000)
ggplot(df,
       aes(x, y,
           color = as.factor(d))) +
  geom_point() + 
  labs(color = "Treatment status")

fit <- lm(y ~ d, data = df)
summary(fit)

#####
df <- tibble(
  # x is our confounder
  x = runif(n, -1, 4),
  # it affects the probability of receiving the treatment
  # in a NON-LINEAR way (step function)
  prob_d = ifelse(x > 0.5 & x < 2.5, 0.1, 0.9),
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
ggplot(df,
       aes(x, y,
           color = as.factor(d))) +
  geom_point() + 
  labs(color = "Treatment status")

fit <- lm(y ~ d, data = df)
summary(fit)


##### Confounder
# n = 1000 points for the simulation
n <- 1000

# create variables
# J is random draws from standard normal (mean = 0, stdev = 1)
confounding_var_J <- rnorm(n)

# J is used in creation of A since it is a cause of A (confounder)
independent_var_A <- 1.1 * confounding_var_J + rnorm(n)

# J is used in creation of X since it is a cause of X (confounder)
dependent_var_X <- 1.9 * confounding_var_J + 0.5 * independent_var_A + rnorm(n)

fit <- lm(dependent_var_X ~ independent_var_A)
summary(fit)

fit <- lm(dependent_var_X ~ independent_var_A + confounding_var_J)
summary(fit)

##### Treatment with Confounder

# create variables
# J is random draws from standard normal (mean = 0, stdev = 1)
confounding_var_J <- rnorm(n)

# J is used in creation of A since it is a cause of A (confounder)
prob <- rescale(confounding_var_J + rnorm(n))
independent_var_A <- rbinom(n, 1, prob)

# J is used in creation of X since it is a cause of X (confounder)
dependent_var_X <- 1.9 * confounding_var_J + 0.5 * independent_var_A + rnorm(n)

fit <- lm(dependent_var_X ~ independent_var_A)
summary(fit)

fit <- lm(dependent_var_X ~ independent_var_A + confounding_var_J)
summary(fit)

# Check if works
# and non linear transform

simul <- SimulateRegressionTreatment(n=1000)
df <- data.frame(y=simul$ydata[,1], treat=simul$treat, simul$xdata)
fit <- lm(y ~ treat, df)
summary(fit)

fit <- lm(y ~ ., df)
summary(fit)

fm <- as.formula(paste0("treat ~ ", paste0("var", c(1:10), collapse="+"), collapse=""))
m.out <- matchit(fm, df, method="quick", distance="mahalanobis")
summary(m.out)
m.dat <- match.data(m.out)

fit <- lm(y ~ treat, m.dat, weights = m.dat$weights)
summary(fit)

fm <- as.formula(paste0("y ~ treat + ", paste0("var", c(1:10), collapse="+"), collapse=""))
fit <- lm(fm, m.dat, weights = m.dat$weights)
summary(fit)

#### different expected number of treated
simul <- SimulateRegressionTreatment(n=10000, pk = 100, treat_p = 0.8, confounded = T, nu_conf = 0.5, nu_xy = 0.5)
df <- data.frame(y=simul$ydata[,1], treat=simul$treat, simul$xdata)
fit <- lm(y ~ treat, df)
summary(fit)

fit <- lm(y ~ ., df)
summary(fit)

fm <- as.formula(paste0("treat ~ ", paste0("var", c(1), collapse="+"), collapse=""))
m.out <- matchit(fm, df, method="quick", distance="mahalanobis", estimand = "ATE")
summary(m.out)
m.dat <- match.data(m.out)

fit <- lm(y ~ treat, m.dat, weights = m.dat$weights)
summary(fit)

fm <- as.formula(paste0("y ~ treat + ", paste0("var", c(1:10), collapse="+"), collapse=""))
fit <- lm(fm, m.dat, weights = m.dat$weights)
summary(fit)

#### non linear
set.seed(2)
start <- Sys.time()
simul <- SimulateRegressionTreatment(n=10000, pk = 1000, treat_p = 0.8, confounded = T, nu_conf = 0.5, nu_xy = 0.5, non_linear = "sq", nu_nl=1)
df <- data.frame(y=simul$ydata[,1], treat=simul$treat, simul$xdata)
end <- Sys.time()
print(end-start)

fit <- lm(y ~ treat, df)
summary(fit)

fit <- lm(y ~ ., df)
summary(fit)


### plot
ggplot(df,
       aes(x, y,
           color = as.factor(d))) +
  geom_point() + 
  labs(color = "Treatment status")

ggplot(df,
       aes(var1, y,
           color = as.factor(treat))) +
  geom_point() + 
  labs(color = "Treatment status")
