rm(list = ls())

##########################################################
### Matching on lalonde dataset
##########################################################
lalonde <- MatchIt::lalonde

##########################################################
### MatchIt
##########################################################
library(MatchIt)

# No matching; constructing a pre-match matchit object
m.out0 <- matchit(
  formula(treat ~ age + educ + race + 
            + married + nodegree + re74 + re75, env = lalonde),
  data = data.frame(lalonde),
  method = NULL,
  # assess balance before matching
  distance = "glm" # logistic regression
)

# Checking balance prior to matching
summary(m.out0)

# 1:1 NN PS matching w/o replacement
m.out1 <- matchit(treat ~ age + educ,
                  data = lalonde,
                  method = "nearest",
                  distance = "glm")
summary(m.out1)

# examine visually
plot(m.out1, type = "jitter", interactive = FALSE)

##########################################################
### designMatch
##########################################################


##########################################################
### MatchingFrontier
##########################################################
