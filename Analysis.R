##################################################
## Project:
## Script purpose:
## Date: 2020-08-27
## Author: Daniel Smith
## Contact: dansmi@ceh.ac.uk
## Licence: MIT
##################################################


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(rstanarm)

# Read in this data (Broad hbaitat here) ----------------------------------

data <- read_csv("Oviposition_Summary_Hab1.csv")

# Ordinary GLM -----------------------------------------------------------

# How does habitat type influence vector status?

# glm version using base R 
glm1 <- glm(formula = vector ~ treehole + rockhole + ground + pool + marsh + swamp, 
            data = data, family = binomial())

# Though the diagnostic plots aren't looking great:
plot(glm1)


# A Bayes RStanarm glm ----------------------------------------------------

# Formula design  for just habitat predictors
habs <- data %>% select(contains("hab")) %>% colnames()

# So, can a vector be predicted by all/any habitats?
habformula <- paste("vector ~ ", paste(habs, collapse = " + ")) %>% as.formula()


# The first GLM model - logistic regression with the default priors used in rstanarm
m1 <- stan_glm(formula = habformula,
               family = binomial(link = "logit"),
               data = data,
               cores = 4, seed = 1234)

# Now criticise the fit of the model using shinystan
launch_shinystan(m1)

##################################################################################
##  Criticising the model is complex and requires a good understanding of how   ##
##  the sampling of joint distributions in Stan and Bayes in general works.     ##
##  This is important to understand for good model outputs, essentially, the    ##
##  chains need to converge (Rhat < 1.01), there needs to be few divergent      ##
##  transitions, and ideally roughly 10% of the number of effective samples     ##
##  (So if chains ran for 2000 iterations, you would expect at least 200).      ##
##  Sometimes these are incredibly difficult to fit, using RStanarm means the   ##
##  statistications who wrote Stan have optimised these functions to work well  ##
##  most times. You should check out the rstanarm help field and website        ##
##  (https://mc-stan.org/rstanarm/articles/rstanarm.html) if you want to know   ##
##  more                                                                        ##
##################################################################################

# This gives a plot of the posterior distribution for our model
# (The sum of the prior expectations and the data)
plot(m1, "areas", prob = 0.95, prob_outer = 1)

# Looks like habitat E & A are positive predictors of vectorability:
# (It's credible intervals do not overlap with 0)
posterior_interval(m1, prob = 0.95)

# These are the default priors after covariates are centered 
# (makes joint distribution easier to estimate, and is done automatically 
# in rstanarm and most other stan packages)
prior_summary(m1)

# We can see how our priors look and influence the data by drawing a model just using them and no data:
priorplot <- stan_glm(formula = habformula,
               family = binomial(link = "logit"),
               data = data,
               cores = 4,
               # Here set to true to only draw conclusions from the prior distributions
               prior_PD = TRUE)

# Now look at only the priors. These make sense, we don't expect any of these 
# values to be far from 1 or 0 in real life or have a huge impact on the outcome
# considering how complex the chances of being a vector might be and what it 
# relies on is likely much more than what habitats mosquitoes prefer

plot(priorplot, "areas", prob = 0.95, prob_outer = 1)

### What about the other variables that aren't EUNIS habitat related?

# Formula design for no habitats
nohabs <- data %>% 
  select(-contains("hab"), -mosquito_species, -vector) %>% 
  colnames()

# So, is Vector predicted by any of the other characteristics we recorded?
nohabformula <- paste("vector ~ ", paste(nohabs, collapse = " + ")) %>% as.formula()

# The second GLM model - logistic regression with the default priors used in rstanarm
m2 <- stan_glm(formula = nohabformula,
               family = binomial(link = "logit"),
               data = data,
               cores = 4)

# Now criticise the fit of the model using shinystan
launch_shinystan(m2)

# Look
plot(m2, prob = 0.95, prob_outer = 1)

# Looks like habitat E & A are positive predictors of vectorability:
# It's credible intervals do not overlap with 0. 
posterior_interval(m2, prob = 0.95)

### Another Example: Can we explain what makes a mosquito invasive?

Introformula <- paste("Introduced ~ ", paste(habs, collapse = " + ")) %>% as.formula()

m3 <- stan_glm(formula = Introformula,
               family = binomial(link = "logit"),
               data = data,
               cores = 4)
plot(m3, prob = 0.95, prob_outer = 1)

### Another example on the more specific level 2 habitat characteristics:

datahab2 <- read_csv("Oviposition_Summary_Hab2.csv")

habs2 <- datahab2 %>% select(contains("hab")) %>% colnames()

# Use 2nd level habs?
Introformulahabs2 <- paste("Introduced ~ ", paste(habs2, collapse = " + ")) %>% as.formula()

m4 <- stan_glm(formula = Introformulahabs2,
               family = binomial(link = "logit"),
               data = datahab2,
               cores = 4)

# Remember to criticise the model

plot(m4, prob = 0.95, prob_outer = 1)

#################################################################################
##  A GBM machine model would be interesting and an example is below. But I    ##
##  think in any case the sample size is too small for this method. If we had  ##
##  2-3x as many mosquitoes it would be a nice approach.                       ##
#################################################################################

# https://xgboost.readthedocs.io/en/latest/R-package/xgboostPresentation.html

sample <- sample.int(n = nrow(data), size = floor(.50 * nrow(data)), replace = F)
train <- data[sample,]
test <- data[-sample,]

library(xgboost)
library(Matrix)

sparse_mat <- sparse.model.matrix(vector ~., -1, data = train)

output_vector = train$vector

output_vector = train[,"vector"] == 1

bst <- xgboost(data = data.matrix(train), label = output_vector, max.depth = 4,
               eta = 1, nthread = 10, nrounds = 100, objective = "binary:logistic", verbose = 0)

pred <- predict(bst, data.matrix(test))

prediction <- as.numeric(pred > 0.5)
print(head(prediction))

err <- mean(as.numeric(pred > 0.5) != test$vector)
print(paste("test-error=", err))

##################################################################################
##  LASSO Regression would be a similar case. The sample size is just too       ##
##  small for anything meaningful. But I've included a link that explains what  ##
##  they are below.                                                             ##
##################################################################################

# http://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/


