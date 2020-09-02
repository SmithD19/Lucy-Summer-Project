library(tidyverse)
library(tidybayes)
library(brms)
library(haven)

ThaiEdu_Raw <- haven::read_sav("https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/blob/master/chapter%206/Thaieduc/thaieduc.sav?raw=true")
head(ThaiEdu_Raw)

ThaiEdu_New <- ThaiEdu_Raw %>%
  mutate(SCHOOLID = factor(SCHOOLID),
         SEX = if_else(SEX == 0, "girl", "boy"),
         SEX = factor(SEX, levels = c("girl", "boy")),
         PPED = if_else(PPED == 0, "no", "yes"),
         PPED = factor(PPED, levels = c("no", "yes")))

head(ThaiEdu_New)

ThaiEdu_New %>%
  summarise_each(list(~sum(is.na(.)))) %>%
  gather()

ThaiEdu_New <- ThaiEdu_New %>% na.omit()


Bayes_Model_Binary <- brm(formula = REPEAT ~ SEX * PPED,  
                          data=ThaiEdu_New, 
                          family = bernoulli(link = "logit"),
                          chains = 4, 
                          inits = "0", 
                          cores = 4,
                          seed = 123)

mcmc_plot(Bayes_Model_Binary, type = "areas", prob = 0.95)

exp(fixef(Bayes_Model_Binary))

conditional_effects(Bayes_Model_Binary)

Pred <- predict(Bayes_Model_Binary, type = "response")
Pred <- if_else(Pred[,1] > 0.5, 1, 0)

ConfusionMatrix <- table(Pred, pull(ThaiEdu_New, REPEAT))

sum(diag(ConfusionMatrix))/sum(ConfusionMatrix)

# Compute AUC for predicting Class with the model
library(ROCR)
Prob <- predict(Bayes_Model_Binary, type="response")
Prob <- Prob[,1]
Pred <- prediction(Prob, as.vector(pull(ThaiEdu_New, REPEAT)))
AUC <- performance(Pred, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC


# transform data ----------------------------------------------------------

ThaiEdu_Prop <- ThaiEdu_New %>%
  group_by(SCHOOLID, MSESC) %>%
  summarise(REPEAT = sum(REPEAT),
            TOTAL = n()) %>%
  ungroup()

head(ThaiEdu_Prop)

Bayes_Model_Prop <- brm(REPEAT | trials(TOTAL) ~ MSESC,
                        data = ThaiEdu_Prop, 
                        family = binomial(link = "logit"),
                        warmup = 500, 
                        iter = 2000, 
                        chains = 2, 
                        inits = "0", 
                        cores = 2,
                        seed = 123)

conditional_effects(Bayes_Model_Prop)
