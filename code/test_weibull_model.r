library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(caTools)  # For sample.split() for a stratified train/test split
library(flexsurv)
# pak::pak("mlr-org/mlr3extralearners")
library(mlr3extralearners)
library(mlr3proba)

# load data
dir_data <- "C:/Users/ra56yaf/Downloads/data_sim"
datasets <- readRDS(file.path(dir_data, "datasets_sim_weibull.rds"))
df <- datasets[[2]][[1]]

set.seed(123)
split = sample.split(df$status, SplitRatio = 0.7)
train_set = subset(df, split == TRUE)
test_set  = subset(df, split == FALSE)

# standard implemetnation
mod1 = flexsurvreg(
  formula = Surv(time, status) ~ x1 + x2 + x1:x2,
  anc = list(shape = ~ x2),
  data = train_set,
  dist = "weibullPH"
)
p1 = predict(mod1, newdata = test_set, type = "lp")
concordance(
  Surv(test_set$time, test_set$status) ~ p1$.pred_link, reverse = TRUE
)$concordance

mod2 = flexsurvspline(
  formula = Surv(time, status) ~ x1 + x2 + x1:x2,
  k = 0,
  scale = "hazard",
  anc = list(gamma1 = ~ x2),
  data = train_set
)
p2 = predict(mod2, newdata = test_set, type = "lp")
concordance(
  Surv(test_set$time, test_set$status) ~ p2$.pred_link, reverse = TRUE
)$concordance

# mlr3 implementation
train_task = as_task_surv(train_set, time = "time", event = "status")
test_task = as_task_surv(test_set, time = "time", event = "status")

learner = lrn("surv.flexible", k = 0, anc = list(gamma1 = ~ x2), formula = Surv(time, status) ~ x1 + x2 + x1:x2)
learner$train(train_task)
learner$model

p = learner$predict(test_task)
p$score()