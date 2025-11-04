library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(caTools)  # For sample.split() for a stratified train/test split
library(flexsurv)
pak::pak("mlr-org/mlr3extralearners")
library(mlr3extralearners)
library(mlr3proba)
library(pec) # integrated brier score

# load data
dir_data <- paste0(here::here(), "/data")
datasets <- readRDS(file.path(dir_data, "datasets_sim_weibull.rds"))
df <- datasets[[1]][[1]]

set.seed(123)
split = sample.split(df$status, SplitRatio = 0.7)
train_set = subset(df, split == TRUE)
test_set  = subset(df, split == FALSE)

# specification needed for mlr3
train_task = as_task_surv(train_set, time = "time", event = "status")
test_task = as_task_surv(test_set, time = "time", event = "status")

# correct specification ------------------------------------
  ## base: flexsurv
    rg1 = flexsurvreg(
      formula = Surv(time, status) ~ x1 + x2 + x1:x2,
      anc = list(shape = ~ x2),
      data = train_set,
      dist = "weibullPH"
    )
  
  ## base: flexsurvspline, knots relevant for distributional assumption
    spl1 <- flexsurvspline(
      Surv(time, status) ~ x1 + x2 + x1:x2, 
      anc = list(gamma1 = ~ x2),
      data=train_set, 
      k=0, 
      scale="hazard"
    )
  
  # mlr3 implementation
    lrnr1 = lrn("surv.flexible", k = 0, anc = list(gamma1 = ~ x2), formula = Surv(time, status) ~ x1 + x2 + x1:x2)
    lrnr1$train(train_task)
    
# miss specification: no interaction ------------------------------------
  ## base: flexsurv
    rg2 = flexsurvreg(
      formula = Surv(time, status) ~ x1 + x2,
      anc = list(shape = ~ x2),
      data = train_set,
      dist = "weibullPH"
    )
    
  ## base: flexsurvspline, knots relevant for distributional assumption
    spl2 <- flexsurvspline(
      Surv(time, status) ~ x1 + x2, 
      anc = list(gamma1 = ~ x2),
      data=train_set, 
      k=0, 
      scale="hazard"
    )
    
  # mlr3 implementation
    lrnr2 = lrn("surv.flexible", k = 0, anc = list(gamma1 = ~ x2), formula = Surv(time, status) ~ x1 + x2)
    lrnr2$train(train_task)

# miss specification: lognormal ------------------------------------
  ## base: flexsurv
    rg3 = flexsurvreg(
      formula = Surv(time, status) ~ x1 + x2 + x1:x2,
      anc = list(sdlog = ~ x2),
      data = train_set,
      dist = "lnorm"
    )
    
  ## base: flexsurvspline, knots relevant for distributional assumption
    spl3 <- flexsurvspline(
      Surv(time, status) ~ x1 + x2 + x1:x2, 
      anc = list(gamma1 = ~ x2),
      data=train_set, 
      k=0, 
      scale="normal"
    )
    
  # mlr3 implementation
    lrnr3 = lrn(
      "surv.flexible", 
      k = 0,
      scale = "odds",
      anc = list(gamma1 = ~ x2), 
      formula = Surv(time, status) ~ x1 + x2 + x1:x2
    )
    lrnr3$train(train_task)
    
    
# Compare correct models
  ## flexsurvreg
    rg1_coefs = rg1$coefficients
    # transform to correct level
    rg1_coefs[names(rg1_coefs) == "shape(x2)"] = 
      exp(rg1_coefs[names(rg1_coefs) == "shape"] + rg1_coefs[names(rg1_coefs) == "shape(x2)"]) - 
      exp(rg1_coefs[names(rg1_coefs) == "shape"])
    rg1_coefs[names(rg1_coefs) == "shape"] = exp(rg1_coefs[names(rg1_coefs) == "shape"])
  
  ## flexsurvspline
    spl1_coefs = spl1$coefficients
    names(spl1_coefs) <- c("scale", "shape", "x1", "x2", "x1:x2", "shape(x2)")

  ## mlr3
    mlr1_coefs = learner$model$coefficients
    names(mlr1_coefs) <- c("scale", "shape", "x1", "x2", "x1:x2", "shape(x2)")

  terms <- names(rg1_coefs)

  df_coefs <- data.frame(
    term = terms,
    rg_cor = rg1_coefs[terms],
    spl_cor = spl1_coefs[terms],
    mlr_cor = mlr1_coefs[terms],
    row.names = NULL
  )

  print(df_coefs)
  
# Compare missspecified model: no intercept
  ## flexsurvreg
    rg2_coefs = rg2$coefficients
  # transform to correct level
    rg2_coefs[names(rg2_coefs) == "shape(x2)"] = 
      exp(rg2_coefs[names(rg2_coefs) == "shape"] + rg2_coefs[names(rg2_coefs) == "shape(x2)"]) - 
      exp(rg2_coefs[names(rg2_coefs) == "shape"])
    rg2_coefs[names(rg2_coefs) == "shape"] = exp(rg2_coefs[names(rg2_coefs) == "shape"])
  
  ## flexsurvspline
    spl2_coefs = spl2$coefficients
    names(spl2_coefs) <- c("scale", "shape", "x1", "x2", "shape(x2)")
  
  ## mlr3
    mlr2_coefs = lrnr2$model$coefficients
    names(mlr2_coefs) <- c("scale", "shape", "x1", "x2", "shape(x2)")
    
  terms <- names(rg2_coefs)
  
  df_coefs_noInt <- data.frame(
    term = terms,
    rg_noInt = rg2_coefs[terms],
    spl_noInt = spl2_coefs[terms],
    mlr_noInt = mlr2_coefs[terms],
    row.names = NULL
  )
  
  print(df_coefs_noInt)
  
  
# predictive performance
  ## concordance index
  
  rg1_p = predict(rg1, newdata = test_set, type = "lp")
  concordance(
    Surv(test_set$time, test_set$status) ~ rg1_p$.pred_link, reverse = TRUE
  )$concordance
  
  spl1_p = predict(spl1, newdata = test_set, type = "lp")
  concordance(
    Surv(test_set$time, test_set$status) ~ spl1_p$.pred_link, reverse = TRUE
  )$concordance
  
  mlr1_p = lrnr1$predict(test_task)
  mlr2_p = lrnr2$predict(test_task)
  mlr3_p = lrnr3$predict(test_task)
  
  ## ISBS / Graf score
  scores = rbind(
    mlr1_p$score(msr("surv.graf", p_max = 0.8)),
    mlr2_p$score(msr("surv.graf", p_max = 0.8)),
    mlr3_p$score(msr("surv.graf", p_max = 0.8))
  )

  print(scores)


