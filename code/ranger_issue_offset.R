## ========================================================================== ##
# COMPARE RSFC AND RANGER FOR COMPETING RISKS
# OBJECTIVE: Ranger and rsfc produced very different results.

# Structure:
# 0. INITIALIZE FIGURE + FUNCTIONS FOR AJ
# 1. RFSC
# 2. RANGER
# 3. RANGER + OFFSET
#
## ========================================================================== ##

# load libraries
library(dplyr)
library(tidyr)
library(mgcv)
library(pammtools)
library(ggplot2)
library(mvna)              # pneunomia data set
library(etm)               # aalen johansen estimator
library(cmprsk)            # needed for competing risks in etm, cf. Beyersmann p. 79
library(pseudo)
library(ranger)
library(mlr3)
library(mlr3learners)
library(mlr3pipelines)
library(mlr3proba)
library(mlr3extralearners) # for randomForestSRC learner / pak::pak("mlr-org/mlr3extralearners")
theme_set(theme_bw())

## ========================================================================== ##
## 0. Initialize Figure + Functions. Make sure to adapt relative path if needed
## ========================================================================== ##
source(here::here("code/functions/etm-ci-trafo.R"))

## initialize variables for plotting
# fontsize
label_size = 24
headline_size = 24

# initialize lines and dots
linewidth = 1
pointsize = 2
strokewidth = 1.5
# model_colors <- c("pam" = "firebrick2", "dt" = "steelblue", "pv" = "springgreen4", "aj" = "black")
model_colors <- c(
  "rsfc" = "#D55E00",   # vivid reddish-orange
  "ranger"  = "#0072B2",   # deep sky blue
  "ranger_offset"  = "#009E73",   # bluish green (unchanged)
  "aj"  = "#000000"    # black (unchanged)
)

model_fills  <- rep("darkgrey", 4); names(model_fills) <- names(model_colors)

theme_cif <- theme(
  text = element_text(size = 10),
  axis.title = element_text(size = 11),
  axis.text = element_text(size = 9),
  legend.position = "bottom",
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),
  strip.text = element_text(size = 11, face = "bold")
)

# load pneunomia data
data(sir.adm, package = "mvna")

## ========================================================================== ##
## 1. "Benchmark" Aalen Johansen incl. confidence intervals
## ========================================================================== ##
tra <- matrix(FALSE, ncol = 3, nrow = 3)
dimnames(tra) <- list(c("0", "1", "2"), c("0", "1", "2"))
tra[1, 2:3] <- TRUE

to <- ifelse(sir.adm$status==0
             , "cens"
             , ifelse(sir.adm$status==1,2,1)
)
my.sir.data <- data.frame(id=sir.adm$id
                          , from=0
                          , to
                          , time=sir.adm$time
                          , pneu=sir.adm$pneu
)

my.nelaal.nop <- mvna(my.sir.data[my.sir.data$pneu == 0, ]
                      , c("0", "1", "2"), tra, "cens")
my.nelaal.p <- mvna(my.sir.data[my.sir.data$pneu == 1, ]
                    , c("0", "1", "2"), tra, "cens")

my.sir.cif <- cuminc(my.sir.data$time
                     , my.sir.data$to
                     , group=my.sir.data$pneu
                     , cencode="cens")

my.sir.etm.nop <- etm(my.sir.data[my.sir.data$pneu == 0, ], c("0", "1", "2"), tra, "cens", s = 0)
my.sir.etm.p <- etm(my.sir.data[my.sir.data$pneu == 1, ], c("0", "1", "2"), tra, "cens", s = 0)

# use internal etm function to generate aj with confidence interval
ls_aj_out <- list(nop = list(ci.transfo(my.sir.etm.nop, tr.choice = '0 1', transfo = "cloglog")
                             , ci.transfo(my.sir.etm.nop, tr.choice = '0 2', transfo = "cloglog"))
                  , p = list(ci.transfo(my.sir.etm.p, tr.choice = '0 1', transfo = "cloglog")
                             , ci.transfo(my.sir.etm.p, tr.choice = '0 2', transfo = "cloglog"))
)

# build data frame for plotting
ndf_cr_aj <- rbind(
  data.frame(tend = ls_aj_out$p[[1]]$`0 1`$time
             , cif = ls_aj_out$p[[1]]$`0 1`$P
             , cif_lower = ls_aj_out$p[[1]]$`0 1`$lower
             , cif_upper = ls_aj_out$p[[1]]$`0 1`$upper
             , cause = 1
             , pneu = 0),
  data.frame(tend = ls_aj_out$p[[2]]$`0 2`$time
             , cif = ls_aj_out$p[[2]]$`0 2`$P
             , cif_lower = ls_aj_out$p[[2]]$`0 2`$lower
             , cif_upper = ls_aj_out$p[[2]]$`0 2`$upper
             , cause = 2
             , pneu = 0),
  data.frame(tend = ls_aj_out$nop[[1]]$`0 1`$time
             , cif = ls_aj_out$nop[[1]]$`0 1`$P
             , cif_lower = ls_aj_out$nop[[1]]$`0 1`$lower
             , cif_upper = ls_aj_out$nop[[1]]$`0 1`$upper
             , cause = 1
             , pneu = 1),
  data.frame(tend = ls_aj_out$nop[[2]]$`0 2`$time
             , cif = ls_aj_out$nop[[2]]$`0 2`$P
             , cif_lower = ls_aj_out$nop[[2]]$`0 2`$lower
             , cif_upper = ls_aj_out$nop[[2]]$`0 2`$upper
             , cause = 2
             , pneu = 1)
) |>
  mutate(
    cause = factor(cause, labels = c("Death", "Discharge"))
    , pneu = factor(pneu, labels = c("Pneumonia", "No Pneumonia"))
    , model = "aj"
  )

## ========================================================================== ##
## 2. Random Forest Classifier RSFC
## ========================================================================== ##

# # define equidistant cut points
eventtimes <- unique(sort(sir.adm[sir.adm$status != 0, "time"]))
cut <- sort(union(seq(from=1, to=150, by=1), eventtimes))
# # -> no effect

# create ped data set
# sir.adm$status <- as.factor(sir.adm$status)
ped_cr <- as_ped(sir.adm, Surv(time, status)~ ., combine = TRUE, max_time = 150
                 , cut = cut
) |>
  mutate(cause = as.factor(cause), pneu = as.factor(pneu))

ped_cr_rf <- ped_cr |>
  mutate(ped_status = as.factor(ped_status))

tsk_pneu = TaskClassif$new(
  id = "pneu", 
  target = "ped_status",
  backend = select(ped_cr_rf, ped_status, tend, cause, pneu)) # offset makes the fit good.

## include RF hazard calculation and prediction
lrn_rf_dt = po("encode", method = "treatment") %>>%
  lrn("classif.rfsrc", ntree=1000L) |> as_learner()

lrn_rf_dt$predict_type <- "prob"

set.seed(32168)
# Train the pipeline
lrn_rf_dt$train(tsk_pneu)

ndf_cr <- ped_cr_rf |>
  make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause))

haz_rf_dt = lrn_rf_dt$predict_newdata(ndf_cr)$prob |> as.data.table()
ndf_cr_dt_interim <- ndf_cr |> cbind(hazard = haz_rf_dt$'1')

ndf_cr_dt_wide <- ndf_cr_dt_interim %>%
  pivot_wider(names_from = cause, values_from = hazard, names_prefix = "hazard_", values_fill = 0) %>% # so each row contains all cause-specific hazards
  mutate(hazard_allCause = rowSums(across(starts_with("hazard_")))) %>% # sum of all cause-specific hazards
  arrange(pneu, tend) %>%
  group_by(pneu) %>%
  mutate(
    surv_allCause = cumprod(1 - hazard_allCause), # S_allCause(t) = product over u ≤ t of (1 - hazard_allCause(u))
    surv_allCause_lag = lag(surv_allCause, default = 1),
    cif_1 = cumsum(surv_allCause_lag * hazard_1), # CIF_j(t) = sum over u ≤ t of [S_allCause(u^-)*hazard_j(u)]
    cif_2 = cumsum(surv_allCause_lag * hazard_2)
  ) %>%
  ungroup()

ndf_cr_rfsc <- ndf_cr_dt_wide %>%
  select(tend, pneu, cif_1, cif_2) %>%
  pivot_longer(
    cols = starts_with("cif_"),
    names_to = "cause",
    names_prefix = "cif_",
    values_to = "cif"
  ) %>%
  mutate(
    cif_lower = NA,
    cif_upper = NA,
    cause = factor(cause, levels = c("1", "2"), labels = c("Discharge", "Death")),
    pneu = factor(pneu, labels = c("No Pneumonia", "Pneumonia")),
    model = "rsfc"
  )

## ========================================================================== ##
## 3. Random Forest Classifier // Ranger
## ========================================================================== ##

# # define equidistant cut points
eventtimes <- unique(sort(sir.adm[sir.adm$status != 0, "time"]))
cut <- sort(union(seq(from=1, to=150, by=1), eventtimes))
# # -> no effect

# create ped data set
# sir.adm$status <- as.factor(sir.adm$status)
ped_cr <- as_ped(sir.adm, Surv(time, status)~ ., combine = TRUE, max_time = 150
                 , cut = cut
) |>
  mutate(cause = as.factor(cause), pneu = as.factor(pneu))

ped_cr_rf <- ped_cr |>
  mutate(ped_status = as.factor(ped_status))

tsk_pneu = TaskClassif$new(
  id = "pneu", 
  target = "ped_status",
  backend = select(ped_cr_rf, ped_status, tend, cause, pneu)) # offset makes the fit good.

## include RF hazard calculation and prediction
lrn_rf_dt = po("encode", method = "treatment") %>>%
  lrn("classif.ranger") |> as_learner()

lrn_rf_dt$predict_type <- "prob"

set.seed(32168)
# Train the pipeline
lrn_rf_dt$train(tsk_pneu)

ndf_cr <- ped_cr_rf |>
  make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause))

haz_rf_dt = lrn_rf_dt$predict_newdata(ndf_cr)$prob |> as.data.table()
ndf_cr_dt_interim <- ndf_cr |> cbind(hazard = haz_rf_dt$'1')

ndf_cr_dt_wide <- ndf_cr_dt_interim %>%
  pivot_wider(names_from = cause, values_from = hazard, names_prefix = "hazard_", values_fill = 0) %>% # so each row contains all cause-specific hazards
  mutate(hazard_allCause = rowSums(across(starts_with("hazard_")))) %>% # sum of all cause-specific hazards
  arrange(pneu, tend) %>%
  group_by(pneu) %>%
  mutate(
    surv_allCause = cumprod(1 - hazard_allCause), # S_allCause(t) = product over u ≤ t of (1 - hazard_allCause(u))
    surv_allCause_lag = lag(surv_allCause, default = 1),
    cif_1 = cumsum(surv_allCause_lag * hazard_1), # CIF_j(t) = sum over u ≤ t of [S_allCause(u^-)*hazard_j(u)]
    cif_2 = cumsum(surv_allCause_lag * hazard_2)
  ) %>%
  ungroup()

ndf_cr_ranger <- ndf_cr_dt_wide %>%
  select(tend, pneu, cif_1, cif_2) %>%
  pivot_longer(
    cols = starts_with("cif_"),
    names_to = "cause",
    names_prefix = "cif_",
    values_to = "cif"
  ) %>%
  mutate(
    cif_lower = NA,
    cif_upper = NA,
    cause = factor(cause, levels = c("1", "2"), labels = c("Discharge", "Death")),
    pneu = factor(pneu, labels = c("No Pneumonia", "Pneumonia")),
    model = "ranger"
  )

## ========================================================================== ##
## 4. Random Forest Classifier // Ranger + offset (always 0)
## ========================================================================== ##

# # define equidistant cut points
eventtimes <- unique(sort(sir.adm[sir.adm$status != 0, "time"]))
cut <- sort(union(seq(from=1, to=150, by=1), eventtimes))
# # -> no effect

# create ped data set
# sir.adm$status <- as.factor(sir.adm$status)
ped_cr <- as_ped(sir.adm, Surv(time, status)~ ., combine = TRUE, max_time = 150
                 , cut = cut
) |>
  mutate(cause = as.factor(cause), pneu = as.factor(pneu))

ped_cr_rf <- ped_cr |>
  mutate(ped_status = as.factor(ped_status))

summary(ped_cr_rf$offset) # check that offset carries no information

tsk_pneu = TaskClassif$new(
  id = "pneu", 
  target = "ped_status",
  backend = select(ped_cr_rf, ped_status, tend, cause, pneu, offset)) # offset makes the fit good.

## include RF hazard calculation and prediction
lrn_rf_dt = po("encode", method = "treatment") %>>%
  lrn("classif.ranger") |> as_learner()

lrn_rf_dt$predict_type <- "prob"

set.seed(32168)
# Train the pipeline
lrn_rf_dt$train(tsk_pneu)

ndf_cr <- ped_cr_rf |>
  make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause))

haz_rf_dt = lrn_rf_dt$predict_newdata(ndf_cr)$prob |> as.data.table()
ndf_cr_dt_interim <- ndf_cr |> cbind(hazard = haz_rf_dt$'1')

ndf_cr_dt_wide <- ndf_cr_dt_interim %>%
  pivot_wider(names_from = cause, values_from = hazard, names_prefix = "hazard_", values_fill = 0) %>% # so each row contains all cause-specific hazards
  mutate(hazard_allCause = rowSums(across(starts_with("hazard_")))) %>% # sum of all cause-specific hazards
  arrange(pneu, tend) %>%
  group_by(pneu) %>%
  mutate(
    surv_allCause = cumprod(1 - hazard_allCause), # S_allCause(t) = product over u ≤ t of (1 - hazard_allCause(u))
    surv_allCause_lag = lag(surv_allCause, default = 1),
    cif_1 = cumsum(surv_allCause_lag * hazard_1), # CIF_j(t) = sum over u ≤ t of [S_allCause(u^-)*hazard_j(u)]
    cif_2 = cumsum(surv_allCause_lag * hazard_2)
  ) %>%
  ungroup()

ndf_cr_ranger_offset <- ndf_cr_dt_wide %>%
  select(tend, pneu, cif_1, cif_2) %>%
  pivot_longer(
    cols = starts_with("cif_"),
    names_to = "cause",
    names_prefix = "cif_",
    values_to = "cif"
  ) %>%
  mutate(
    cif_lower = NA,
    cif_upper = NA,
    cause = factor(cause, levels = c("1", "2"), labels = c("Discharge", "Death")),
    pneu = factor(pneu, labels = c("No Pneumonia", "Pneumonia")),
    model = "ranger_offset"
  )


## ========================================================================== ##
## 5. Combine Data Sets and Plot
## ========================================================================== ##

# combine all data sets
ndf_cr_combined <- rbind(ndf_cr_aj
                         # , ndf_cr_pam # exclude due to bias --> prob bug in cr calculation
                         # , ndf_cr_msm # exchanged pamm with xgboost approach
                         , ndf_cr_rfsc 
                         , ndf_cr_ranger
                         , ndf_cr_ranger_offset) %>%
  mutate(cause = factor(cause, levels = c("Discharge", "Death"))) # to ensure correct order in plot


gg_survCurves <- ggplot(ndf_cr_combined, aes(x = tend, y = cif)) +
  geom_line(aes(color = model, linetype = pneu), linewidth = linewidth) +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, linetype = pneu, fill = model), alpha = 0.3) +
  facet_wrap(~cause) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_color_manual(
    name = "Model:",
    values = model_colors,
    breaks = c("rsfc", "ranger", "ranger_offset", "aj"),
    labels = c("RFSC", "Ranger", "Ranger + Offset", "AJ")
  ) +
  scale_linetype_discrete(
    name   = "Pneumonia:",
    labels = c("Pneumonia" = "yes", "No Pneumonia" = "no")
  ) +
  scale_shape_manual(
    values = c("Pneumonia" = 2, "No Pneumonia" = 0),
    name   = "Pneumonia:",
    labels = c("Pneumonia" = "yes", "No Pneumonia" = "no")
  ) +
  scale_fill_manual(values = model_fills) +
  labs(
    x = "Time (in Days)",
    y = "Cumulative Incidence Function"
  ) +
  coord_cartesian(xlim = c(1, 100)) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2, override.aes = list(shape = c(2, 0))),
    shape = guide_legend(order = 2),
    fill = "none"
  )  +
  theme_cif


gg_survCurves




