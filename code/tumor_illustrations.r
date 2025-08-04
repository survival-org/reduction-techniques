# ============================================================================
# CODE FOR REDUCTION TECHNIQUES PAPER TUMOR EXAMPLE
# includes survival curves + rmst
#
# Structure:
# 0. INITIALIZE LIBRARIES AND FIGURES
# 1. KAPLAN-MEIER ESTIMATES (STRATIFIED BY COMPLICATIONS)
# 2. DISCRETE TIME PIPELINE WITH RANDOM FOREST
# 3. PIECEWISE EXPONENTIAL MODEL (PEM) PIPELINE WITH XGBOOST
# 4. PSEUDO OBSERVATION WITH RANDOM FOREST (TBD)
# 5. COMBINE DATA AND PLOT SURVIVAL
# 6. CALCULATE RMST (ALL ESTIMATES)
# 7. COMBINE DATA AND PLOT RMST
# 8. COMBINE SURVIVAL PLOT WITH RMST PLOT
# 
# ============================================================================

# ============================================================================
# 0. INITIALIZE LIBRARIES AND FIGURES
# ============================================================================

# load packages
library(dplyr)
library(tidyr)
library(mgcv)
library(survival)
library(pammtools)
library(geepack)
library(pseudo)
library(ggplot2)
library(randomForestSRC)
library(mlr3learners)
library(mlr3proba)
library(mlr3pipelines)
library(mlr3extralearners) # for randomForestSRC learner / pak::pak("mlr-org/mlr3extralearners")
theme_set(theme_bw())

# initialize variables
source(here::here("code/functions/calc-pv-and-rmst.R"))

# initalize plotting parameters
linewidth = 1
pointsize = 2
strokewidth = 1.25
# model_colors <- c("pam" = "firebrick2", "dt" = "steelblue", "pv" = "springgreen4", "aj" = "black")
model_colors <- c(
  "xgb" = "#D55E00",   # vivid reddish-orange
  "rf"  = "#0072B2",   # deep sky blue
  "pv"  = "#009E73",   # bluish green (unchanged)
  "km"  = "#000000"    # black (unchanged)
)

model_fills  <- rep("darkgrey", 4); names(model_fills) <- names(model_colors)

theme_reduction <- theme(
  text = element_text(size = 10),
  axis.title = element_text(size = 11),
  axis.text = element_text(size = 9),
  legend.position = "bottom",
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),
  strip.text = element_text(size = 11, face = "bold")
)

# Shared color scale
shared_color <- scale_color_manual(
  name = "Model:",
  values = model_colors,
  breaks = c("xgb", "rf", "pv", "km"),
  labels = c("XGB PEM", "RF DT", "RF PV", "KM")
)

# Shared linetype
shared_linetype <- scale_linetype_discrete(
  name = "Complications:",
  labels = c("Yes" = "yes", "No" = "no")
)

# Shared shape
shared_shape <- scale_shape_manual(
  values = c("yes" = 2, "no" = 0),
  name   = "Complications:",
  labels = c("Yes" = "yes", "No" = "no")
)

# Shared fill (to disable fill legend)
shared_fill <- scale_fill_manual(values = c("xgb" = "darkgrey", "rf" = "darkgrey", "pv" = "darkgrey", "km" = "darkgrey"))

# Shared guides
shared_guides <- guides(
  color = guide_legend(order = 1),
  linetype = guide_legend(order = 2, override.aes = list(shape = c(2, 0))),
  shape = guide_legend(order = 2),
  fill = "none"
)


# Load data
data("tumor", package = "pammtools")
head(tumor)

# Create survival task
tsk_tumor = TaskSurv$new( 
  id = "tumor", 
  backend = select(tumor, days, status, complications), 
  time = "days", 
  event = "status")

# Define cut points for discrete time models
cut = seq(0, max(tumor$days), length.out = 100) |> round()

# ============================================================================
# 1. KAPLAN-MEIER ESTIMATES (STRATIFIED BY COMPLICATIONS)
# ============================================================================

# Fit Kaplan-Meier for each group
km_complications = survival::survfit(Surv(days, status) ~ complications, data = tumor)
km_tidy = broom::tidy(km_complications) |> 
  mutate(
    model = "km", 
    complications = factor(ifelse(strata=="complications=yes", "yes", "no"), levels = c("yes", "no"))
    ) |>  
  select(tend=time, prob=estimate, km_low=conf.low, km_high=conf.high, complications, model)


# ============================================================================
# 2. DISCRETE TIME PIPELINE WITH RANDOM FOREST
# ============================================================================

# Define discrete time pipeline with random forest
pipeline_rf = ppl(
  "survtoclassif_disctime",
  learner = lrn("classif.rfsrc"),
  cut = cut,
  graph_learner = TRUE)


# Train the pipeline
pipeline_rf$train(tsk_tumor)

# Get predictions
pred_rf = pipeline_rf$predict(tsk_tumor)

# Extract survival curves for each group
rf_no_complications <- data.frame(
  tend = as.numeric(names(pred_rf$data$distr[1, ])),
  prob = as.numeric(pred_rf$data$distr[1, ]),
  km_low = NA,
  km_high = NA,
  complications = "no",
  model = "rf"
)

rf_yes_complications <- data.frame(
  tend = as.numeric(names(pred_rf$data$distr[2, ])),
  prob = as.numeric(pred_rf$data$distr[2, ]),
  km_low = NA,
  km_high = NA,
  complications = "yes",
  model = "rf"
)

rf_tidy <- rbind(rf_no_complications, rf_yes_complications) 


# ============================================================================
# 3. PIECEWISE EXPONENTIAL MODEL (PEM) PIPELINE WITH XGBOOST
# ============================================================================

# Define PEM pipeline with XGBoost
pipeline_pem_xgb = po("encode", method = "treatment") %>>% 
  ppl(
    "survtoregr_pem",
    learner = lrn("regr.xgboost", 
                  nrounds = 1000, 
                  eta = 0.12, 
                  max_depth = 2, 
                  base_score = 1, # very important for poisson loss
                  objective = "count:poisson", 
                  lambda = 0),
    cut = cut,
    graph_learner = TRUE) |> as_learner()

pipeline_pem_xgb$train(tsk_tumor)

# Get predictions
pred_xgb_PEM = pipeline_pem_xgb$predict(tsk_tumor)

# Extract survival curves for each group
xgb_no_complications <- data.frame(
  tend = as.numeric(names(pred_xgb_PEM$data$distr[1, ])),
  prob = as.numeric(pred_xgb_PEM$data$distr[1, ]),
  km_low = NA,
  km_high = NA,
  complications = "no",
  model = "xgb"
)

xgb_yes_complications <- data.frame(
  tend = as.numeric(names(pred_xgb_PEM$data$distr[2, ])),
  prob = as.numeric(pred_xgb_PEM$data$distr[2, ]),
  km_low = NA,
  km_high = NA,
  complications = "yes",
  model = "xgb"
)

xgb_tidy <- rbind(xgb_no_complications, xgb_yes_complications)

# ============================================================================
# 4. PSEUDO OBSERVATION WITH RANDOM FOREST (TBD)
# ============================================================================
tumor$patID <- 1:length(tumor$days)
cutpoints <- seq(0, 3800, by=400)

pv_tidy <- data.frame(
  tend = rep(cutpoints[2:length(cutpoints)],2),
  complications = c(rep('yes', length(cutpoints)-1), rep('no', length(cutpoints)-1)),
  km_low = NA,
  km_high = NA,
  model = "pv"
  )

pv_tidy$prob = mapply(predict_pv, pv_tidy$tend, pv_tidy$complications) #may take a few seconds

# ============================================================================
# 5. COMBINE DATA AND PLOT SURVIVAL
# ============================================================================

surv_tidy <- rbind(
  km_tidy
  , xgb_tidy 
  , rf_tidy
  , pv_tidy
  ) |>
  mutate(complications = factor(complications, levels = c("yes", "no")))

surv_lines = subset(surv_tidy, model != "pv")
surv_points = subset(surv_tidy, model == "pv")

gg_survCurves <- ggplot(surv_lines, aes(x = tend, y = prob)) +
  geom_line(aes(color = model, linetype = complications), linewidth = linewidth) +
  geom_ribbon(aes(ymin = km_low, ymax = km_high, linetype = complications, fill = model), alpha = 0.3, show.legend = FALSE) +
  geom_point(data=surv_points, aes(color=model, shape=complications), size=pointsize, stroke=strokewidth) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  shared_color + shared_linetype + shared_shape + shared_fill + shared_guides +
  labs(x = "Time", y = "Survival Probability") +
  theme_reduction

# plot and save
gg_survCurves
ggsave("figures/tumor_survivalCurves.png", gg_survCurves, width = 85, height = 100, units = "mm", dpi = 300)
  

# ============================================================================
# 6. CALCULATE RMST (ALL ESTIMATES)
# ============================================================================

# KM, XGB, RF --> plotted as lines
rmst_lines <- surv_lines |>
  group_by(model, complications) |>
  mutate(rmst = cumsum(prob * diff(c(0, tend)))) |>
  select(-c(prob, km_low, km_high))

rmst_points = data.frame(tend = rep(cutpoints[2:length(cutpoints)],2), 
                          complications = c(rep('yes', length(cutpoints)-1), rep('no', length(cutpoints)-1)), 
                          model = rep('pv', 2*length(cutpoints)-2))

rmst_points$rmst = mapply(predict_pv_RMST, rmst_points$tend, rmst_points$complications) #may take a few seconds

## ========================================================================== ##
## 7. COMBINE DATA AND PLOT RMST
## ========================================================================== ##

rmst_tidy <- rbind(rmst_lines, rmst_points) |>
  mutate(complications = factor(complications, levels = c("yes", "no")))

# needed for legend to be combined
rmst_lines = subset(rmst_tidy, model != "pv")
rmst_points = subset(rmst_tidy, model == "pv")


# plot rmst over time
gg_rmst_time <- ggplot(rmst_lines, aes(x = tend, y = rmst)) +
  geom_line(aes(color = model, linetype = complications), linewidth = linewidth) +
  geom_point(data=rmst_points, aes(color=model, shape=complications), size=pointsize, stroke=strokewidth) +
  shared_color + shared_linetype + shared_shape + shared_fill + shared_guides + guides(fill = "none") +
  labs(x = "Time", y = "Restricted Mean Survival Time") +
  theme_reduction

gg_rmst_time
ggsave("figures/tumor_rmst.png", gg_rmst_time, width = 85, height = 100, units = "mm", dpi = 300)


## ========================================================================== ##
## 8. COMBINE SURVIVAL PLOT WITH RMST PLOT
## ========================================================================== ##

# horizontal alignment
library(patchwork)
gg_tumor <- gg_survCurves + theme(legend.position = "none") + gg_rmst_time + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

gg_tumor
ggsave("figures/tumor_combined.png", gg_tumor, width = 170, height = 100, units = "mm", dpi = 300)
