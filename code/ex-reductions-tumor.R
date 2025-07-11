library(mlr3)
library(mlr3learners)
library(mlr3extralearners)
library(mlr3pipelines)
library(mlr3proba)
library(survival)
library(pammtools)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

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
km_tidy = broom::tidy(km_complications)

# ============================================================================
# 2. DISCRETE TIME PIPELINE WITH RANDOM FOREST
# ============================================================================

# Define discrete time pipeline with random forest
pipeline_rf = ppl(
  "survtoclassif_disctime",
  learner = lrn("classif.ranger", num.trees = 1000L),
  cut = cut,
  graph_learner = TRUE)

# Train the pipeline
pipeline_rf$train(tsk_tumor)

# Get predictions
pred_rf = pipeline_rf$predict(tsk_tumor)

# Extract survival curves for each group
rf_no_complications = pred_rf$data$distr[1, ]
rf_yes_complications = pred_rf$data$distr[2, ]

# ============================================================================
# 3. PIECEWISE EXPONENTIAL MODEL (PEM) PIPELINE WITH XGBOOST
# ============================================================================

# Define PEM pipeline with XGBoost
lrn_xgb_regr_depth2 = lrn(
    id = 'regr.xgboost', 
    nrounds = 1000, 
    eta = 0.12, 
    max_depth = 2, 
    base_score = 1, # very important for poisson loss
    objective = "count:poisson", 
    lambda = 0)

po_xgb_pem = po("encode", method = "treatment") %>>% 
    lrn_xgb_regr_depth2 |>as_learner()
ppl_xgb_pem = ppl(
  "survtoregr_pem",
  learner = po_xgb_pem,
  cut = cut,
  graph_learner = TRUE)

ppl_xgb_pem$train(tsk_tumor)
pred_xgb_PEM = ppl_xgb_pem$predict(tsk_tumor, row_ids = c(1,2))

# Get predictions
# Extract survival curves for each group
pem_xgb_no_complications = pred_xgb_PEM$data$distr[1, ]
pem_xgb_yes_complications = pred_xgb_PEM$data$distr[2, ]

# ============================================================================
# CREATE COMPARISON PLOT
# ============================================================================

# Create data frame for plotting
plot_data = data.frame(
  time = as.numeric(names(rf_no_complications)),
  # Random Forest predictions
  rf_no = rf_no_complications,
  rf_yes = rf_yes_complications,
  # PEM XGBoost predictions  
  pem_xgb_no = pem_xgb_no_complications,
  pem_xgb_yes = pem_xgb_yes_complications
)

# Create the comparison plot
p_comparison = ggplot() +
  # Kaplan-Meier curves
  geom_surv(data = km_tidy, 
             aes(x = time, y = estimate, col = "Kaplan-Meier", lty = strata)) +
  # Random Forest curves
  geom_surv(data = plot_data, 
             aes(x = time, y = rf_no, col = "Discrete Time RF", lty = "no")) +
  geom_surv(data = plot_data, 
             aes(x = time, y = rf_yes, col = "Discrete Time RF", lty = "yes")) +
  # PEM XGBoost curves
  geom_surv(data = plot_data, 
             aes(x = time, y = pem_xgb_no, col = "PEM XGBoost", lty = "no")) +
  geom_surv(data = plot_data, 
             aes(x = time, y = pem_xgb_yes, col = "PEM XGBoost", lty = "yes")) +
  ylim(c(0, 1)) +
  labs(y = "Survival Probability", 
       x = "Time (days)",
       title = "Survival Curves by Complications Status",
       subtitle = "Comparison of Kaplan-Meier, Discrete Time RF, and PEM XGBoost") +
  scale_color_discrete(name = "Method") +
  scale_linetype_discrete(name = "Complications") +
  theme(legend.position = "bottom")

# Display the plot
p_comparison














