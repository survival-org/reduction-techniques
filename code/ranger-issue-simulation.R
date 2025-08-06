library(future)
library(furrr)
library(ggplot2)
 
source(here::here("code/functions/calc-rfsc-ranger.R"))

# Set parallel strategy
plan(multisession, workers = parallel::detectCores() - 1)

n_sim <- 500
eval_times <- seq(1, 150, by = 1)

# Get AJ reference CIs once
aj_ref <- ndf_cr_aj %>%
  filter(tend %in% eval_times) %>%
  select(tend, pneu, cause, cif_lower, cif_upper)

# Simulation function for one iteration
simulate_once <- function(i) {
  library(dplyr)
  library(tidyr)
  library(ranger)
  library(mlr3)
  library(mlr3learners)
  library(mlr3pipelines)
  library(mlr3proba)
  library(mlr3extralearners) # for randomForestSRC learner / pak::pak("mlr-org/mlr3extralearners")
  
  set.seed(1000 + i)  # Ensure reproducibility
  boot_data <- sir.adm[sample(nrow(sir.adm), replace = TRUE), ]
  
  # 🔧 Fix: Assign new unique IDs
  boot_data$id <- seq_len(nrow(boot_data))
  
  # Compute all models
  rsfc_df <- compute_rsfc(boot_data)
  ranger_df <- compute_ranger(boot_data)
  ranger_offset_df <- compute_ranger_offset(boot_data)
  
  # Combine and tag models
  bind_rows(
    rsfc_df %>% mutate(model = "rsfc"),
    ranger_df %>% mutate(model = "ranger"),
    ranger_offset_df %>% mutate(model = "ranger_offset")
  ) %>%
    filter(tend %in% eval_times)
}

# Run all simulations in parallel
simulation_results <- future_map(
  1:n_sim
  , simulate_once
  , .progress = TRUE
  , .options = furrr_options(seed = TRUE)
  )

# Post-process: compute coverage per simulation
coverage_results <- lapply(simulation_results, function(sim_df) {
  sim_df %>%
    left_join(aj_ref, by = c("tend", "pneu", "cause")) %>%
    mutate(covered = (cif >= cif_lower) & (cif <= cif_upper)) %>%
    group_by(model, cause, pneu, tend) %>%
    summarise(coverage = mean(covered), .groups = "drop")
})

# Aggregate results across simulations
coverage_summary <- bind_rows(coverage_results) %>%
  group_by(model, cause, pneu, tend) %>%
  filter(!all(is.na(coverage))) %>%
  summarise(mean_coverage = mean(coverage, na.rm = TRUE), .groups = "drop")

# Plot
gg_coverage <- ggplot(coverage_summary, aes(x = tend, y = mean_coverage, color = model)) +
  geom_line(linewidth = 1) +
  facet_grid(cause ~ pneu) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  labs(
    title = "Coverage Probability over Time",
    x = "Time (Days)",
    y = "Coverage Probability",
    color = "Model"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 14)

gg_coverage
