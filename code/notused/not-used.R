
## -------------------------------------------------------------------------- ##
## Backup: Further plots for tumor illustrations
##  NOT USED
## -------------------------------------------------------------------------- ##

## PAM

# calculate rmst (integrate step function)
rmst_pam <- pred_pam %>%
  group_by(complications) %>%
  summarise(rmst = sum(prob_pam * diff(c(0, tend))))

# transform into wide format for plotting area between survival curves
pred_pam_wide <- pred_pam |> select(tend
                                    , complications
                                    , prob_pam) |>
  pivot_wider(names_from = complications
              , values_from = prob_pam) |>
  mutate(diff = abs(yes - no))

# plot survival curves and area corresponding to rmst
gg_rmst_pam <- ggplot(pred_pam_wide, aes(x = tend)) +
  geom_line(aes(y = yes), col="firebrick2") +
  geom_line(aes(y = no), col="steelblue") +
  # Ribbon for the shaded area
  geom_ribbon(aes(ymin = yes, ymax = no), fill = "grey", alpha = 0.2) +
  labs(y = "Survival Probability", x = "Time",
       title = "RMST complications 'yes' and 'no'",
       subtitle = "Grey area represents RMST difference") +
  theme_minimal() +
  theme(legend.position = "none")


## DT

rmst_dt <- pred_dt %>%
  group_by(complications) %>%
  summarise(rmst = sum(prob_dt * diff(c(0, tend))))

# transform into wide format for plotting area between survival curves
pred_dt_wide <- pred_dt |> select(tend
                                  , complications
                                  , prob_dt) |>
  pivot_wider(names_from = complications
              , values_from = prob_dt) |>
  mutate(diff = abs(yes - no))

# plot survival curves and area corresponding to rmst
gg_rmst_dt <- ggplot(pred_dt_wide, aes(x = tend)) +
  geom_line(aes(y = yes), col="firebrick2") +
  geom_line(aes(y = no), col="steelblue") +
  # Ribbon for the shaded area
  geom_ribbon(aes(ymin = yes, ymax = no), fill = "grey", alpha = 0.2) +
  labs(y = "Survival Probability", x = "Time",
       title = "RMST complications 'yes' and 'no'",
       subtitle = "Grey area represents RMST difference") +
  theme_minimal() +
  theme(legend.position = "none")

## -------------------------------------------------------------------------- ##
# # vertical alignment
## -------------------------------------------------------------------------- ##

gg_survCurves_adj <- gg_survCurves +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

gg_rmst_time_adj <- gg_rmst_time +
  theme(
    legend.position = "bottom",
    # legend.box = "horizontal", # uncomment to have both legend items next to each other
    legend.direction = "horizontal",
    legend.spacing.y = unit(0, "cm") # decreases vertical gap between lines
  )

gg_tumor <- gg_survCurves_adj + gg_rmst_time_adj +
  plot_layout(guides = "collect")