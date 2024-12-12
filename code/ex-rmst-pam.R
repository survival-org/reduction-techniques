## Reduction techniques zucker method

# load libraries
library(dplyr)
library(tidyr)
library(mgcv)
library(pammtools)
library(ggplot2)
theme_set(theme_bw())

# build ped
ped <- tumor %>% as_ped(Surv(days, status)~ complications, id = "id")

# fit gam to ped
pam <- gam(ped_status ~ s(tend, by = complications) + complications, data = ped,
           family = poisson(), offset = offset)

# build new data frame
ped_df <- ped %>%
  make_newdata(tend = unique(tend)
               , complications = unique(complications)) |>
  group_by(complications) |>
  add_surv_prob(pam)

ggplot(ped_df, aes(x = tend, y = surv_prob, group = complications)) +
  geom_ribbon(aes(ymin = surv_lower, ymax = surv_upper, fill = complications), alpha = 0.3) +
  geom_surv(aes(col = complications)) +
  ylab(expression(hat(S)(t))) + xlab(expression(t)) +
  facet_wrap(~complications, labeller = label_both)

# calculate rmst
rmst_df <- ped_df %>%
  group_by(complications) %>%
  summarise(rmst = sum(surv_prob * diff(c(0, tend))))

rmst_df

ped_df_wide <- ped_df |> select(tstart
                                , tend
                                , intlen
                                , interval
                                , offset
                                , complications
                                , surv_prob) |>
  pivot_wider(names_from = complications
              , values_from = surv_prob) |>
  mutate(diff = abs(yes - no))

ggplot(ped_df_wide, aes(x = tend)) +
  geom_line(aes(y = yes), col="firebrick2") +
  geom_line(aes(y = no), col="steelblue") +
  # Ribbon for the shaded area
  geom_ribbon(aes(ymin = yes, ymax = no), fill = "grey", alpha = 0.2) +
  labs(y = "Survival Probability", x = "Time",
       title = "RMST complications 'yes' and 'no'",
       subtitle = "Grey area represents RMST difference") +
  theme_minimal() +
  theme(legend.position = "none")

