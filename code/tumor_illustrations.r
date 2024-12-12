# load packages
library(dplyr)
library(tidyr)
library(mgcv)
library(survival)
library(pammtools)
library(pseudo)
library(geepack)
library(ggplot2)

# load data
data("tumor", package = "pammtools")

# transform data
cutpoints <- seq(0, 3800, by=100)
tumor_partitioned <- tumor %>% as_ped(Surv(days, status) ~ ., cut = cutpoints)

# modeling ----

## KM
km <- survfit(Surv(days, status) ~ complications, data = tumor, se.fit = TRUE)

## Cox
cox <- coxph(
    formula = Surv(days, status) ~ strata(complications),
    data = tumor)

## DT
dt <- gam(
  formula = ped_status ~ s(tend, by = complications) + complications,
  data = tumor_partitioned,
  family = binomial(link = "logit"))

## PAM
pam <- pamm(
  formula = ped_status ~ s(tend, by = complications) + complications,
  data = tumor_partitioned)

## PV
potsurv <- pseudosurv(tumor$days, tumor$status, tmax = seq(1000, 3800, by=1000)) # can't we do the same seq(0, 3800, by=100) here?
longpbc3 <- NULL
for(it in 1:length(potsurv$time)){
    longpbc3 <- rbind(
        longpbc3,
        cbind(
            tumor,
            pseudo = 1-potsurv$pseudo[,it],
            tpseudo = potsurv$time[it],
            id = 1:nrow(tumor)
        )
    )
}
longpbc3.3 <- longpbc3[order(longpbc3$id),]
pv <- geese(
    formula = pseudo ~ as.factor(tpseudo), # how to add stratification according to complications?
    id=id,
    data=longpbc3.3,
    mean.link="cloglog",
    corstr="independence",
    family = gaussian(link = "identity"))

# survival curves ----
new_data <- tumor_partitioned %>%
  make_newdata(tend = unique(tend), complications = unique(complications))

## KM
km_df <- broom::tidy(km) %>%
  mutate(
    model = "km",
    complications = ifelse(strata=="complications=yes", "yes", "no")
  ) %>%
  select(tend=time, prob=estimate, km_low=conf.low, km_high=conf.high, complications, model)

## Cox
cox_df <- basehaz(cox) %>%
  group_by(strata) %>%
  mutate(prob_km = exp(-hazard)) %>%
  ungroup() %>%
  rename(complications = strata)

## PAM
pred_pam <- new_data %>%
    group_by(complications) %>%
    add_surv_prob(pam) %>%
    ungroup() %>%
    rename(prob_pam = surv_prob) %>%
    select(tend, complications, prob_pam)

## DT
pred_dt <- new_data %>%
    mutate(
        eta_dt = predict(dt, newdata = ., type = "link"),
        hazard_dt = exp(eta_dt) / (1 + exp(eta_dt))
    ) %>%
    group_by(complications) %>%
    mutate(prob_dt = cumprod(1 - hazard_dt)) %>%
    ungroup() %>%
    select(tend, complications, prob_dt)

## combine results
pred <- pred_pam %>%
    left_join(pred_dt, by = c("tend", "complications"))

pred_long <- pred %>%
  pivot_longer(
    cols = c(prob_pam, prob_dt),
    names_to = "model",
    values_to = "prob",
    names_pattern = "prob_(.*)"
  )

## plot
gg_survCurves <- ggplot(pred_long, aes(x = tend, y = prob)) +
  geom_line(aes(color = model, linetype = complications), linewidth = linewidth) +
  geom_stephazard(data = km_df, aes(linetype = complications, color = "km"), linewidth = linewidth) +
  geom_ribbon(data = km_df, aes(ymin = km_low, ymax = km_high, fill = complications), alpha = .3, color = NA) +
  scale_color_manual(
    name = "model",
    values = c("pam" = "blue", "dt" = "red", "km" = "black"),
    breaks = c("pam", "dt", "km"),
    labels = c("PAM", "DT", "KM")
  ) +
  scale_fill_manual(values = c("no" = "darkgrey", "yes" = "darkgrey")) +
  labs(
    x = "Time",
    y = "Survival Probability"
  ) +
  theme_minimal(base_size = label_size) +
  theme(
    axis.title = element_text(size = label_size),
    axis.text = element_text(size = label_size),
    legend.text = element_text(size = label_size),
    legend.position = "right"
  )

gg_survCurves
ggsave("figures/survival_curves.png", gg_survCurves, width = 10, height = 6, dpi = 300) # TBD: add PV (whole curves or only points?) as dark orange/brown

# RMST ----

## PEM / PAM

# calculate rmst (integrate step function)
rmst_pam <- pred_pam %>%
  group_by(complications) %>%
  summarise(rmst = sum(prob_pam * diff(c(0, tend))))

rmst_pam

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
# calculate rmst (integrate step function)
rmst_dt <- pred_dt %>%
  group_by(complications) %>%
  summarise(rmst = sum(prob_dt * diff(c(0, tend))))

rmst_dt

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

## merge results
tbd
