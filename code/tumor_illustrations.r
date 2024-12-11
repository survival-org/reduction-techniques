# load packages
library(dplyr)
library(tidyr)
library(mgcv)
library(survival)
library(pammtools)
library(pseudo)
library(geepack)
library(ggplot2)

# helpers


# load data
data("tumor", package = "pammtools")

# transform data
cutpoints <- seq(0, 3800, by=100)
tumor_partitioned <- tumor %>% as_ped(Surv(days, status) ~ ., cut = cutpoints)

# modeling
## Cox
cox <- coxph(
    formula = Surv(days, status) ~ strata(complications),
    data = tumor)

base_df <- basehaz(cox) %>% rename(nelson_aalen = hazard)
ggplot(base_df, aes(x = time, y = nelson_aalen)) +
  geom_stephazard() +
  ylab(expression(hat(Lambda)(t))) + xlab("t") +
  ggtitle("Nelson-Aalen estimate of the cumulative hazard")

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
potsurv <- pseudosurv(tumor$days, tumor$status, tmax = 1:3)
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

# survival curves
cox_df <- broom::tidy(survfit(cox))
ggplot(cox_df, aes(x=time,y=estimate,col=strata)) + geom_step()


new_data <- expand.grid(tend = cutpoints[-1], complications = unique(tumor$complications)) %>%
    mutate(intlen = 100) # hard coded!!!

pred_pam <- new_data %>%
    group_by(complications) %>%
    add_surv_prob(pam) %>%
    ungroup() %>%
    rename(prob_pam = surv_prob, lower_pam = surv_lower, upper_pam = surv_upper) %>%
    select(-c(intlen, lower_pam, upper_pam))

pred_dt <- new_data %>%
    mutate(
        eta_dt = predict(dt, newdata = ., type = "link"),
        hazard_dt = exp(eta_dt) / (1 + exp(eta_dt))
    ) %>%
    group_by(complications) %>%
    mutate(prob_dt = cumprod(1 - hazard_dt)) %>%
    ungroup() %>%
    select(-c(intlen, eta_dt, hazard_dt))

pred <- pred_dt %>%
    left_join(pred_pam, by = c("tend", "complications"))

pred_long <- pred %>%
  pivot_longer(
    cols = c(prob_pam, prob_dt),
    names_to = "model",
    values_to = "prob",
    names_pattern = "prob_(.*)"
  )

gg_baseline <- ggplot(pred_long, aes(x = tend, y = prob, color = interaction(complications, model))) +
  geom_line(linewidth = 1) +
  geom_stephazard(data = base_df, aes(x=time, y = nelson_aalen, col = "Nelson-Aalen")) # tbd
  labs(
    title = "Survival Probabilities by Model and Complications",
    x = "Time (tend)",
    y = "Survival Probability",
    color = "Model & Complications"
  ) +
  theme_minimal()


# RMST
# PEM / PAM
  
  # calculate rmst
  rmst_df <- pred_pam %>%
    group_by(complications) %>%
    summarise(rmst = sum(prob_pam * diff(c(0, tend))))
  
  rmst_df
  
  pred_pam_wide <- pred_pam |> select(tend
                                  , complications
                                  , prob_pam) |> 
    pivot_wider(names_from = complications
                , values_from = prob_pam) |>
    mutate(diff = abs(yes - no))
  
  ggplot(pred_pam_wide, aes(x = tend)) +
    geom_line(aes(y = yes), col="firebrick2") +
    geom_line(aes(y = no), col="steelblue") +
    # Ribbon for the shaded area
    geom_ribbon(aes(ymin = yes, ymax = no), fill = "grey", alpha = 0.2) +
    labs(y = "Survival Probability", x = "Time",
         title = "RMST complications 'yes' and 'no'",
         subtitle = "Grey area represents RMST difference") +
    theme_minimal() +
    theme(legend.position = "none")
  
TBD