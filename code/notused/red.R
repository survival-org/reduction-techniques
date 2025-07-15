# load packages
library(dplyr)
library(tidyr)
library(mgcv)
library(survival)
library(pammtools)
library(pseudo)
library(geepack)
library(ggplot2)
library(doBy)

# initiagee# initialize variables
#setwd("C:/Users/ra63liw/Documents/98_git/reduction-techniques")                                   # !! set wd individually !!
linewidth = 1
label_size = 24



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


## -------------------------------------------------------------------------- ##
## PV - Léa
## -------------------------------------------------------------------------- ##

### add a patID variable in the dataset 
tumor$patID <- 1:length(tumor$days)

predict_pv <- function (time, complication){
  #compute pseudo-value
  pv <- pseudosurv(tumor$days, tumor$status, tmax = time)
  
  #One pseudo-value per patients
  tumor$pv<-as.vector(pv$pseudo)
  tumor$ipv<-as.vector(1-pv$pseudo) # needed because the cloglog function implemented in geepack is log(log(1-x))
  
  ### Data analysis
  ### Univariate analysis
  fit <- geese(ipv ~ complications, 
               data = tumor, id = patID,  mean.link="cloglog",
               corstr="independence", family = gaussian())
  return(as.numeric(exp(-exp(c(1,1*(complication == 'yes'))%*%fit$beta))))
}
### test
predict_pv(time = 100, complication = 'yes')


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

## PV
pred_pv =  data.frame(tend = rep(cutpoints[2:length(cutpoints)],2), 
                      complications = c(rep('yes', length(cutpoints)-1), rep('no', length(cutpoints)-1)))

pred_pv$prob_pv = mapply(predict_pv, pred_pv$tend, pred_pv$complications) #may take a few seconds

## combine results
pred <- pred_pam %>%
  left_join(pred_dt, by = c("tend", "complications"))%>%
  left_join(pred_pv, by = c("tend", "complications"))

pred_long <- pred %>%
  pivot_longer(
    cols = c(prob_pam, prob_dt, prob_pv),
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
    values = c("pam" = "firebrick2", "dt" = "steelblue",'pv'='springgreen4', "km" = "black"),
    breaks = c("pam", "dt", 'pv', "km"),
    labels = c("PAM", "DT", 'PV', "KM")
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

## -------------------------------------------------------------------------- ##
# RMST
## -------------------------------------------------------------------------- ##

## -------------------------------------------------------------------------- ##
## KM
## -------------------------------------------------------------------------- ##
rmst_km_time <- km_df |>
  group_by(complications) |>
  mutate(rmst = cumsum(prob * diff(c(0, tend)))
         , model = "km") |>
  select(-c(prob, km_low, km_high))

## -------------------------------------------------------------------------- ##
## PEM / PAM
## -------------------------------------------------------------------------- ##
# calculate rmst (integrate step function)
rmst_pam_time <- pred_pam %>%
  group_by(complications) %>%
  mutate(rmst = cumsum(prob_pam * diff(c(0, tend)))
         , model = "pam") |>
  select(-prob_pam)

## -------------------------------------------------------------------------- ##
## DT
## -------------------------------------------------------------------------- ##

# calculate rmst (integrate step function)
rmst_dt_time <- pred_dt %>%
  group_by(complications) %>%
  mutate(rmst = cumsum(prob_dt * diff(c(0, tend)))
         , model = "dt") |>
  select(-prob_dt)

## -------------------------------------------------------------------------- ##
## PV
## -------------------------------------------------------------------------- ##
predict_pv_RMST <- function (time, complication){
  #compute pseudo-value
  pv <- pseudomean(tumor$days, tumor$status, tmax = time)
  #One pseudo-value per patients
  tumor$rmst<-as.vector(pv)
  
  ### Data analysis
  fit <- geese(rmst ~ complications, data =tumor, id = patID, mean.link = "identity")
  return(as.numeric(c(1,1*(complication == 'yes'))%*%fit$beta))
}
### test
predict_pv_RMST(time = 1000, complication = 'yes')
predict_pv_RMST(time = 1000, complication = 'no')

## ---
## difference in rmst can be taken directly from estimates in summary 
time = 1000
complication = 'yes'

pv <- pseudomean(tumor$days, tumor$status, tmax = time)
#One pseudo-value per patients
tumor$rmst<-as.vector(pv)
fit <- geese(rmst ~ complications + age + resection, data =tumor, id = patID, mean.link = "identity")
fit_lm <- lm(rmst ~ complications + age + resection, data =tumor)
summary(fit)
summary(fit_lm)
## ---


# TBD
# table must contain
# - time
# - complications
# - rmst
# - model = "pv"

rmst_pv_time = data.frame(tend = rep(cutpoints[2:length(cutpoints)],2), 
                                     complications = c(rep('yes', length(cutpoints)-1), rep('no', length(cutpoints)-1)), 
                          model = rep('pv', 2*length(cutpoints)-2))

rmst_pv_time$rmst = mapply(predict_pv_RMST, rmst_pv_time$tend, rmst_pv_time$complications) #may take a few seconds

## merge results
rmst_time <- rbind(rmst_km_time
                   , rmst_pam_time
                   , rmst_dt_time
                   , rmst_pv_time
)

# plot rmst over time
gg_rmst_time <- ggplot(rmst_time, aes(x = tend, y = rmst)) +
  geom_line(aes(color = model, linetype = complications), linewidth = linewidth) +
  scale_color_manual(
    name = "model",
    values = c("km" = "black", "pam" = "firebrick2", "dt" = "steelblue", "pv" = "springgreen4"),
    breaks = c("km" ,"pam", "dt", "pv"),
    labels = c("KM" ,"PAM", "DT", "PV")
  ) +
  theme_minimal() +
  theme(legend.position = "right")

gg_rmst_time
ggsave("figures/rmst_over_time.png", gg_rmst_time, width = 10, height = 6, dpi = 300)


## -------------------------------------------------------------------------- ##
## Backup: Further plots
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

## PV

rmst_pv <- pred_pv %>%
  group_by(complications) %>%
  summarise(rmst = sum(prob_pv * diff(c(0, tend))))

# transform into wide format for plotting area between survival curves
pred_pv_wide <- pred_pv |> select(tend
                                  , complications
                                  , prob_pv) |>
  pivot_wider(names_from = complications
              , values_from = prob_pv) |>
  mutate(diff = abs(yes - no))

# plot survival curves and area corresponding to rmst
gg_rmst_pv <- ggplot(pred_pv_wide, aes(x = tend)) +
  geom_line(aes(y = yes), col="firebrick2") +
  geom_line(aes(y = no), col="steelblue") +
  # Ribbon for the shaded area
  geom_ribbon(aes(ymin = yes, ymax = no), fill = "grey", alpha = 0.2) +
  labs(y = "Survival Probability", x = "Time",
       title = "RMST complications 'yes' and 'no'",
       subtitle = "Grey area represents RMST difference") +
  theme_minimal() +
  theme(legend.position = "none")

