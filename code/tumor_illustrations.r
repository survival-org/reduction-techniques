# load packages
library(dplyr)
library(mgcv)
library(survival)
library(pammtools)
library(pseudo)
library(geepack)
library(ggplot2)

# load data
data("tumor", package = "pammtools")

# transform data
tumor_partitioned <- tumor %>% as_ped(Surv(days, status) ~ ., cut = seq(0, 3800, by=100))

# modeling
## Cox
cox <- coxph(
    formula = Surv(days, status) ~ strata(complications),
    data = tumor)

## DT
dt <- gam(
  formula = ped_status ~ s(tend, by = complications),
  data = tumor_partitioned,
  family = binomial(link = "logit"))

## PAM
pam <- pamm(
  formula = ped_status ~ s(tend, by = complications),
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
geese(
    formula = pseudo ~ as.factor(tpseudo), # how to add stratification according to complications?
    id=id,
    data=longpbc3.3,
    mean.link="cloglog",
    corstr="independence",
    family = gaussian(link = "identity"))

# survival curves
TBD

# RMST
TBD