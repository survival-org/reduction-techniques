## Reduction techniques competing risk

# load libraries
library(dplyr)
library(tidyr)
library(mgcv)
library(pammtools)
library(ggplot2)
library(mvna)              # pneunomia data set
library(etm)               # aalen johansen estimator
library(cmprsk)            # needed for competing risks in etm, cf. Beyersmann p. 79
theme_set(theme_bw())

#setwd("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_reductionTechniques/reduction-techniques/code") # set wd here!!!
setwd("C:/Users/ra63liw/Documents/98_git/reduction-techniques")
source("code/functions/etm-ci-trafo.R")

# required: multi state branch for transition probability calculation
setwd("C:/Users/ra63liw/Documents/98_git/pammtools-multi-state/pammtools") # set pammtools multi state branch here!!!
devtools::load_all()

# load pneunomia data
data(sir.adm, package = "mvna")
head(sir.adm)

## -------------------------------------------------------------------------- ##
## Aalen-Johansen, following Beyersmann chapter 4.3
## -------------------------------------------------------------------------- ##
tra <- matrix(FALSE, ncol = 3, nrow = 3)
dimnames(tra) <- list(c("0", "1", "2"), c("0", "1", "2"))
tra[1, 2:3] <- TRUE
tra

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
    , model = "AJ"
  )

## -------------------------------------------------------------------------- ##
## PAM, ndf_cr can be used for DT as well
## -------------------------------------------------------------------------- ##

# create ped data set
ped_cr <- as_ped(sir.adm, Surv(time, status)~ ., combine = TRUE) |>
  mutate(cause = as.factor(cause), pneu = as.factor(pneu))

# estimate pam
pam <- pamm(ped_status ~ s(tend, by = interaction(cause, pneu)) + cause*pneu, data = ped_cr)

# build new data frame including cif for pam
ndf_cr <- ped_cr |>
  make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause))

ndf_cr_pam <- ndf_cr |>
  group_by(cause, pneu) |> # important!
  add_cif(pam, ci=F, time_var = "tend") |> ungroup() |>
  mutate(
    cif_lower = NA,
    cif_upper = NA,
    cause = factor(cause, labels = c("Discharge", "Death")),
    pneu = factor(pneu, labels = c("No Pneumonia", "Pneumonia")),
    model = "pam"
  ) |>
  select(tend, pneu, cause, cif, cif_lower, cif_upper, model)

## compare with transition probabilities
msm.sir.adm <- sir.adm |> mutate(from = 0
                                 , to = ifelse(status != 0, status, 1)
                                 , transition = ifelse(status != 0, paste0("0->", status), "0->1")
                                 , status = ifelse(status==0, 0, 1)
                                 , tstart = 0
                                 , tstop = time)

msm.sir.adm <- msm.sir.adm |> add_counterfactual_transitions()

ped_msm <- as_ped_multistate(
  data       = msm.sir.adm,
  formula    = Surv(tstart, tstop, status)~ .,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

pam_msm <- gam(
  formula = ped_status ~ s(tend, by = interaction(transition, as.factor(pneu))) + transition*as.factor(pneu),
  data = ped_msm,
  family = poisson(link = "log"),
  offset = offset)

summary(pam_msm)

ndf_msm <- ped_msm |>
  make_newdata(tend = unique(tend)
               , transition = unique(transition)
               , pneu = as.factor(unique(pneu))
  ) |>
  group_by(transition, pneu) |>
  add_trans_prob(pam_msm)

ndf_cr_msm <- ndf_msm |>
  rename(cif = trans_prob) |>
  mutate(cif_lower = NA
         , cif_upper = NA
         , cause = ifelse(transition == "0->1", "Discharge", "Death")
         , model = "pam"
         , pneu = ifelse(pneu == "0", "No Pneumonia", "Pneumonia")) |>
  ungroup() |>
  select(tend
         , cif
         , cif_lower
         , cif_upper
         , cause
         , pneu
         , model)

## -------------------------------------------------------------------------- ##
## Discrete Time
## -------------------------------------------------------------------------- ##

dt <- gam(
  formula = ped_status ~ s(tend, by = interaction(cause, pneu)) + cause*pneu,
  data = ped_cr,
  family = binomial(link = "logit"))

ndf_cr_dt_interim <- ndf_cr %>%
    mutate(eta = predict(dt, newdata = ., type = "link"), hazard = exp(eta) / (1 + exp(eta))) # predict cause-specific hazards

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

ndf_cr_dt <- ndf_cr_dt_wide %>%
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
    model = "dt"
  )

## -------------------------------------------------------------------------- ##
## Pseudo Obs
## -------------------------------------------------------------------------- ##

# visualize effect

# combine all data sets
ndf_cr_combined <- rbind(ndf_cr_aj
                         # , ndf_cr_pam # exclude due to bias --> prob bug in cr calculation
                         , ndf_cr_msm
                         , ndf_cr_dt)

# tbd: include dt example with color "firebrick2" to be consistent
ggplot(ndf_cr_combined, aes(x = tend, y = cif, col = model)) + geom_line(aes(linetype = pneu, col = model)) +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, linetype = pneu, fill = model), alpha = .3) +
  facet_wrap(~cause) +
  labs(y = "CIF", x = "time", color = "Model", fill = "Model") +
  scale_color_manual(
    name = "model",
    values = c("pam" = "firebrick2", "dt" = "steelblue", "km" = "black"),
    breaks = c("pam", "dt", "km"),
    labels = c("PAM", "DT", "KM")
  ) +
  scale_fill_manual(
    name = "model",
    values = c("pam" = "firebrick2", "dt" = "steelblue", "km" = "black"),
    breaks = c("pam", "dt", "km"),
    labels = c("PAM", "DT", "KM")) +
  labs(
    x = "Time",
    y = "Cumulative Incidence Function"
  ) +
  theme_minimal(base_size = label_size) +
  theme(
    axis.title = element_text(size = label_size),
    axis.text = element_text(size = label_size),
    legend.text = element_text(size = label_size),
    legend.position = "right"
  ) +
  xlim(c(0,100)) + 
  ylim(c(0,1))



