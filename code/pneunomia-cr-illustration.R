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
library(pseudo)
library(mlr3)
library(mlr3learners)
library(mlr3pipelines)
library(mlr3proba)
theme_set(theme_bw())

# setwd("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_reductionTechniques/reduction-techniques")
setwd("C:/Users/ra63liw/Documents/98_git/reduction-techniques")
source("code/functions/etm-ci-trafo.R")


## initialize variables for plotting
# fontsize
label_size = 24
headline_size = 24

# initialize lines and dots
linewidth = 1
pointsize = 2
strokewidth = 1.5
# model_colors <- c("pam" = "firebrick2", "dt" = "steelblue", "pv" = "springgreen4", "aj" = "black")
model_colors <- c(
  "xgboost" = "#D55E00",   # vivid reddish-orange
  "randforest"  = "#0072B2",   # deep sky blue
  "pv"  = "#009E73",   # bluish green (unchanged)
  "aj"  = "#000000"    # black (unchanged)
)

model_fills  <- rep("darkgrey", 4); names(model_fills) <- names(model_colors)

# theme
# theme_cif <- theme_minimal(base_size = label_size) +
#   theme(
#     strip.text   = element_text(size = headline_size, face = "bold"),
#     axis.title   = element_text(size = label_size),
#     axis.text    = element_text(size = label_size),
#     legend.text  = element_text(size = label_size),
#     legend.position = "right",
#     panel.spacing.x = unit(2, "lines")
#   )

theme_cif <- theme(
  text = element_text(size = 10),
  axis.title = element_text(size = 11),
  axis.text = element_text(size = 9),
  legend.position = "bottom",
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),
  strip.text = element_text(size = 11, face = "bold")
)

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
    , model = "aj"
  )

## -------------------------------------------------------------------------- ##
## PAM, ndf_cr can be used for DT as well
## -------------------------------------------------------------------------- ##
# 
# # create ped data set
# ped_cr <- as_ped(sir.adm, Surv(time, status)~ ., combine = TRUE) |>
#   mutate(cause = as.factor(cause), pneu = as.factor(pneu))
# 
# # estimate pam
# pam <- pamm(ped_status ~ s(tend, by = interaction(cause, pneu)) + cause*pneu, data = ped_cr)
# 
# # build new data frame including cif for pam
# ndf_cr <- ped_cr |>
#   make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause))
# 
# ndf_cr_pam <- ndf_cr |>
#   group_by(cause, pneu) |> # important!
#   add_cif(pam, ci=F, time_var = "tend") |> ungroup() |>
#   mutate(
#     cif_lower = NA,
#     cif_upper = NA,
#     cause = factor(cause, labels = c("Discharge", "Death")),
#     pneu = factor(pneu, labels = c("No Pneumonia", "Pneumonia")),
#     model = "pam"
#   ) |>
#   select(tend, pneu, cause, cif, cif_lower, cif_upper, model)

## compare with transition probabilities

# create ped data set
ped_cr <- as_ped(sir.adm, Surv(time, status)~ ., combine = TRUE, max_time = 150) |>
  mutate(cause = as.factor(cause), pneu = as.factor(pneu))

# estimate pam -> not used really, only placeholder for add_trans_prob to work.
pam_msm <- pamm(ped_status ~ s(tend, by = interaction(cause, pneu)) + cause*pneu, data = ped_cr)

# pem xgb 
lrn_xgb_pem = po("encode", method = "treatment") %>>% 
  lrn(
    'regr.xgboost', 
    nrounds = 1000, 
    eta = 0.025, 
    max_depth = 3, 
    base_score = 1, # very important for poisson loss
    objective = "count:poisson") |> 
  as_learner()

tsk_pneu = TaskRegr$new(
  id = "pneu", 
  backend = select(ped_cr, ped_status, tend, cause, pneu, offset), 
  target = "ped_status")
tsk_pneu$set_col_roles("offset", roles = "offset")
lrn_xgb_pem$train(tsk_pneu)

# build new data frame including cif for pam
ndf_cr_xgboost <- ped_cr |>
  make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause))

# make hazard prediction
haz_xgb_pem = lrn_xgb_pem$predict_newdata(ndf_cr_xgboost) |> as.data.table()
ndf_cr_xgb_pem <- ndf_cr_xgboost |> cbind(haz_xgb_pem = haz_xgb_pem[["response"]])

# calculate cumu hazard
ndf_cr_xgb_pem <- ndf_cr_xgb_pem |>
  group_by(cause, pneu) |> 
  mutate(cumu_hazard = cumsum(haz_xgb_pem * intlen))

# restructure and calculate cifs
ndf_cr_xgb_pem <- ndf_cr_xgb_pem |>
  mutate(transition = ifelse(cause == 1, "1->2", "1->3")) |>
  group_by(pneu, transition) |> 
  add_trans_prob(pam, ci = FALSE) |> # object not used because ci = FaLSE and cumu_hazard in data set
  ungroup() |>
  rename(cif = trans_prob) |>
  mutate(cif_lower = NA
    , cif_upper = NA
    , cause = factor(cause, labels = c("Discharge", "Death"))
    , pneu = factor(pneu, labels = c("No Pneumonia", "Pneumonia"))
    , model = "xgboost"
  ) |>
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

# define equidistant cut points
eventtimes <- unique(sort(sir.adm[sir.adm$status != 0, "time"]))
cut <- sort(union(seq(from=1, to=150, by=10), eventtimes))
# -> no effect

# create ped data set
ped_cr <- as_ped(sir.adm, Surv(time, status)~ ., combine = TRUE, max_time = 150, cut = cut) |>
  mutate(cause = as.factor(cause), pneu = as.factor(pneu))

dt <- gam(
  formula = ped_status ~ s(tend, by = interaction(cause, pneu)) + cause*pneu,
  data = ped_cr,
  family = binomial(link = "logit"))

ped_cr_rf <- ped_cr |>
  mutate(ped_status = as.factor(ped_status))

tsk_pneu = TaskClassif$new(
  id = "pneu", 
  target = "ped_status",
  backend = select(ped_cr_rf, ped_status, tend, cause, pneu, offset)) # offset makes the fit perfect.

## include RF hazard calculation and prediction
lrn_rf_dt = po("encode", method = "treatment") %>>% 
  lrn("classif.ranger", num.trees=1000L) |> as_learner()
lrn_rf_dt$predict_type <- "prob"

# Train the pipeline
lrn_rf_dt$train(tsk_pneu)

ndf_cr <- ped_cr_rf |>
  make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause))

haz_rf_dt = lrn_rf_dt$predict_newdata(ndf_cr)$prob |> as.data.table()
ndf_cr_dt_interim <- ndf_cr |> cbind(hazard = haz_rf_dt$'1')

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

ndf_cr_rf_dt <- ndf_cr_dt_wide %>%
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
    model = "randforest"
  )


## -------------------------------------------------------------------------- ##
## Pseudo Obs
## -------------------------------------------------------------------------- ##

## PV
ci_pv_fun <- function (time, pneumonia, cause = 'Death'){
  if(cause == 'Discharge'){
    #compute pseudo-value
    cipo <- pseudoci(sir.adm$time, sir.adm$status, tmax = time)
    
    #One pseudo-value per patients
    sir.adm$pv<-as.vector(cipo$pseudo[[1]])
    sir.adm$ipv<-as.vector(1-sir.adm$pv) # needed because the cloglog function implemented in geepack is log(log(1-x))
    
    ### Data analysis
    ### Univariate analysis
    fit <- geese(ipv ~ pneu, 
                 data = sir.adm, id = id,  mean.link="cloglog",
                 corstr="independence", family = gaussian())
    return(as.numeric(exp(-exp(c(1,1*(pneumonia == 'Pneumonia'))%*%fit$beta))))
  }else{
    #compute pseudo-value
    cipo <- pseudoci(sir.adm$time, sir.adm$status, tmax = time)
    
    #One pseudo-value per patients
    sir.adm$pv<-as.vector(cipo$pseudo[[2]])
    sir.adm$ipv<-as.vector(1-sir.adm$pv) # needed because the cloglog function implemented in geepack is log(log(1-x))
    
    ### Data analysis
    ### Univariate analysis
    fit <- geese(ipv ~ pneu, 
                 data = sir.adm, id = id,  mean.link="cloglog",
                 corstr="independence", family = gaussian())
    return(as.numeric(exp(-exp(c(1,1*(pneumonia == 'Pneumonia'))%*%fit$beta))))
  }
}
### test
# ci_pv_fun(time = 50, pneumonia = 'No Pneumonia', cause = 'Discharge')
# ci_pv_fun(time = 50, pneumonia = 'Pneumonia', cause = 'Death')

mapply(ci_pv_fun, c(5, 130), 'yes')

cutpoints = seq(5, 130, 5)
ci_pv_data =  data.frame(tend = rep(cutpoints[2:length(cutpoints)],4), 
                         pneumonia = c(rep('yes', length(cutpoints)-1), rep('no', length(cutpoints)-1),
                                       rep('yes', length(cutpoints)-1), rep('no', length(cutpoints)-1)), 
                         cause = c(rep('Discharge', 2*length(cutpoints) - 2), rep('no', 2*length(cutpoints) - 2)))

ndf_cr_pv = ndf_cr_aj |> filter(tend %in% cutpoints)
ndf_cr_pv$model = "pv"
ndf_cr_pv$cif_lower = NA
ndf_cr_pv$cif_upper = NA
ndf_cr_pv$cif = NA
ndf_cr_pv$cif = mapply(ci_pv_fun, ndf_cr_pv$tend, pneumonia = ndf_cr_pv$pneu, cause = ndf_cr_pv$cause) #may take a few seconds

# visualize effect

# combine all data sets
ndf_cr_combined <- rbind(ndf_cr_aj
                         # , ndf_cr_pam # exclude due to bias --> prob bug in cr calculation
                         # , ndf_cr_msm # exchanged pamm with xgboost approach
                         , ndf_cr_xgb_pem 
                         , ndf_cr_rf_dt
                         , ndf_cr_pv) %>%
  mutate(cause = factor(cause, levels = c("Discharge", "Death"))) # to ensure correct order in plot


# prepare plotting
ndf_lines = subset(ndf_cr_combined, model != "pv")
ndf_points = subset(ndf_cr_combined, model == "pv")



# tbd: include dt example with color "firebrick2" to be consistent
gg_survCurves <- ggplot(ndf_lines, aes(x = tend, y = cif)) +
  geom_line(aes(color = model, linetype = pneu), linewidth = linewidth) +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, linetype = pneu, fill = model), alpha = 0.3) +
  geom_point(data=ndf_points, aes(color=model, shape=pneu), size=pointsize, stroke=strokewidth) +
  facet_wrap(~cause) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_color_manual(
    name = "Model:",
    values = model_colors,
    breaks = c("xgboost", "randforest", "pv", "aj"),
    labels = c("XGBoost PEM", "Random Forest DT", "PV", "AJ")
  ) +
  scale_linetype_discrete(
    name   = "Pneumonia:",
    labels = c("Pneumonia" = "yes", "No Pneumonia" = "no")
  ) +
  scale_shape_manual(
    values = c("Pneumonia" = 2, "No Pneumonia" = 0),
    name   = "Pneumonia:",
    labels = c("Pneumonia" = "yes", "No Pneumonia" = "no")
  ) +
  scale_fill_manual(values = c("xgboost" = "darkgrey", "randforest" = "darkgrey", "pv" = "darkgrey", "aj" = "darkgrey")) +
  labs(
    x = "Time (in Days)",
    y = "Cumulative Incidence Function"
  ) +
  coord_cartesian(xlim = c(1, 100)) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2, override.aes = list(shape = c(2, 0))),
    shape = guide_legend(order = 2),
    fill = "none"
  )  +
  theme_cif


gg_survCurves
#ggsave("figures/sir.adm_survivalCurves.png", gg_survCurves, width = 10, height = 6, dpi = 300) # TBD: add PV (whole curves or only points?) as dark orange/brown
ggsave("figures/sir.adm_survivalCurves.png", gg_survCurves, width = 170, height = 100, units = "mm", dpi = 300) # TBD: add PV (whole curves or only points?) as dark orange/brown




