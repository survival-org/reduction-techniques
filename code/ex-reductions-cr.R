## Reduction techniques competing risk

# load libraries
library(dplyr)
library(tidyr)
library(mgcv)
library(pammtools)
library(xgboost)
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
  "pam" = "#D55E00",   # vivid reddish-orange
  "dt"  = "#0072B2",   # deep sky blue
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

# create ped data set
ped_cr <- as_ped(sir.adm, Surv(time, status)~ ., combine = TRUE) |>
  mutate(cause = as.factor(cause), pneu = as.factor(pneu))

# estimate pam
pam <- pamm(ped_status ~ s(tend, by = interaction(cause, pneu)) + cause*pneu, data = ped_cr)

# pem xgb 
lrn_xgb_pem = po("encode", method = "treatment") %>>% 
    lrn(
        'regr.xgboost', 
        nrounds = 1000, 
        eta = 0.12, 
        max_depth = 3, 
        base_score = 1, # very important for poisson loss
        objective = "count:poisson", 
        lambda = 0) |> 
    as_learner()

tsk_pneu = TaskRegr$new(
    id = "pneu", 
    backend = select(ped_cr, ped_status, tend, cause, pneu), 
    target = "ped_status")
lrn_xgb_pem$train(tsk_pneu)


# build new data frame including cif for pam
ndf_cr <- ped_cr |>
  make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause))


haz_xgb_pem = lrn_xgb_pem$predict_newdata(ndf_cr) |> as.data.table()
ndf_cr_xgb_pem <- ndf_cr |> cbind(haz_xgb_pem = haz_xgb_pem[["response"]])

ndf_cr_xgb_pem <- ndf_cr_xgb_pem |>
  group_by(pneu) |> 
  get_cif_from_haz(haz_col = "haz_xgb_pem") |>
  ungroup()
