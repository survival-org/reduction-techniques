## Reduction techniques competing risk

# load libraries
library(dplyr)
library(mgcv)
library(pammtools)
library(ggplot2)
library(mvna)              # pneunomia data set
theme_set(theme_bw())

# load pneunomia data
data(sir.adm, package = "mvna")
head(sir.adm)

## -------------------------------------------------------------------------- ##
## Aalen-Johansen
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## PAM, ndf_cr can be used for DT as well
## -------------------------------------------------------------------------- ##

# create ped data set
ped_cr <- as_ped(sir.adm, Surv(time, status)~ pneu, combine = TRUE) |>
  mutate(cause = as.factor(cause))

# estimate pam
pam <- pamm(ped_status ~ s(tend, by = cause) + cause*pneu, data = ped_cr)

# build new data frame including cif for pam
ndf_cr <- ped_cr |> 
  make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause)) 

ndf_cr_pam <- ndf_cr |>
  group_by(cause, pneu) |> # important!
  add_cif(pam) |> ungroup() |>
  mutate(
    cause = factor(cause, labels = c("Discharge", "Death")),
    pneu = factor(pneu, labels = c("No Pneumonia", "Pneumonia")),
    model = "pam"
  )

## -------------------------------------------------------------------------- ##
## Discrete Time
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## Pseudo Obs
## -------------------------------------------------------------------------- ##

# visualize effect
# tbd: include dt example with color "firebrick2" to be consistent
ggplot(ndf_cr_pam, aes(x = tend, y = cif, col = model)) + geom_line(aes(linetype = pneu, col = model)) +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, linetype = pneu, fill = model), alpha = .3) +
  facet_wrap(~cause) + labs(y = "CIF", x = "time", color = "Model", fill = "Model") +
  scale_color_manual(values = c("firebrick2", "steelblue", "tan4"))
  
  

