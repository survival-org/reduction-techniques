predict_pv<- function (time, complication){
  #compute pseudo-value
  tumor <- tumor%>%
    mutate(pseudo=pseudosurv(time = days, event = status, tmax = time)$pseudo)
  
  ### Data analysis
  rforest <- ranger(formula = pseudo~complications, data = tumor, replace = FALSE)
  return(as.numeric(predict(rforest, data = tumor %>% mutate(complications = "no"))$predictions[1]*(complication == 'no')+
                      predict(rforest, data = tumor %>% mutate(complications = "yes"))$predictions[1]*(complication == 'yes')))
}

predict_pv_RMST <- function (time, complication){
  #compute pseudo-value
  tumor <- tumor%>%
    mutate(pseudo=pseudomean(time = days, event = status, tmax = time))
  
  ### Data analysis
  rforest <- ranger(formula = pseudo~complications, data = tumor, replace = FALSE)
  return(as.numeric(predict(rforest, data = tumor %>% mutate(complications = "no"))$predictions[1]*(complication == 'no')+
                      predict(rforest, data = tumor %>% mutate(complications = "yes"))$predictions[1]*(complication == 'yes')))
}
