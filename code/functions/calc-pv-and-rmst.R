predict_pv<- function (time, complication){
  #compute pseudo-value

  #FIXME Tumor shouldn't be used in the function
  pseudo_mat = pseudosurv(time = tumor$days, event = tumor$status, tmax = time)
  pseudo_df = data.frame(pseudo = pseudo_mat$pseudo)

  tmp <- cbind(tumor, pseudo_df)

  # FIXME: sometimes, two values are predicted which are very close to each other
  ### Data analysis
  rforest <- rfsrc(formula = pseudo~complications, data = tmp, replace = FALSE)
  tmp <- cbind(tmp, rforest$predicted)
  
  return(as.numeric(rforest$predicted[tmp$complications == "no"][1]*(complication == 'no')+
                      rforest$predicted[tmp$complications == "yes"][1]*(complication == 'yes')))
}

predict_pv_RMST <- function (time, complication){
  #compute pseudo-value
  pseudo_df <- data.frame(pseudo=pseudomean(time = tumor$days, event = tumor$status, tmax = time))

  tmp <- cbind(tumor, pseudo_df)
  
  ### Data analysis
  rforest <- rfsrc(formula = pseudo~complications, data = tmp, replace = FALSE)
  return(as.numeric(rforest$predicted[tmp$complications == "no"][1]*(complication == 'no')+
                      rforest$predicted[tmp$complications == "yes"][1]*(complication == 'yes')))
}
