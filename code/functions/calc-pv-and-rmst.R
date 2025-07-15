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

predict_pv_RMST <- function (time, complication){
  #compute pseudo-value
  pv <- pseudomean(tumor$days, tumor$status, tmax = time)
  #One pseudo-value per patients
  tumor$rmst<-as.vector(pv)
  
  ### Data analysis
  fit <- geese(rmst ~ complications, data =tumor, id = patID, mean.link = "identity")
  return(as.numeric(c(1,1*(complication == 'no'))%*%fit$beta))
}
