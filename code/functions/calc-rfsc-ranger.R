compute_rsfc <- function(data) {
  eventtimes <- unique(sort(data[data$status != 0, "time"]))
  cut <- sort(union(seq(from = 1, to = 150, by = 1), eventtimes))
  
  ped_cr <- as_ped(data, Surv(time, status) ~ ., combine = TRUE, max_time = 150, cut = cut) %>%
    mutate(cause = as.factor(cause), pneu = as.factor(pneu)) %>%
    mutate(ped_status = as.factor(ped_status))
  
  tsk_pneu <- TaskClassif$new(
    id = "rsfc",
    target = "ped_status",
    backend = select(ped_cr, ped_status, tend, cause, pneu)
  )
  
  learner <- po("encode", method = "treatment") %>>%
    lrn("classif.rfsrc", ntree = 1000L) |> as_learner()
  
  learner$predict_type <- "prob"
  learner$train(tsk_pneu)
  
  ndf_cr <- ped_cr %>%
    make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause))
  
  pred_probs <- learner$predict_newdata(ndf_cr)$prob |> as.data.table()
  ndf_cr <- cbind(ndf_cr, hazard = pred_probs$`1`)
  
  df <- ndf_cr %>%
    pivot_wider(names_from = cause, values_from = hazard, names_prefix = "hazard_", values_fill = 0) %>%
    mutate(hazard_allCause = rowSums(across(starts_with("hazard_")))) %>%
    arrange(pneu, tend) %>%
    group_by(pneu) %>%
    mutate(
      surv_allCause = cumprod(1 - hazard_allCause),
      surv_allCause_lag = lag(surv_allCause, default = 1),
      cif_1 = cumsum(surv_allCause_lag * hazard_1),
      cif_2 = cumsum(surv_allCause_lag * hazard_2)
    ) %>%
    ungroup() %>%
    select(tend, pneu, cif_1, cif_2) %>%
    pivot_longer(cols = starts_with("cif_"),
                 names_to = "cause", names_prefix = "cif_", values_to = "cif") %>%
    mutate(
      cause = factor(cause, levels = c("1", "2"), labels = c("Discharge", "Death")),
      pneu = factor(pneu, labels = c("No Pneumonia", "Pneumonia"))
    )
  return(df)
}

compute_ranger <- function(data) {
  eventtimes <- unique(sort(data[data$status != 0, "time"]))
  cut <- sort(union(seq(from = 1, to = 150, by = 1), eventtimes))
  
  ped_cr <- as_ped(data, Surv(time, status) ~ ., combine = TRUE, max_time = 150, cut = cut) %>%
    mutate(cause = as.factor(cause), pneu = as.factor(pneu)) %>%
    mutate(ped_status = as.factor(ped_status))
  
  tsk_pneu <- TaskClassif$new(
    id = "ranger",
    target = "ped_status",
    backend = select(ped_cr, ped_status, tend, cause, pneu)
  )
  
  learner <- po("encode", method = "treatment") %>>%
    lrn("classif.ranger", num.trees=1000L) |> as_learner()
  
  learner$predict_type <- "prob"
  learner$train(tsk_pneu)
  
  ndf_cr <- ped_cr %>%
    make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause))
  
  pred_probs <- learner$predict_newdata(ndf_cr)$prob |> as.data.table()
  ndf_cr <- cbind(ndf_cr, hazard = pred_probs$`1`)
  
  df <- ndf_cr %>%
    pivot_wider(names_from = cause, values_from = hazard, names_prefix = "hazard_", values_fill = 0) %>%
    mutate(hazard_allCause = rowSums(across(starts_with("hazard_")))) %>%
    arrange(pneu, tend) %>%
    group_by(pneu) %>%
    mutate(
      surv_allCause = cumprod(1 - hazard_allCause),
      surv_allCause_lag = lag(surv_allCause, default = 1),
      cif_1 = cumsum(surv_allCause_lag * hazard_1),
      cif_2 = cumsum(surv_allCause_lag * hazard_2)
    ) %>%
    ungroup() %>%
    select(tend, pneu, cif_1, cif_2) %>%
    pivot_longer(cols = starts_with("cif_"),
                 names_to = "cause", names_prefix = "cif_", values_to = "cif") %>%
    mutate(
      cause = factor(cause, levels = c("1", "2"), labels = c("Discharge", "Death")),
      pneu = factor(pneu, labels = c("No Pneumonia", "Pneumonia"))
    )
  return(df)
}

compute_ranger_offset <- function(data) {
  eventtimes <- unique(sort(data[data$status != 0, "time"]))
  cut <- sort(union(seq(from = 1, to = 150, by = 1), eventtimes))
  
  ped_cr <- as_ped(data, Surv(time, status) ~ ., combine = TRUE, max_time = 150, cut = cut) %>%
    mutate(cause = as.factor(cause), pneu = as.factor(pneu)) %>%
    mutate(ped_status = as.factor(ped_status))
  
  # Use offset even if it's 0 — preserved for structure
  tsk_pneu <- TaskClassif$new(
    id = "ranger_offset",
    target = "ped_status",
    backend = select(ped_cr, ped_status, tend, cause, pneu, offset)
  )
  
  learner <- po("encode", method = "treatment") %>>%
    lrn("classif.ranger", num.trees=1000L) |> as_learner()
  
  learner$predict_type <- "prob"
  learner$train(tsk_pneu)
  
  ndf_cr <- ped_cr %>%
    make_newdata(tend = unique(tend), pneu = unique(pneu), cause = unique(cause))
  
  pred_probs <- learner$predict_newdata(ndf_cr)$prob |> as.data.table()
  ndf_cr <- cbind(ndf_cr, hazard = pred_probs$`1`)
  
  df <- ndf_cr %>%
    pivot_wider(names_from = cause, values_from = hazard, names_prefix = "hazard_", values_fill = 0) %>%
    mutate(hazard_allCause = rowSums(across(starts_with("hazard_")))) %>%
    arrange(pneu, tend) %>%
    group_by(pneu) %>%
    mutate(
      surv_allCause = cumprod(1 - hazard_allCause),
      surv_allCause_lag = lag(surv_allCause, default = 1),
      cif_1 = cumsum(surv_allCause_lag * hazard_1),
      cif_2 = cumsum(surv_allCause_lag * hazard_2)
    ) %>%
    ungroup() %>%
    select(tend, pneu, cif_1, cif_2) %>%
    pivot_longer(cols = starts_with("cif_"),
                 names_to = "cause", names_prefix = "cif_", values_to = "cif") %>%
    mutate(
      cause = factor(cause, levels = c("1", "2"), labels = c("Discharge", "Death")),
      pneu = factor(pneu, labels = c("No Pneumonia", "Pneumonia"))
    )
  return(df)
}

