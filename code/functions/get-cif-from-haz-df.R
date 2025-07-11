add_cif_from_haz <- function(
  newdata,
  haz_col = "hazard",
  cause_var = "cause",
  time_var  = "tend") {

  purrr::map_dfr(
    split(newdata, dplyr::group_indices(newdata)),
    ~get_cif_from_haz(
      newdata = .x, haz_col = haz_col, cause_var = cause_var, time_var = time_var)
  )

}

get_cif_from_haz <- function(
  newdata,
  haz_col,
  time_var = "tend",
  cause_var = "cause") {


  causes_model <- as.factor(levels(newdata[[cause_var]]))
  cause_data   <- unique(newdata[[cause_var]])

  browser()
  hazards <- purrr::map(
    causes_model,
    ~ {
        .df <- newdata |> filter(cause == .x) |>
          dplyr::arrange(.data[[time_var]], .by_group = TRUE)
        .df[[haz_col]]
      }
    )
  overall_survival <- exp(-cumsum(Reduce("+", hazards) * newdata[["intlen"]]))
  names(hazards) <- causes_model
  # calculate cif
  hazard           <- hazards[[cause_data]]
  # Value of survival just prior to time-point
  survival         <- overall_survivals - 1e-20
  hps              <- hazard * survival
  cifs             <- apply(hps, 2, function(z) cumsum(z * newdata[["intlen"]]))
  newdata[["cif"]] <- rowMeans(cifs)

  newdata

}
