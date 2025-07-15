add_cif_from_haz <- function(
  newdata,
  haz_col = "hazard",
  cause_var = "cause",
  time_var  = "tend") {

  browser()
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
        
        list(
          haz = .df[[haz_col]],
          intlen = .df[["intlen"]]
        )
      }
    )
  
  cumu_hazards <- lapply(hazards, function(x) x$haz * x$intlen)
  overall_survival <- exp(-cumsum(Reduce("+", cumu_hazards)))
  # calculate cif
  # Value of survival just prior to time-point
  
  newdata[["survival"]] <- rep(overall_survival, times = length(unique(causes_model)))
  newdata <- newdata |> dplyr::arrange(.data[[time_var]], .data[[cause_var]], .by_group = TRUE)
  
  newdata[["cif"]] <- cumsum(newdata[[haz_col]] * newdata[["survival"]] * newdata[["intlen"]])

  newdata

}
