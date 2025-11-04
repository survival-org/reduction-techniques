library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

dir_data <- paste0(here::here(), "/data")

# helpers ----
sim_weibull <- function(
  formula,
  n = 500,
  cut = seq(0, 10, by = 0.1),
  lambda_cens = 0.1
) {

  df <- tibble::tibble(
    id = 1:n,
    time = NA,
    status = NA,
    x1 = rnorm(n, 0, 1),
    x2 = rbinom(n, 1, 0.5),
    x3 = runif(n, -1, 1)
  )

  # logh(t)=log(α)+(α−1)log(t)−αlog(s) (so base R / English wiki parametrization of Weibull distribution)
  # logh(t)=(time terms) + linear_predictor
  # linear_predictor = -alpha*log(s) --> s = exp(-linear_predictor/alpha)
  linear_predictor <- with(df, eval(parse(text = as.character(formula)[2])))
  shape <- ifelse(df$x2 == 1, 1.0, 3.0)
  scale <- exp(-linear_predictor / shape)

  # inverse-transform sampling
  u <- runif(n)
  survival_time <- scale * (-log(u))^(1 / shape)

  # censoring
  censoring_time <- pmin(rexp(n, rate = lambda_cens), max(cut))
  observed_time <- pmin(survival_time, censoring_time)
  status <- as.integer(survival_time <= censoring_time)

  df$time <- observed_time
  df$status <- status

  return(df)
}

calib_sum <- function(df) {
  n_events <- sum(df$status)
  n_cens <- nrow(df) - n_events
  n_cens_1 <- sum(df$time == 10)
  n_cens_3 <- n_cens - n_cens_1
  sum_events <- df %>% filter(status == 1) %>% pull(time) %>% summary()
  sum_cens <- df %>% filter(status == 0) %>% pull(time) %>% summary()
  hist_events <- df %>%
    filter(status == 1) %>%
    ggplot(aes(x = time)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black") +
    labs(title = paste("Histogram of event times")) +
    xlab("Event time")
  hist_cens <- df %>%
    filter(status == 0) %>%
    ggplot(aes(x = time)) +
    geom_histogram(bins = 50, fill = "lightgreen", color = "black") +
    labs(title = paste("Histogram of censoring times")) +
    xlab("Censoring time")

  out <- list(
    n_events = n_events,
    n_cens = n_cens,
    n_cens_1 = n_cens_1,
    n_cens_3 = n_cens_3,
    sum_events = sum_events,
    sum_cens = sum_cens,
    hist_events = hist_events,
    hist_cens = hist_cens
  )

  return(out)
}

# calibration ----
b0 <- -3.1
b1 <- 1.0
b2 <- 0.9
b3 <- -0.3
b4 <- 0
formula <- as.formula("~ b0 + b1*x1 + b2*x2 + b3*x1*x2 + b4*x3")
cut <- seq(0, 10, by = 0.1)
lambda_cens <- 0.07
n <- 750

instance <- sim_weibull(
  formula = formula,
  n = n,
  cut = cut,
  lambda_cens = lambda_cens
)

calib <- calib_sum(instance)
calib$n_events

n_total = calib$n_events + calib$n_cens
calib$n_cens / n_total
calib$n_cens_1 / n_total
calib$n_cens_3 / n_total
calib$sum_events
calib$sum_cens
calib$hist_events
calib$hist_cens

# parameters ----
b0 <- -3.1
b1 <- 1.0
b2 <- 0.9
b3 <- -0.3
b4 <- 0
formula <- as.formula("~ b0 + b1*x1 + b2*x2 + b3*x1*x2 + b4*x3")
cut <- seq(0, 10, by = 0.1)
lambda_cens <- 0.07
ns <- c(750, 10000)
nsim <- 100

# run ----
datasets <- list()
set.seed(11022022)

start.time <- Sys.time()
for(n in ns) {
  key <- paste0("n", n)
  datasets[[key]] <- list()

  for(i in seq_len(nsim)) {
    instance <- sim_weibull(
      formula = formula,
      n = n,
      cut = cut,
      lambda_cens = lambda_cens
    )
    datasets[[key]][[i]] <- instance
  }
}
end.time <- Sys.time()
run.time <- end.time - start.time
print(run.time)

# save ----
saveRDS(datasets, file = file.path(dir_data, "datasets_sim_weibull.rds"))
