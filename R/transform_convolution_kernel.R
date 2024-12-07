#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param kernel_daily
#' @param max_daily_diff
#' @param timeperiod_days
#' @return
#' @author Nick Golding
#' @export
# Create a kernel function operating on an arbitrary timestep.
#

# Given a function kernel_daily: the kernel function on a daily timestep, which
# take a vector of delays (e.g. from infection to test, or reporting) in days
# and returns a vector of kernel densities; an integer max_diff_days: the
# maximum delay considered by kernel_daily (all larger values return 0 density);
# and an integer timeperiod_days: the length of the new modelling timeperiod in
# days, return a new function: the kernel of the equivalent convolution on the
# new timeperiod. This takes a vector of delays in units of modelling
# timeperiods and the output is the integralfo the density over these new
# periods.
transform_convolution_kernel <- function(kernel_daily, max_diff_days, timeperiod_days) {

  # maximum number of time periods to consider
  max_timeperiods <- 1 + ceiling(max_diff_days / timeperiod_days)

  # Integral of the function (average days detectable) if not truncated
  total_integral <- sum(kernel_daily(0:max_diff_days))

  # compute the fractions of q for each timeperiod difference

  table <- expand_grid(
    infection_day = seq_len(max_timeperiods * timeperiod_days),
    test_day = seq_len(timeperiod_days)
  ) %>%
    mutate(
      day_diff = infection_day - test_day
    ) %>%
    filter(
      day_diff >= 0
    ) %>%
    mutate(
      kernel_val_day = kernel_daily(day_diff),
      timeperiod_diff = (infection_day - 1) %/% timeperiod_days
    ) %>%
    group_by(
      test_day,
      timeperiod_diff
    ) %>%
    summarise(
      partial_integral = sum(kernel_val_day),
      .groups = "drop"
    ) %>%
    mutate(
      fraction = partial_integral / total_integral
    ) %>%
    group_by(timeperiod_diff) %>%
    summarise(
      fraction = mean(fraction),
      .groups = "drop"
    ) %>%
    mutate(
      kernel_val = fraction * total_integral
    ) %>%
    select(-fraction)

  # return the discrete kernel on the new timeperiod
  function(time_period_difference) {
    idx <- match(time_period_difference, table$timeperiod_diff)
    result <- table$kernel_val[idx]
    result[is.na(result)] <- NA
    result
  }

}
