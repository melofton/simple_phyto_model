library(deSolve)
library(tidyverse)
source("build_output_df.R")
source("depth_phyto_model.R")

run_datetimes <- seq(lubridate::as_date("2018-08-06"), lubridate::as_date("2021-12-31"), by = "1 day")

#Parameters
#Currency = mmolC/m3
#Time scale = day

parms <- c(
  -0.001, #w_p (negative is down, positive is up)
  2, #R_growth
  1.08, #theta_growth
  1, #light_extinction
  10, #I_K
  0.0, #N_o
  2, #K_N
  0.0, #P_o
  0.0001, #K_P
  0.1, #f_pr
  0.16, #R_resp
  1.08, #theta_resp
  10, #T_std
  28, #T_opt
  35, #T_max
  0.02, #N_C_ratio
  0.002, #P_C_ratio
  0, #phyto_flux_top
  1, #area (not used)
  9.5,# lake_depth
  38,# num_boxes
  0.005,#KePHYTO
  0.01, #D_temp
  20) #Xcc

# INPUTS

obs <- readr::read_csv("observations.csv", show_col_types = FALSE)
num_boxes <- parms[21]

#create a list of the inputs for each day

morphometry_H <- c(497.683, 497.983, 498.283, 498.683, 498.983, 499.283, 499.583, 499.883, 500.183, 500.483, 500.783, 501.083, 501.383, 501.683, 501.983, 502.283, 502.583, 502.883, 503.183, 503.483, 503.783, 504.083, 504.383, 504.683, 505.083, 505.383, 505.683, 505.983, 506.283, 506.583, 506.983)
morphometry_A <- c(10, 61.408883, 494.615572, 1201.23579, 2179.597283, 3239.620513, 4358.358439, 5637.911458, 6929.077352, 8228.697419, 9469.324081, 10811.30792, 12399.67051, 14484.22802, 16834.20941, 19631.05422, 22583.1399, 25790.70893, 28442.99667, 31155.95008, 36269.3312, 42851.13714, 51179.89109, 59666.85885, 68146.39437, 76424.14457, 85430.25429, 95068.47603, 103030.4489, 111302.1604, 119880.9164)
morphometry_depth <- abs(max(morphometry_H) - morphometry_H)

inputs <- list()
for(i in 1:length(run_datetimes)){

  day_inflow <- obs |> dplyr::filter(date == run_datetimes[i],
                                     variable %in% c("inflow_rate", "n_load", "p_load")) |>
    tidyr::pivot_wider(names_from = variable, values_from = observation)

  inflow_n <- rep(0,num_boxes)
  inflow_n[1] <- day_inflow$n_load

  inflow_p <- rep(0,num_boxes)
  inflow_p[1] <- day_inflow$p_load

  outflow <- rep(0,num_boxes)
  outflow[1] <- day_inflow$inflow_rate

  inflow_area <- 5
  outflow_area <- 5

  temps <- obs |> dplyr::filter(date == run_datetimes[i],
                                variable %in% c("temperature")) |>
    tidyr::pivot_wider(names_from = variable, values_from = observation)

  if(length(which(!is.na(temps$temperature))) < 2){
    j = 1
    while(length(which(!is.na(temps$temperature))) < 2){
      temps <- obs |> dplyr::filter(date == run_datetimes[i-j],
                                    variable %in% c("temperature")) |>
        tidyr::pivot_wider(names_from = variable, values_from = observation)
      j <- j + 2
    }
  }

  temp_func <- approxfun(x = temps$depth, y = temps$temperature, rule = 2)

  light <- obs |> dplyr::filter(date == run_datetimes[i],
                                variable %in% c("shortwave_radiation")) |>
    dplyr::pull(observation)

  delx <- parms[20] / parms[21]
  depths_mid <- seq(from = delx / 2, by = delx, length.out = parms[21])
  depths_interface <- seq(from = 0, to = parms[20], by = delx)

  # use this if you don't know morphometry (assumes 1 meter area per depth)
  areas_mid <- rep(1, length(depths_mid))
  areas_interface <- rep(1, length(depths_interface))

  #use this if you have your morphometry defined
  areas_mid <- approx(x = morphometry_depth, y = morphometry_A, xout = depths_mid, rule = 2)$y
  areas_interface <- approx(x = morphometry_depth, y = morphometry_A, xout = depths_interface, rule = 2)$y

  inputs[[i]] <- list(datetime = run_datetimes[i],
                      par = light,
                      temp = temp_func(depths_mid),
                      inflow_n = inflow_n,
                      inflow_p = inflow_p,
                      outflow = outflow,
                      areas_mid = areas_mid,
                      areas_interface = areas_interface,
                      inflow_area = inflow_area,
                      outflow_area = outflow_area)

}

names(inputs) <- run_datetimes

# Initial conditions

yini <- rep(0, 3*parms[21])
delx <- parms[20] / parms[21]
depths <- seq(from = delx / 2, by = delx, length.out = parms[21])

initial_day <- obs |> dplyr::filter(date == run_datetimes[1])

init_phyto <- initial_day |>
  dplyr::filter(variable == "chla") |>
  dplyr::mutate(observation = observation / (1 / parms[24] * 12))

# Assign a value for each depth
yini[1:num_boxes] <- approx(x = init_phyto$depth, y = init_phyto$observation, xout = depths, rule = 2)$y

init_din <- initial_day |>
  dplyr::filter(variable == "din")
yini[(num_boxes+1):(num_boxes*2)] <- approx(x = init_din$depth, y = init_din$observation, xout = depths, rule = 2)$y

init_frp <- initial_day |>
  dplyr::filter(variable == "frp")
yini[(2*num_boxes+1):(num_boxes*3)] <- approx(x = init_frp$depth, y = init_frp$observation, xout = depths, rule = 2)$y

simulation_time <- length(run_datetimes) #DAYS
dt <- 1
times <- seq(1, simulation_time, by = dt)

output <- ode(y = yini,
              times = times,
              func = phyto_depth_model,
              parms = parms,
              inputs = inputs,
              method = "ode45")

output_df <- build_output_df(output, obs, parms, run_datetimes)

# Visualize output

# assess model run
output_df |>
  filter(depth == 1.5 | is.na(depth)) |>
  ggplot(aes(x = datetime, y = prediction)) +
  geom_point(aes(y = observation)) +
  geom_line(color = "lightblue3") +
  facet_wrap(~variable, scale = "free") +
  theme_bw()

# temporal-spatial plot of the concentrations
par(oma = c(0, 0, 3, 0))   # set margin size (oma) so that the title is included
col <- topo.colors
lake_depth <- parms[20]

mat <- output_df |>
  filter(variable == "phyto") |>
  select(datetime, depth, prediction) |>
  pivot_wider(names_from = depth, values_from = prediction) |>
  select(-datetime)

filled.contour(x = run_datetimes,
               y = depths,
               z = as.matrix(mat),
               color = col,
               ylim = c(lake_depth, 0),
               zlim = range(c(mat)),
               xlab = "time, days",
               ylab = "Depth, m",
               main = "Concentration, mmolC/m3")

mat <- output_df |>
  filter(variable == "din") |>
  select(datetime, depth, prediction) |>
  pivot_wider(names_from = depth, values_from = prediction) |>
  select(-datetime)

filled.contour(x = run_datetimes,
               y = depths,
               z = as.matrix(mat),
               color = col,
               ylim = c(lake_depth, 0),
               zlim = range(c(mat)),
               xlab = "time, days",
               ylab = "Depth, m",
               main = "Concentration, mmolN/m3")

mat <- output_df |>
  filter(variable == "frp") |>
  select(datetime, depth, prediction) |>
  pivot_wider(names_from = depth, values_from = prediction) |>
  select(-datetime)

filled.contour(x = run_datetimes,
               y = depths,
               z = as.matrix(mat),
               color = col,
               ylim = c(lake_depth, 0),
               zlim = range(c(mat)),
               xlab = "time, days",
               ylab = "Depth, m",
               main = "Concentration, mmolP/m3")

mat <- output_df |>
  filter(variable == "PAR") |>
  select(datetime, depth, prediction) |>
  pivot_wider(names_from = depth, values_from = prediction) |>
  select(-datetime)

filled.contour(x = run_datetimes,
               y = depths,
               z = as.matrix(mat),
               color = col,
               ylim = c(lake_depth, 0),
               zlim = range(c(mat)),
               xlab = "time, days",
               ylab = "Depth, m",
               main = "PAR, (umol/m2/s)")

output_df |>
  filter(variable == "temperature") |>
  ggplot(aes(x = datetime, y = observation)) +
  geom_line() +
  facet_wrap(~depth, scale = "free")

output_df |>
  filter(variable == "frp") |>
  ggplot(aes(x = datetime, y = prediction)) +
  geom_line() +
  geom_point(aes(y = observation)) +
  facet_wrap(~depth, scale = "free")

output_df |>
  filter(variable == "din") |>
  ggplot(aes(x = datetime, y = prediction)) +
  geom_line() +
  geom_point(aes(y = observation)) +
  facet_wrap(~depth, scale = "free")

output_df |> filter(depth %in% c(1,8),
                   variable == "din") |>
  ggplot(aes(x = datetime, y = prediction, color = factor(depth))) +
  geom_line()
