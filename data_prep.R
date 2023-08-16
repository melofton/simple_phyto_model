catwalk_data <- readr::read_csv("https://pasta.lternet.edu/package/data/eml/edi/271/7/71e6b946b751aa1b966ab5653b01077f")

daily_temp_data <- catwalk_data |>
  dplyr::select(DateTime, ThermistorTemp_C_surface, ThermistorTemp_C_1, ThermistorTemp_C_2,
                ThermistorTemp_C_3, ThermistorTemp_C_4, ThermistorTemp_C_5, ThermistorTemp_C_6,
                ThermistorTemp_C_7, ThermistorTemp_C_8, ThermistorTemp_C_9, EXOTemp_C_1) |>
  dplyr::rename("0.1" = ThermistorTemp_C_surface,
                "1" = ThermistorTemp_C_1,
                "1.6" = EXOTemp_C_1,
                "2" = ThermistorTemp_C_2,
                "3" = ThermistorTemp_C_3,
                "4" = ThermistorTemp_C_4,
                "5" = ThermistorTemp_C_5,
                "6" = ThermistorTemp_C_6,
                "7" = ThermistorTemp_C_7,
                "8" = ThermistorTemp_C_8,
                "9" = ThermistorTemp_C_9) |>
  tidyr::pivot_longer(-DateTime, names_to = "depth", values_to = "observation") |>
  dplyr::mutate(date = lubridate::as_date(DateTime)) |>
  dplyr::summarize(observation = mean(observation, na.rm = TRUE), .by = c(depth,date))  |>
  dplyr::mutate(variable = "temperature") |>
  dplyr::mutate(depth = as.numeric(depth)) |>
  dplyr::select(date, depth, variable, observation)


exo_chla_data <- catwalk_data |>
  dplyr::select(DateTime, EXOChla_ugL_1) |>
  dplyr::rename("1.6" = EXOChla_ugL_1) |>
  tidyr::pivot_longer(-DateTime, names_to = "depth", values_to = "chla") |>
  dplyr::mutate(date = lubridate::as_date(DateTime)) |>
  dplyr::summarize(observation = mean(chla, na.rm = TRUE), .by = c(depth,date))  |>
  dplyr::mutate(variable = "chla_exo") |>
  dplyr::mutate(depth = as.numeric(depth)) |>
  dplyr::select(date, depth, variable, observation)

ctd <- readr::read_csv("https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf")

ctd <- ctd |> dplyr::filter(Reservoir == "FCR" & Site == 50) |>
  dplyr::rename(depth = Depth_m) |>
  dplyr::select(DateTime, depth, Chla_ugL) |>
  dplyr::mutate(date = lubridate::as_date(DateTime)) |>
  dplyr::summarise(observation = mean(Chla_ugL, na.rm = TRUE), .by = c(date, depth)) |>
  dplyr::mutate(variable = "chla_ctd") |>
  dplyr::mutate(depth = as.numeric(depth)) |>
  dplyr::select(date, depth, variable, observation)

modeled_depths <- seq(0,20, 0.25)

cuts <- tibble::tibble(cuts = as.integer(factor(seq(0,20, 0.25))),
                       depth = seq(0,20, 0.25))

daily_chla_data <- dplyr::bind_rows(exo_chla_data, ctd)

binned_chla <- daily_chla_data |>
dplyr::mutate(cuts = cut(depth, breaks = modeled_depths, include.lowest = TRUE, right = FALSE, labels = FALSE)) |>
  dplyr::group_by(cuts, variable, date) |>
  dplyr::summarize(observation = mean(observation, na.rm = TRUE), .groups = "drop") |>
  dplyr::left_join(cuts, by = "cuts") |>
  dplyr::select(date, variable, depth, observation)

regression_data <- binned_chla |>
  dplyr::filter(depth == 1.5) |>
  tidyr::pivot_wider(names_from = variable, values_from = observation) |>
  na.omit()

fit <- lm(chla_ctd ~ chla_exo, regression_data)

daily_chla_data <- binned_chla |>
  dplyr::mutate(observation = ifelse(variable == "chla_exo", fit$coefficients[1] + observation * fit$coefficients[2], observation)) |>
  dplyr::summarise(observation = mean(observation, na.rm = TRUE), .by = c(date, depth)) |>
  dplyr::mutate(variable = "chla") |>
  dplyr::select(date, depth, variable, observation)

inflow_data <- readr::read_csv("https://s3.flare-forecast.org/targets/fcre_v2/fcre/fcre-targets-inflow.csv")

inflow_data <- inflow_data |>
  dplyr::filter(variable %in% c("FLOW","NIT_amm", "NIT_nit", "PHS_frp")) |>
  tidyr::pivot_wider(names_from = variable, values_from = observation) |>
  dplyr::mutate(n_load = FLOW * (NIT_amm + NIT_nit),
         p_load = FLOW * (PHS_frp)) |>
  dplyr::rename(date = datetime,
                inflow_rate = FLOW) |>
  dplyr::select(date, inflow_rate, n_load, p_load) |>
  tidyr::pivot_longer(-date, names_to = "variable", values_to = "observation") |>
  dplyr::mutate(depth = NA)  |>
  dplyr::mutate(depth = as.numeric(depth)) |>
  dplyr::select(date, depth, variable, observation)



met_data <- readr::read_csv("https://pasta.lternet.edu/package/data/eml/edi/389/7/02d36541de9088f2dd99d79dc3a7a853")

sw_met_data <- met_data |>
  dplyr::select(DateTime, ShortwaveRadiationUp_Average_W_m2) |>
  dplyr::mutate(date = lubridate::as_date(DateTime)) |>
  dplyr::summarize(sw = mean(ShortwaveRadiationUp_Average_W_m2, na.rm = TRUE), .by = c(date))

sw_met_data$sw <- imputeTS::na_interpolation(sw_met_data$sw, option = "linear")

sw_met_data <- sw_met_data |>
  dplyr::rename(observation = sw) |>
  dplyr::mutate(variable = "shortwave_radiation",
                depth = NA)  |>
  dplyr::mutate(depth = as.numeric(depth)) |>
  dplyr::select(date, depth, variable, observation)



nutrients <- readr::read_csv("https://pasta.lternet.edu/package/data/eml/edi/199/11/509f39850b6f95628d10889d66885b76")

nutrients <- nutrients |>
  dplyr::filter(Reservoir == "FCR" & Site == 50) |>
  dplyr::rename(depth = Depth_m) |>
  dplyr::select(DateTime, depth, NO3NO2_ugL, NH4_ugL, SRP_ugL) |>
  dplyr::mutate(date = lubridate::as_date(DateTime),
                NH4 = NH4_ugL * 1000 * 0.001 * (1 / 18.04),
                NO3NO2 = NO3NO2_ugL * 1000 * 0.001 * (1/62.00),
                frp = SRP_ugL * 1000 * 0.001 * (1/95),
                din = NH4 + NO3NO2) |>
  dplyr::select(date, depth, frp, din) |>
  tidyr::pivot_longer(-c(date, depth), names_to = "variable", values_to = "observation") |>
  dplyr::mutate(depth = as.numeric(depth)) |>
  dplyr::select(date, depth, variable, observation) |>
  dplyr::mutate(cuts = cut(depth, breaks = modeled_depths, include.lowest = TRUE, right = FALSE, labels = FALSE)) |>
  dplyr::group_by(cuts, variable, date) |>
  dplyr::summarize(observation = mean(observation, na.rm = TRUE), .groups = "drop") |>
  dplyr::left_join(cuts, by = "cuts") |>
  dplyr::select(date, variable, depth, observation)

secchi <- readr::read_csv("https://pasta.lternet.edu/package/data/eml/edi/198/11/81f396b3e910d3359907b7264e689052")

secchi <- secchi |> dplyr::filter(Reservoir == "FCR" & Site == 50) |>
  dplyr::select(DateTime, Secchi_m) |>
  dplyr::mutate(date = lubridate::as_date(DateTime)) |>
  dplyr::summarise(Secchi_m = mean(Secchi_m, na.rm = TRUE), .by = "date") |>
  dplyr::mutate(variable = "secchi",
                depth = NA,
                depth = as.numeric(depth)) |>
  dplyr::rename(observation = Secchi_m) |>
  dplyr::select(date, depth, variable, observation)



df <- dplyr::bind_rows(daily_temp_data, daily_chla_data, nutrients, inflow_data, sw_met_data, secchi)

readr::write_csv(df, file = "observations.csv")
