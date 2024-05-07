build_output_df <- function(output, obs, parms, run_datetimes){
  #Note that the first column in time and the other columns are the different depths that the derivative are calculated
  #Initial conditions
  lake_depth <- parms[20]
  num_boxes <- parms[21]
  xcc <- parms[24]
  delx <- lake_depth / num_boxes
  output <- as.data.frame(output)

  #Rename the columns to match the depths
  depths <- seq(from = delx / 2, by = delx, length.out = num_boxes)  # sequence, 1 m intervals
  names(output) <- c("time", rep(depths, 9))
  #output
  #Grab the columns that correspond to the different depths
  PHYTO_index <- 1:num_boxes
  NIT_index <- (num_boxes+1):(2*num_boxes)
  PHS_index <- (2*num_boxes+1):(3*num_boxes)
  fN_index <- (3*num_boxes+1):(4*num_boxes)
  fP_index <- (4*num_boxes+1):(5*num_boxes)
  fI_index <- (5*num_boxes+1):(6*num_boxes)
  fResources_index <- (6*num_boxes+1):(7*num_boxes)
  fT_index <- (7*num_boxes+1):(8*num_boxes)
  light_extinction_index <- (8*num_boxes+1):(9*num_boxes)
  layer_PAR_index <- (9*num_boxes+1):(10*num_boxes)

  PHYTO <- output[,PHYTO_index+1 ]
  NIT <- output[, NIT_index + 1]
  PHS <- output[, PHS_index + 1]
  fN <- output[, fN_index + 1]
  fP <- output[, fP_index + 1]
  fI <- output[, fI_index + 1]
  fResources <- output[, fResources_index + 1]
  fT <- output[, fT_index + 1]
  light_extinction <- output[, light_extinction_index + 1]
  layer_PAR <- output[, layer_PAR_index + 1]

  df_PHYTO <- PHYTO |>
    dplyr::mutate(datetime = run_datetimes) |>
    tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |>
    dplyr::mutate(depth = as.numeric(depth) - 0.125,
                  variable = "phyto")

  df_chla <- df_PHYTO |>
    dplyr::mutate(prediction = prediction * (1/xcc) * 12,
                  variable = "chla")

  names(NIT) <- names(PHYTO)

  df_NIT <- NIT |>
    dplyr::mutate(datetime = run_datetimes) |>
    tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |>
    dplyr::mutate(depth = as.numeric(depth) - 0.125,
                  variable = "din")

  names(PHS) <- names(PHYTO)

  df_PHS <- PHS |>
    dplyr::mutate(datetime = run_datetimes) |>
    tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |>
    dplyr::mutate(depth = as.numeric(depth) - 0.125,
                  variable = "frp")

  names(fResources) <- names(PHYTO)

  df_fResources <- fResources |>
    dplyr::mutate(datetime = run_datetimes) |>
    tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |>
    dplyr::mutate(depth = as.numeric(depth) - 0.125,
                  variable = "fResources")

  names(fN) <- names(PHYTO)

  df_fN <- fN |>
    dplyr::mutate(datetime = run_datetimes) |>
    tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |>
    dplyr::mutate(depth = as.numeric(depth) - 0.125,
                  variable = "fN")

  names(fP) <- names(PHYTO)

  df_fP <- fP|>
    dplyr::mutate(datetime = run_datetimes) |>
    tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |>
    dplyr::mutate(depth = as.numeric(depth) - 0.125,
                  variable = "fP")

  names(fI) <- names(PHYTO)

  df_fI <- fI|>
    dplyr::mutate(datetime = run_datetimes) |>
    tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |>
    dplyr::mutate(depth = as.numeric(depth) - 0.125,
                  variable = "fI")

  names(fT) <- names(PHYTO)

  df_fT <- fT|>
    dplyr::mutate(datetime = run_datetimes) |>
    tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |>
    dplyr::mutate(depth = as.numeric(depth) - 0.125,
                  variable = "fT")

  names(light_extinction) <- names(PHYTO)

  df_light_extinction <- light_extinction|>
    dplyr::mutate(datetime = run_datetimes) |>
    tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |>
    dplyr::mutate(depth = as.numeric(depth) - 0.125,
                  variable = "light_extinction")

  df_seechi <- df_light_extinction |>
    filter(depth == 1) |>
    mutate(prediction = 1.7 / prediction,
           depth = NA,
           variable = "secchi")

  names(layer_PAR) <- names(PHYTO)

  df_layer_PAR <- layer_PAR |>
    dplyr::mutate(datetime = run_datetimes) |>
    tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |>
    dplyr::mutate(depth = as.numeric(depth) - 0.125,
                  variable = "PAR")

  combined <- bind_rows(df_PHYTO, df_NIT, df_PHS,
                        df_fResources, df_fN, df_fP,
                        df_fI, df_fT, df_light_extinction,
                        df_seechi, df_chla, df_layer_PAR) |>
    dplyr::rename(date = datetime)

  df_temp <- obs |>
    filter(variable == "temperature") |>
    add_column(prediction = NA) |>
    select(date, depth, variable, prediction, observation)

  combined <- left_join(combined, obs, by = c("date", "depth", "variable")) |>
    bind_rows(df_temp) |>
    dplyr::rename(datetime = date) |>
    dplyr::select(datetime, depth, variable, prediction,  observation)

  return(combined)
}
