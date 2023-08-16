phyto_depth_model <- function(t, state, parms, inputs) {
  
  
  #Calculate the environment (note that this depends on time)
  layer_temp <- inputs[[t]]$temp
  PAR_surface <- inputs[[t]]$par
  inflow_n <-  inputs[[t]]$inflow_n
  inflow_p <-  inputs[[t]]$inflow_p
  outflow <- inputs[[t]]$outflow
  areas_mid <- inputs[[t]]$areas_mid
  areas_interface <- inputs[[t]]$areas_interface
  inflow_area <- inputs[[t]]$inflow_area
  outflow_area<- inputs[[t]]$outflow_area
  
  #Unpack parameters
  w_p <- parms[1]
  R_growth <- parms[2]
  light_extinction <- parms[3]
  ksPAR <- parms[4]
  N_o <- parms[5]
  K_N <- parms[6]
  P_o <- parms[7]
  K_P <- parms[8]
  f_pr <- parms[9]
  R_resp <- parms[10]
  theta_resp <- parms[11]
  Tmin <- parms[12]
  Topt <- parms[13]
  Tmax <- parms[14]
  N_C_ratio <- parms[15]
  P_C_ratio <- parms[16]
  phyto_flux_top <- parms[17]
  #area <- parms[18]
  lake_depth <- parms[19]
  num_boxes <- parms[20]
  KePHY <-parms[21]
  D_temp <- parms[22]
  
  #unpack states
  #note that PHYTOS is a vector where each cell is different depth
  PHYTO <- state[1:num_boxes]
  NIT <- state[(num_boxes+1):(2 * num_boxes)]
  PHS <- state[(2* num_boxes+1):(3 * num_boxes)]
  
  #Calcuate the thickness of each layer
  delx <- lake_depth / num_boxes
  
  # Reaction
  #calculate the depth of the middle of the box
  layer_mid_depths <- seq(delx, lake_depth, by = delx) 
  #calculate the light extinction due to the phytos above
  layer_light_extinction <- light_extinction + cumsum(as.vector(PHYTO))*KePHY
  #calculate the PAR at each layer
  layer_PAR <- PAR_surface * exp(-layer_light_extinction * layer_mid_depths)
  
  # Temperature regulation of photosynthesis 
  fT <-  ((layer_temp - Tmin) / (Topt - Tmin)) *((Tmax - layer_temp) / (Tmax - Topt)) ^((Tmax - Topt) / (Topt - Tmin))
  fT[fT < 0] <- 0.0 
  
  #Nitrogen limitation
  fN <- (NIT - N_o) / (NIT - N_o + K_N)
  
  #Phosphorus limitation
  fP <- (PHS - P_o) / (PHS - P_o + K_P)
  
  #Light limitation
  fI <- (layer_PAR / (layer_PAR + ksPAR)) 
  
  #Combined resource limitation
  fResources <- apply(data.frame(fN = fN, fP = fP, fI = fI), 1, FUN = min)
  #fResources <- apply(data.frame(fN = fN, fP = fP), 1, FUN = min) * fI
  
  #primary productivity
  prim_prod <- R_growth * fT * fResources
  
  #Photoexudation
  exudation <- prim_prod * f_pr
  
  #Temperature regulation of respiration
  fT_respiration <- theta_resp^(layer_temp - 20.0)
  
  #Respiration
  respiration <- PHYTO * R_resp * fT_respiration
  
  #Nitrogen uptake associated with primary productivity
  NIT_uptake <- (prim_prod - exudation) * N_C_ratio
  
  #Phorphorus uptake associated with primary productivity
  PHS_uptake <- (prim_prod - exudation) * P_C_ratio
  
  #nutrient turnover associated with respiration
  NIT_turnover <- respiration * N_C_ratio
  PHS_turnover <- respiration * P_C_ratio
  
  #Stream inflow and outflow
  PHYTO_horizontal <- -PHYTO * outflow * (outflow_area / areas_mid)
  NIT_horizontal <- inflow_n * (inflow_area / areas_mid) - NIT * outflow * (outflow_area / areas_mid)
  PHS_horizontal <-  inflow_p * (inflow_area / areas_mid) - PHS * outflow * (outflow_area / areas_mid)
  
  #Net reaction of the states
  PHYTO_reaction <- prim_prod - respiration - exudation + PHYTO_horizontal
  NIT_reaction <- -NIT_uptake + NIT_turnover + NIT_horizontal
  PHS_reaction <- -PHS_uptake + PHS_turnover + PHS_horizontal
  
  # Advection calculation (assume only PHYTOs advect)
  PHYTO_advection_flux <- c(phyto_flux_top, -w_p * PHYTO) * areas_interface
  PHYTO_advection <- -(1/areas_mid) * (diff(PHYTO_advection_flux) / delx)
  
  w_p_nut <- 0.001
  NIT_advection_flux <- c(0, -w_p_nut * NIT) * areas_interface
  NIT_advection <- -(1/areas_mid) * (diff(NIT_advection_flux) / delx)
  
  PHS_advection_flux <- c(0, -w_p_nut * PHS) * areas_interface
  PHS_advection <- -(1/areas_mid) * (diff(PHS_advection_flux) / delx)
  
  
  #Diffusion (assume proportional to temperature gradient)
  
  temp_diff <- diff(layer_temp)
  temp_diff[which(abs(temp_diff) < 0.01)] <- 0.01
  D <- D_temp * c(0, 1/abs(temp_diff), 0)
  
  #Nitrogen
  gradient_middle_boxes <- diff(NIT) 
  gradient <- c(0, gradient_middle_boxes, 0) / delx
  diffusion_flux <- areas_interface * D * gradient
  NIT_diffusion <- (1/areas_mid) * (diff(diffusion_flux) / delx)
  
  #Phorphorus
  gradient_middle_boxes <- diff(PHS) 
  gradient <- c(0, gradient_middle_boxes, 0) / delx
  diffusion_flux <- areas_interface * D * gradient
  PHS_diffusion <- (1/areas_mid) * (diff(diffusion_flux) / delx)
  
  #Net change for each box
  dPHYTO_dt <- PHYTO_advection + PHYTO_reaction
  dNIT_dt <- NIT_advection + NIT_diffusion + NIT_reaction
  dPHS_dt <- PHS_advection + PHS_diffusion + PHS_reaction
  
  list(c(dPHYTO_dt, dNIT_dt, dPHS_dt),  #This returns the vector of derivatives for each layer box
       c(fN = fN, fP = fP, fI = fI, fResources = fResources, fT = fT, layer_light_extinction = layer_light_extinction))  #This returns the vector of diagnostics
}


#Parameters
#Currency = mmolC/m3
#Time scale = day

parms <- c(
  -0.05, #w_p (negative is down, positive is up)
  30, #R_growth
  0.5, #light_extinction
  50, #ksSW
  0.0, #N_o
  2.0, #K_N
  0.0, #P_o
  0.1, #K_P
  0.1, #f_pr
  0.08, #R_resp
  1.08, #theta_resp
  18, #Tmin
  25, #Topt
  35, #Tmax
  0.02, #N_C_ratio
  0.002, #P_C_ratio
  0, #phyto_flux_top
  1, #area (not used)
  9.5,# lake_depth
  38,# num_boxes
  0.005,#KePHYTO
  0.01) #D_temp

# Initial conditions

yini <- rep(0, 3*parms[20])
delx <- parms[19] / parms[20]
depths <- seq(from = delx / 2, by = delx, length.out = parms[20]) 

# Assign a value for each depth
num_boxes <- parms[20]
yini[1:num_boxes] <- approx(x = c(0, 2, parms[19]), y = c(30,0,0), xout = depths, rule = 2)$y
yini[(num_boxes+1):(num_boxes*2)] <- approx(x = c(0, 4, parms[19]), y = c(0,2, 2), xout = depths, rule = 2)$y
yini[(2*num_boxes+1):(num_boxes*3)] <- approx(x = c(0, parms[19]), y = c(0.05, 0.05), xout = depths, rule = 2)$y

# INPUTS

#a quick temporary light function (not used anymore because we have data)

light_function <- function(x){
  time = x
  0.5*(540+440*sin(2*pi*time/365-1.4))
}

#create a list of the inputs for each day
datetime <- seq(lubridate::as_date("2021-01-01"), lubridate::as_date("2022-12-31"), by = "1 day")
num_boxes <- parms[20]

morphometry_H <- c(497.683, 497.983, 498.283, 498.683, 498.983, 499.283, 499.583, 499.883, 500.183, 500.483, 500.783, 501.083, 501.383, 501.683, 501.983, 502.283, 502.583, 502.883, 503.183, 503.483, 503.783, 504.083, 504.383, 504.683, 505.083, 505.383, 505.683, 505.983, 506.283, 506.583, 506.983)
morphometry_A <- c(10, 61.408883, 494.615572, 1201.23579, 2179.597283, 3239.620513, 4358.358439, 5637.911458, 6929.077352, 8228.697419, 9469.324081, 10811.30792, 12399.67051, 14484.22802, 16834.20941, 19631.05422, 22583.1399, 25790.70893, 28442.99667, 31155.95008, 36269.3312, 42851.13714, 51179.89109, 59666.85885, 68146.39437, 76424.14457, 85430.25429, 95068.47603, 103030.4489, 111302.1604, 119880.9164)

morphometry_depth <- abs(max(morphometry_H) - morphometry_H)

temp_data <- readr::read_csv("https://pasta.lternet.edu/package/data/eml/edi/271/7/71e6b946b751aa1b966ab5653b01077f")

daily_temp_data <- temp_data |> 
  dplyr::select(DateTime, ThermistorTemp_C_surface, ThermistorTemp_C_1, ThermistorTemp_C_2,
                ThermistorTemp_C_3, ThermistorTemp_C_4, ThermistorTemp_C_5, ThermistorTemp_C_6, 
                ThermistorTemp_C_7, ThermistorTemp_C_8, ThermistorTemp_C_9) |> 
  dplyr::rename("0.1" = ThermistorTemp_C_surface,
                "1" = ThermistorTemp_C_1,
                "2" = ThermistorTemp_C_2,
                "3" = ThermistorTemp_C_3,
                "4" = ThermistorTemp_C_4,
                "5" = ThermistorTemp_C_5,
                "6" = ThermistorTemp_C_6,
                "7" = ThermistorTemp_C_7,
                "8" = ThermistorTemp_C_8,
                "9" = ThermistorTemp_C_9) |> 
  tidyr::pivot_longer(-DateTime, names_to = "depth", values_to = "temperature") |> 
  dplyr::mutate(date = lubridate::as_date(DateTime)) |> 
  dplyr::summarize(temperature = mean(temperature, na.rm = TRUE), .by = c(depth,date))

inflow_data <- readr::read_csv("https://s3.flare-forecast.org/targets/fcre_v2/fcre/fcre-targets-inflow.csv")

inflow_data <- inflow_data |> filter(variable %in% c("FLOW","NIT_amm", "NIT_nit", "PHS_frp")) |> 
  pivot_wider(names_from = variable, values_from = observation) |> 
  mutate(n_load = FLOW * (NIT_amm + NIT_nit),
         p_load = FLOW * (PHS_frp)) 


met_data <- readr::read_csv("https://pasta.lternet.edu/package/data/eml/edi/389/7/02d36541de9088f2dd99d79dc3a7a853")

sw_met_data <- met_data |> 
  dplyr::select(DateTime, ShortwaveRadiationUp_Average_W_m2) |> 
  dplyr::mutate(date = lubridate::as_date(DateTime)) |> 
  dplyr::summarize(sw = mean(ShortwaveRadiationUp_Average_W_m2, na.rm = TRUE), .by = c(date))

sw_met_data$sw_no_na <- imputeTS::na_interpolation(sw_met_data$sw, option = "linear")  

inputs <- list()
for(i in 1:length(datetime)){
  
  day_inflow <- inflow_data |> filter(datetime == datetime[i])
  
  inflow_n <- rep(0,num_boxes)
  inflow_n[1] <- day_inflow$n_load
  
  inflow_p <- rep(0,num_boxes)
  inflow_p[1] <- day_inflow$p_load
  
  outflow <- rep(0,num_boxes)
  outflow[1] <- day_inflow$FLOW
  
  inflow_area <- 5
  outflow_area <- 5
  
  temps <- daily_temp_data |> dplyr::filter(date == datetime[i])
  if(length(which(!is.na(temps$temperature))) < 2){
    j = 1
    while(length(which(!is.na(temps$temperature))) < 2){
      temps <- daily_temp_data |> dplyr::filter(date == datetime[i-j])
      j <- j + 2
    }
  }
  
  temp_func <- approxfun(x = temps$depth, y = temps$temperature, rule = 2)
  
  light <- sw_met_data |> dplyr::filter(date == datetime[i]) |> 
    dplyr::pull(sw_no_na)
  
  delx <- parms[19] / parms[20]
  depths_mid <- seq(from = delx / 2, by = delx, length.out = parms[20]) 
  depths_interface <- seq(from = 0, to = parms[19], by = delx) 
  
  # use this if you don't know morphometry (assumes 1 meter area per depth)
  areas_mid <- rep(1, length(depths_mid))
  areas_interface <- rep(1, length(depths_interface))
  
  #use this if you have your morphometry defined 
  areas_mid <- approx(x = morphometry_depth, y = morphometry_A, xout = depths_mid, rule = 2)$y
  areas_interface <- approx(x = morphometry_depth, y = morphometry_A, xout = depths_interface, rule = 2)$y
  
  inputs[[i]] <- list(datetime = datetime[i],
                      par = light, 
                      temp = temp_func(depths),
                      inflow_n = inflow_n,
                      inflow_p = inflow_p,
                      outflow = outflow,
                      areas_mid = areas_mid,
                      areas_interface = areas_interface,
                      inflow_area = inflow_area,
                      outflow_area = outflow_area)
  
}

names(inputs) <- datetime

#Use DeSolve to intergrate 

library(deSolve)

simulation_time <- 2 * 365 #DAYS
dt <- 1
times <- seq(1, simulation_time, by = dt)

output <- ode(y = yini, 
              times = times, 
              func = phyto_depth_model, 
              parms = parms,
              inputs = inputs,
              method = "ode45")

#Note that the first column in time and the other columns are the different depths that the derivative are calculated
#Initial conditions
area <- parms[18]
lake_depth <- parms[19]
num_boxes <- parms[20]
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

PHYTO <- output[,PHYTO_index+1 ]
NIT <- output[, NIT_index + 1]
PHS <- output[, PHS_index + 1]
fN <- output[, fN_index + 1]
fP <- output[, fP_index + 1]
fI <- output[, fI_index + 1]
fResources <- output[, fResources_index + 1] 
fT <- output[, fT_index + 1] 
light_extinction <- output[, light_extinction_index + 1] 

# temporal-spatial plot of the concentrations
par(oma = c(0, 0, 3, 0))   # set margin size (oma) so that the title is included
col <- topo.colors

filled.contour(x = times, 
               y = depths, 
               z = as.matrix(PHYTO), 
               color = col, 
               ylim = c(lake_depth, 0), 
               zlim = c(0,50), #range(c(PHYTO)),
               xlab = "time, days", 
               ylab = "Depth, m", 
               main = "Concentration, mmolC/m3")

filled.contour(x = times, 
               y = depths, 
               z = as.matrix(NIT), 
               color = col, 
               ylim = c(lake_depth, 0), 
               zlim = c(0,3), #range(c(NIT)),
               xlab = "time, days", 
               ylab = "Depth, m", 
               main = "Concentration, mmolN/m3")

filled.contour(x = times, 
               y = depths, 
               z = as.matrix(PHS), 
               color = col, 
               ylim = c(lake_depth, 0), 
               zlim = range(c(PHS)),
               xlab = "time, days", 
               ylab = "Depth, m", 
               main = "Concentration, mmolP/m3")

filled.contour(x = times, 
               y = depths, 
               z = as.matrix(fResources), 
               color = col, 
               ylim = c(lake_depth, 0), 
               zlim = range(c(fResources)),
               xlab = "time, days", 
               ylab = "Depth, m", 
               main = "Resource limitation factor")

mtext(outer = TRUE, side = 3, "Vertical phytoplankton model", cex = 1.5)

library(tidyverse)

df_PHYTO <- PHYTO |> 
  dplyr::mutate(datetime = datetime) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "phyto")

df_chla <- df_PHYTO |> 
  dplyr::mutate(prediction = prediction * 0.4,
                variable = "chla")

names(NIT) <- names(PHYTO)


df_NIT <- NIT |> 
  dplyr::mutate(datetime = datetime) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "nit")

names(PHS) <- names(PHYTO)

df_PHS <- PHS |> 
  dplyr::mutate(datetime = datetime) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "phs")

names(fResources) <- names(PHYTO)

df_fResources <- fResources |> 
  dplyr::mutate(datetime = datetime) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "fResources")

names(fN) <- names(PHYTO)

df_fN <- fN |> 
  dplyr::mutate(datetime = datetime) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "fN")

names(fP) <- names(PHYTO)

df_fP <- fP|> 
  dplyr::mutate(datetime = datetime) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "fP")

names(fI) <- names(PHYTO)

df_fI <- fI|> 
  dplyr::mutate(datetime = datetime) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "fI")

names(fT) <- names(PHYTO)
df_fT <- fT|> 
  dplyr::mutate(datetime = datetime) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "fT")

names(light_extinction) <- names(PHYTO)
df_light_extinction <- light_extinction|> 
  dplyr::mutate(datetime = datetime) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "light_extinction")

df_seechi <- df_light_extinction |> 
  filter(depth == 1) |> 
  mutate(prediction = 1.7 / prediction,
         depth = NA, 
         variable = "seechi")


combined <- bind_rows(df_PHYTO, df_NIT, df_PHS, df_fResources, df_fN, df_fP, df_fI, df_fT, df_light_extinction, df_seechi, df_chla)

combined |> 
  filter(depth == 1.0 | is.na(depth)) |> 
  ggplot(aes(x = datetime, y = prediction)) +
  geom_line() +
  facet_wrap(~variable, scale = "free")

combined |> filter(depth %in% c(1,8),
                   variable == "nit") |>
  ggplot(aes(x = datetime, y = prediction, color = factor(depth))) +
  geom_line()
