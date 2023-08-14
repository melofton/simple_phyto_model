phyto_depth_model <- function(t, state, parms, inputs) {
  
  
  #Calculate the environment (note that this depends on time)
  layer_temp <- inputs[[t]]$temp
  PAR_surface <- inputs[[t]]$par
  inflow_n <-  inputs[[t]]$inflow_n
  inflow_p <-  inputs[[t]]$inflow_p
  outflow <- inputs[[t]]$outflow

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
  area <- parms[18]
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
  PHYTO_horizontal <- -PHYTO * outflow
  NIT_horizontal <- n_inflow - NIT * outflow
  PHS_horizontal <-  inflow_p - PHS * outflow 
  
  #Net reaction of the states
  PHYTO_reaction <- prim_prod - respiration - exudation + PHYTO_horizontal
  NIT_reaction <- -NIT_uptake + NIT_turnover + NIT_horizontal
  PHS_reaction <- -PHS_uptake + PHS_turnover + PHS_horizontal
  
  # Advection calculation (assume only PHYTOs advect)
  PHYTO_advection_flux <- c(phyto_flux_top, w_p * PHYTO) * area
  PHYTO_advection <- -(1/area) * (diff(PHYTO_advection_flux) / delx)
  
  NIT_advection <- 0.0
  PHS_advection <- 0.0
  
  #Diffusion (assume proportional to temperature gradient)
  
  temp_diff <- diff(layer_temp)
  D <- D_temp * abs(temp_diff)
  
  #Nitrogen
  gradient_middle_boxes <- diff(NIT) 
  gradient <- c(0, gradient_middle_boxes, 0) / delx
  diffusion_flux <- area * D * gradient
  NIT_diffusion <- (1/area) * (diff(diffusion_flux) / delx)

  #Phorphorus
  gradient_middle_boxes <- diff(PHS) 
  gradient <- c(0, gradient_middle_boxes, 0) / delx
  diffusion_flux <- area * D * gradient
  PHS_diffusion <- (1/area) * (diff(diffusion_flux) / delx)
  
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
  0.05, #w_p
  10, #R_growth
  0.5, #light_extinction
  120, #ksPAR
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
  1, #area
  12,# lake_depth
  48,# num_boxes
  0.005,#KePHYTO
  0.001) #D_temp
  
# Initial conditions

yini <- rep(0, 3*parms[20])
delx <- parms[19] / parms[20]
depths <- seq(from = delx / 2, by = delx, length.out = parms[20]) 

# Assign a value for each depth
yini[1:numboxes] <- approx(x = c(0, 2, parms[19]), y = c(30,0,0), xout = depths, rule = 2)$y
yini[(numboxes+1):(numboxes*2)] <- approx(x = c(0, parms[19]), y = c(10, 10), xout = depths, rule = 2)$y
yini[(2*numboxes+1):(numboxes*3)] <- approx(x = c(0, parms[19]), y = c(0.1, 0.1), xout = depths, rule = 2)$y

# INPUTS

#a quick temperary light function

light_function <- function(x){
  time = x
  0.5*(540+440*sin(2*pi*time/365-1.4))
}

#create a list of the inputs for each day
datetime <- seq(lubridate::as_date("2021-01-01"), lubridate::as_date("2022-12-31"), by = "1 day")
num_boxes <- parms[20]

inputs <- list()
for(i in 1:length(datetime)){
  
  inflow_n <- rep(0,num_boxes)
  inflow_n[3] <- 0.1
  
  inflow_p <- rep(0,num_boxes)
  inflow_p[3] <- 0.01
  
  outflow <- rep(0,num_boxes)
  outflow[3] <- 0.001 # 0.01
  
  if(month(datetime[i]) %in% c(10,11,12,1,2,3)){
    temp_func <- approxfun(x = c(1, 5, 30), y = c(18, 18, 18), rule = 2)
  }else{
    temp_func <- approxfun(x = c(1, 5, 30), y = c(25, 18, 10), rule = 2)
  }
  
  delx <- parms[19] / parms[20]
  depths <- seq(from = delx / 2, by = delx, length.out = parms[20]) 
  
  inputs[[i]] <- list(datetime = datetime[i],
                      par = light_function(i), 
                      temp = temp_func(depths),
                      inflow_n = inflow_n,
                      inflow_p = inflow_p,
                      outflow = outflow)
  
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
              inputs = inputs)

#Note that the first column in time and the other columns are the different depths that the derivative are calculated
#Initial conditions
area <- parms[18]
lake_depth <- parms[19]
num_boxes <- parms[20]
delx <- lake_depth / numboxes
output <- as.data.frame(output)

#Rename the columns to match the depths
depths <- seq(from = delx / 2, by = delx, length.out = numboxes)  # sequence, 1 m intervals
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
               zlim = range(c(PHYTO)),
               xlab = "time, days", 
               ylab = "Depth, m", 
               main = "Concentration, mmolC/m3")

filled.contour(x = times, 
               y = depths, 
               z = as.matrix(NIT), 
               color = col, 
               ylim = c(lake_depth, 0), 
               zlim = range(c(NIT)),
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
  dplyr::mutate(datetime = 1:nrow(PHYTO)) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "phyto")

names(NIT) <- names(PHYTO)


df_NIT <- NIT |> 
  dplyr::mutate(datetime = 1:nrow(NIT)) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "nit")

names(PHS) <- names(PHYTO)

df_PHS <- PHS |> 
  dplyr::mutate(datetime = 1:nrow(PHS)) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "phs")

names(fResources) <- names(PHYTO)

df_fResources <- fResources |> 
  dplyr::mutate(datetime = 1:nrow(fResources)) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "fResources")

names(fN) <- names(PHYTO)

df_fN <- fN |> 
  dplyr::mutate(datetime = 1:nrow(fN)) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "fN")

names(fP) <- names(PHYTO)

df_fP <- fP|> 
  dplyr::mutate(datetime = 1:nrow(fP)) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "fP")

names(fI) <- names(PHYTO)

df_fI <- fI|> 
  dplyr::mutate(datetime = 1:nrow(fI)) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "fI")

names(fT) <- names(PHYTO)
df_fT <- fT|> 
  dplyr::mutate(datetime = 1:nrow(fT)) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "fT")

names(light_extinction) <- names(PHYTO)
df_light_extinction <- light_extinction|> 
  dplyr::mutate(datetime = 1:nrow(light_extinction)) |> 
  tidyr::pivot_longer(cols = -datetime, names_to = "depth", values_to = "prediction") |> 
  dplyr::mutate(depth = as.numeric(depth) - 0.125,
                variable = "light_extinction")

df_seechi <- df_light_extinction |> 
  filter(depth == 1) |> 
  mutate(prediction = 1.7 / prediction,
         depth = NA, 
         variable = "seechi")
  

combined <- bind_rows(df_PHYTO, df_NIT, df_PHS, df_fResources, df_fN, df_fP, df_fI, df_fT, df_light_extinction, df_seechi)

combined |> 
  filter(depth == 1.0 | is.na(depth)) |> 
  ggplot(aes(x = datetime, y = prediction)) +
  geom_line() +
  facet_wrap(~variable, scale = "free")
