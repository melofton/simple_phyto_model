phyto_depth_model <- function(t, state, parms, inputs) {
  
  
  #Calculate the environment (note that this depends on time)
  #This the par at the top of the water column
  #NOTE: this will be an input after we get it working
  temp_func <- approxfun(x = inputs[[t]]$depth, y = inputs[[t]]$temp, rule = 2)
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
  
  #Boundry conditions
  phyto_flux_top <- parms[17]
  
  #Physical description
  area <- parms[18]
  lake_depth <- parms[19]
  num_boxes <- parms[20]
  
  KePHY <-parms[21]
  
  #unpack states
  #note that PHYTOS is a vector where each cell is different depth
  PHYTO <- state[1:num_boxes]
  NIT <- state[(num_boxes+1):(2 * num_boxes)]
  PHS <- state[(2* num_boxes+1):(3 * num_boxes)]
  
  #Calcuate the thickness of each layer
  delx <- lake_depth / num_boxes
  
  # Reaction
  #calculate the depth of the middle of the box
  layer_light_extinction <- light_extinction + cumsum(as.vector(PHYTO))*KePHY
  layer_mid_depths <- seq(delx, lake_depth, by = delx) 
  layer_PAR <- PAR_surface * exp(-layer_light_extinction * layer_mid_depths)
  
  layer_temp <- temp_func(layer_mid_depths)
  
  
  
  # Temperature regulation of photosynthesis 
  fT <-  ((layer_temp - Tmin) / (Topt - Tmin)) *((Tmax - layer_temp) / (Tmax - Topt)) ^((Tmax - Topt) / (Topt - Tmin))
  fT[fT < 0] <- 0.0
  
  fN <- (NIT - N_o) / (NIT - N_o + K_N)
  
  fP <- (PHS- P_o) / (PHS- P_o + K_P)
  
  fI <- (layer_PAR / (layer_PAR + ksPAR)) 
  
  fReources <- apply(data.frame(fN = fN, fP = fP, fI = fI), 1, FUN = min)
  
  primprod <- R_growth * fT * fReources
  
  fT_respiration <- theta_resp^(layer_temp-20.0)
  
  respiration <- PHYTO * R_resp * fT_respiration
  
  exudation <- primprod*f_pr
  
  PHYTO_reaction <- primprod - respiration - exudation

  NIT_uptake <- (primprod - exudation) * N_C_ratio
  
  PHS_uptake <- (primprod - exudation) * P_C_ratio
  
  NIT_turnover <- respiration * N_C_ratio
  PHS_turnover <- respiration * P_C_ratio
  
  PHYTO_horizontal <- -PHYTO * outflow
  NIT_horizontal <- n_inflow - NIT * outflow
  PHS_horizontal <-  inflow_p - PHS * outflow 
  
  PHYTO_reaction <- primprod - respiration - exudation + PHYTO_horizontal
  NIT_reaction <- -NIT_uptake + NIT_turnover + inflow_n + NIT_horizontal
  PHS_reaction <- -PHS_uptake + PHS_turnover + inflow_p + PHS_horizontal
  
  # Advection calculation
  PHYTO_advection_flux <- c(phyto_flux_top, w_p * PHYTO) * area
  PHYTO_advection <- -(1/area) * (diff(PHYTO_advection_flux) / delx)
  
  PHYTO_advection_flux <- c(phyto_flux_top, w_p * PHYTO) * area
  PHYTO_advection <- -(1/area) * (diff(PHYTO_advection_flux) / delx)
  
  NIT_advection <- 0.0
  
  PHS_advection <- 0.0
  
  #Net change for each box
  dPHYTO_dt <- PHYTO_advection + PHYTO_reaction
  dNIT_dt <- NIT_advection + NIT_reaction
  dPHS_dt <- PHS_advection + PHS_reaction
  
  list(c(dPHYTO_dt, dNIT_dt, dPHS_dt), c(fN = fN, fP = fP, fI = fI, fReources = fReources, fT = fT)) #This returns the vector of derivatives for each layer box
  
}


#Parameters
#Currency = mmolC/m3
#Time scale = day

parms <- c(
  0.01, #w_p
  10, #R_growth
  0.1, #light_extinction
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
  0.0002) #KePHYTO
  
numboxes <- parms[20]

#First set all layers equal to zero
yini <- rep(0, numboxes)

#Then initialize the layers where the Phytos are starting
yini[1] <- 3
yini[(numboxes+1):(numboxes*2)] <- 10
yini[(2*numboxes+1):(numboxes*3)] <- 0.1

#Write the csv file with the PAR and temperature data
light_function <- function(x){
  time = x
  0.5*(540+440*sin(2*pi*time/365-1.4))
}
#Create the PAR data
inputs <- list()
for(i in 1:730){
  
  inflow_n <- rep(0,num_boxes)
  inflow_n[1] <- 0.01
  
  inflow_p <- rep(0,num_boxes)
  inflow_p[1] <- 0.001
  
  outflow <- rep(0,num_boxes)
  outflow[1] <- 0.01
  
  inputs[[i]] <- list(par = light_function(i), 
                      temp = c(20,18, 10), 
                      depth = c(1, 5, 30),
                      inflow_n = inflow_n,
                      inflow_p = inflow_p,
                      outflow = outflow
                      )
}

library(deSolve)

#Use DeSolve to intergrate 
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
names(output) <- c("time", depths)
#output

#Grab the columns that correspond to the different depths
PHYTO_index <- 1:num_boxes
NIT_index <- (num_boxes+1):(2*num_boxes)
PHS_index <- (2*num_boxes+1):(3*num_boxes)

PHYTO <- output[,PHYTO_index+1 ]
NIT <- output[, NIT_index + 1]
PHS <- output[, PHS_index + 1]

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

combined <- bind_rows(df_PHYTO, df_NIT, df_PHS)

combined |> 
  filter(depth == 1.5) |> 
  ggplot(aes(x = datetime, y = prediction)) +
  geom_line() +
  facet_wrap(~variable, scale = "free")
