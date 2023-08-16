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
