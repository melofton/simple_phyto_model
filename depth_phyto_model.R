get_T_parms <- function(group_parms){

  # unpack parms
  T_std = group_parms$T_std
  T_opt = group_parms$T_opt
  T_max = group_parms$T_max

  theta = group_parms$theta_growth

  t20 = 20
  tol   = 0.05
  inn = 1
  a0 = theta^(T_std-t20)
  a1 = theta^(T_opt-t20)
  a2 = theta^(T_max-t20)

  # Perform the iteration to find the constants.
  # First approximation of k.
  k = 6.0
  i = 0
  G = tol + 1.0
  curvef = TRUE
  # Do the iterations until -tol < G < tol

  repeat{

    i=i+1

    if(i == 100){ #increases the tolerance if more than 100 iterations performed
      i=0
      tol=tol+0.01
    }

    if(curvef == TRUE){ # Use the condition f(T)=v**(T-20) at T=Tsta
      G = k * theta^(k * T_opt) * a2 - a1 * (theta^(k * T_max) - theta^(k * T_std))
      devG = theta^(k * T_opt) * a2 * (inn + k * T_opt * log(theta)) - a1 * log(theta) * (T_max * theta^(k * T_max) - T_std * theta^(k * T_std))
    } else { # Use the condition f(T)=1 at T=Tsta
      G = k * theta^(k * T_opt) * (a0 - a2 - inn) - a1 * (theta^(k * T_std) - theta^(k * T_max))
      devG = (a0 - a2 - inn) * theta^(k * T_opt) * (inn + k * T_opt * log(theta)) - a1 * log(theta) * (T_std * theta^(k * T_std) - T_max * theta^(k * T_max))
    }

    # Find the next iteration of k
    k = k - G / devG

    if((G <= -tol) | (G >= tol)){
      break
    }
  }

  # Get the remaining model constants
  if(k != 0.0){
    a=-log(a1/(k*theta^(k*T_opt)))/(k*log(theta))
    if(curvef == TRUE){
      b=theta^(k*(T_std-a))
    } else {
      b=inn+theta^(k*(T_std-a))-a0
    }
  } else {
    a=0.0
    b=0.0
  }

  # Set the model constants to the calculated values
  kTn = k
  aTn = a
  bTn = b

  return(list(kTn = kTn, aTn = aTn, bTn = bTn))

}

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
  theta_growth <- parms[3]
  light_extinction <- parms[4]
  I_K <- parms[5]
  N_o <- parms[6]
  K_N <- parms[7]
  P_o <- parms[8]
  K_P <- parms[9]
  f_pr <- parms[10]
  R_resp <- parms[11]
  theta_resp <- parms[12]
  T_std <- parms[13]
  T_opt <- parms[14]
  T_max <- parms[15]
  N_C_ratio <- parms[16]
  P_C_ratio <- parms[17]
  phyto_flux_top <- parms[18]
  #area <- parms[19]
  lake_depth <- parms[20]
  num_boxes <- parms[21]
  KePHY <-parms[22]
  D_temp <- parms[23]

  Tparms <- get_T_parms(group_parms = list(T_std = parms[13], T_opt = parms[14], T_max = parms[15], theta_growth = parms[3]))


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
  fT = NULL
  tp = 20
  kTn = Tparms$kTn
  aTn = Tparms$aTn
  bTn = Tparms$bTn
  for(i in 1:length(layer_temp)){
    fT[i] = 1
    if(layer_temp[i] > T_max){
      fT[i] = 0
    } else if(layer_temp[i] < T_std){
      fT[i] = theta_growth^(layer_temp[i]-tp)
    } else {
      fT[i] = theta_growth^(layer_temp[i]-tp) - theta_growth^(kTn*(layer_temp[i] - aTn)) + bTn
    }
  }
  fT[fT < 0] <- 0.0

  #Nitrogen limitation
  fN <- (NIT - N_o) / (NIT - N_o + K_N)

  #Phosphorus limitation
  fP <- (PHS - P_o) / (PHS - P_o + K_P)

  #Light limitation
  fI = (layer_PAR/I_K) / (1 + (layer_PAR/I_K))

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
       c(fN = fN, fP = fP, fI = fI, fResources = fResources, fT = fT, layer_light_extinction = layer_light_extinction, layer_PAR = layer_PAR))  #This returns the vector of diagnostics
}
