
cmb <- function(f, d) function(p) f(p, d)

# ADMB model 20.1 bridged to RTMB
# no input parameter or data changes
base <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  
  RTMB::getAll(pars, data)
  
  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M = exp(log_M)            # Natural mortality
  a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
  a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
  q = exp(log_q)            # Survey catchability
  F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
  F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
  F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass
  
  # Spawning adjustments
  spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
  spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning
  
  # Index values and dimensions
  A = nrow(age_error)                         # Number of ages in model
  A1 = length(ages)                           # Number of ages in comps
  T = sum(catch_ind)                          # Number of fishery years
  Ts = sum(srv_ind)                           # Number of survey years
  Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
  Tsa = sum(srv_age_ind)                      # Number of survey age comp years
  Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
  L = length(length_bins)                     # Number of length bins
  g = 0.00001                                 # Small number to avoid division by zero
  
  # Containers for model outputs
  Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat = rep(0, A)                                              # Initial numbers-at-age
  catch_pred = rep(0, T)                                           # Predicted catch
  srv_pred = rep(0, Ts)                                            # Predicted survey index
  srv_var = rep(0,Ts)                                              # Survey variance
  fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
  srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
  fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
  spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
  N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit
  
  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
  nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
  nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability
  
  # selectivity ----
  slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
  slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
  slx_block[,3] = sel_double_logistic(1:A, a50C[2], deltaC[2], adj=0) # Block 3: double logistic selectivity
  slx_block[,2] = (slx_block[,1] + slx_block[,3]) * 0.5          # Block 2: average of blocks 1 and 3
  slx_block[,4] = to_one(sel_double_logistic(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
  slx_block[,3] = to_one(slx_block[,3])                          # Normalize block 3 selectivity - must be done after block 2 to match ADMB
  
  for(t in 1:T) {
    slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
  }
  
  slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)
  
  
  # mortality ----
  # Calculate fishing mortality for each year
  Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
  for(t in 1:T){
    Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
    Zat[,t] = Fat[,t] + M           # Total mortality at age and year
  }
  Sat = exp(-Zat)                     # Survivorship at age and year
  
  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)
  
  # Remaining columns: survivors from previous year
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
  }
  
  # Calculate recruits and biomasses
  recruits = Nat[1,]                          # Recruitment time series
  spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
  tot_bio = colSums(Nat * waa)                # Total biomass by year
  
  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj = Sat[,T]^(spawn_fract)
  spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
  
  ## catch ----
  # Calculate predicted catch at age and year
  Cat = Fat / Zat * Nat * (1-Sat)
  catch_pred = colSums(Cat * waa) # Predicted catch biomass
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  
  ## survey biomass ----
  isrv = 1
  srv_like = 0.0
  
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
      # Survey likelihood (lognormal, using observed and predicted survey biomass)
      srv_like = srv_like + sum((log(srv_obs[isrv]) - log(srv_pred[isrv]))^2 / 
                                  (2 * (srv_sd[isrv] / srv_obs[isrv])^2))
      isrv = isrv + 1
    }
  }
  
  like_srv = srv_like * srv_wt # Weighted survey likelihood
  
  ## fishery age comp ----
  fish_age_lk = 0.0
  offset = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      # Offset for multinomial likelihood
      offset = offset - fish_age_iss[icomp] * 
        sum((fish_age_obs[,icomp] + g) * 
              log(fish_age_obs[,icomp] + g))
      # Multinomial likelihood for age composition
      fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * 
                                        (fish_age_obs[,icomp] + g) * 
                                        log(fish_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_age_lk = fish_age_lk - offset
  like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood
  
  ## survey age comp ----
  srv_age_lk = 0.0
  offset_sa = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
      # Offset for multinomial likelihood
      offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  srv_age_lk = srv_age_lk - offset_sa
  like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood
  
  ## fishery size comp ----
  icomp = 1
  fish_size_lk = 0.0
  offset_fs = 0.0
  
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
      # Offset for multinomial likelihood
      offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_size_lk = fish_size_lk - offset_fs
  like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood
  
  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  n_rec = nrow(df)
  yrs_rec = df$year
  pred_rec = mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD
  
  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F50 * slx_fish[A,T])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F40 * slx_fish[A,T])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F35 * slx_fish[A,T])))
  
  # Calculate spawning biomass per recruit for reference points
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
  }
  
  # Calculate reference point spawning biomasses
  SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
  SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
  SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
  SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35
  
  # SPR penalties to enforce reference point constraints
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  # Scale reference points by mean recruitment
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # likelihood/penalties --------------------
  like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  
  # nll ----
  nll = ssqcatch
  nll = nll + like_srv
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec 
  nll = nll + f_regularity
  nll = nll + nll_M
  nll = nll + nll_q
  nll = nll + nll_sigmaR
  nll = nll + sprpen    
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)
  
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

# RTMB model with:
# survey likelihood bias correction
# no input parameter or data changes
srv_like <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  
  RTMB::getAll(pars, data)
  
  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M = exp(log_M)            # Natural mortality
  a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
  a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
  q = exp(log_q)            # Survey catchability
  F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
  F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
  F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass
  
  # Spawning adjustments
  spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
  spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning
  
  # Index values and dimensions
  A = nrow(age_error)                         # Number of ages in model
  A1 = length(ages)                           # Number of ages in comps
  T = sum(catch_ind)                          # Number of fishery years
  Ts = sum(srv_ind)                           # Number of survey years
  Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
  Tsa = sum(srv_age_ind)                      # Number of survey age comp years
  Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
  L = length(length_bins)                     # Number of length bins
  g = 0.00001                                 # Small number to avoid division by zero
  
  # Containers for model outputs
  Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat = rep(0, A)                                              # Initial numbers-at-age
  catch_pred = rep(0, T)                                           # Predicted catch
  srv_pred = rep(0, Ts)                                            # Predicted survey index
  srv_var = rep(0,Ts)                                              # Survey variance
  fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
  srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
  fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
  spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
  N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit
  
  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
  nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
  nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability
  
  # selectivity ----
  slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
  slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
  slx_block[,3] = sel_double_logistic(1:A, a50C[2], deltaC[2], adj=0) # Block 3: double logistic selectivity
  slx_block[,2] = (slx_block[,1] + slx_block[,3]) * 0.5          # Block 2: average of blocks 1 and 3
  slx_block[,4] = to_one(sel_double_logistic(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
  slx_block[,3] = to_one(slx_block[,3])                          # Normalize block 3 selectivity - must be done after block 2 to match ADMB
  
  for(t in 1:T) {
    slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
  }
  
  slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)
  
  
  # mortality ----
  # Calculate fishing mortality for each year
  Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
  for(t in 1:T){
    Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
    Zat[,t] = Fat[,t] + M           # Total mortality at age and year
  }
  Sat = exp(-Zat)                     # Survivorship at age and year
  
  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)
  
  # Remaining columns: survivors from previous year
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
  }
  
  # Calculate recruits and biomasses
  recruits = Nat[1,]                          # Recruitment time series
  spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
  tot_bio = colSums(Nat * waa)                # Total biomass by year
  
  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj = Sat[,T]^(spawn_fract)
  spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
  
  ## catch ----
  # Calculate predicted catch at age and year
  Cat = Fat / Zat * Nat * (1-Sat)
  catch_pred = colSums(Cat * waa) # Predicted catch biomass
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  
  ## survey biomass - bias corrected ----
  isrv = 1
  srv_like = 0.0
  
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
      # survey likelihood (lognormal)
      CV = srv_sd[isrv] / srv_obs[isrv]
      log_sd = sqrt(log(1 + CV^2))
      mu = log(srv_pred[isrv]) - 0.5 * log_sd^2
      srv_like = srv_like - dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
      isrv = isrv + 1
    }
  }
  
  like_srv = srv_like * srv_wt # Weighted survey likelihood
  
  ## fishery age comp ----
  fish_age_lk = 0.0
  offset = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      # Offset for multinomial likelihood
      offset = offset - fish_age_iss[icomp] * 
        sum((fish_age_obs[,icomp] + g) * 
              log(fish_age_obs[,icomp] + g))
      # Multinomial likelihood for age composition
      fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * 
                                        (fish_age_obs[,icomp] + g) * 
                                        log(fish_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_age_lk = fish_age_lk - offset
  like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood
  
  ## survey age comp ----
  srv_age_lk = 0.0
  offset_sa = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
      # Offset for multinomial likelihood
      offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  srv_age_lk = srv_age_lk - offset_sa
  like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood
  
  ## fishery size comp ----
  icomp = 1
  fish_size_lk = 0.0
  offset_fs = 0.0
  
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
      # Offset for multinomial likelihood
      offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_size_lk = fish_size_lk - offset_fs
  like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood
  
  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  n_rec = nrow(df)
  yrs_rec = df$year
  pred_rec = mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD
  
  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F50 * slx_fish[A,T])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F40 * slx_fish[A,T])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F35 * slx_fish[A,T])))
  
  # Calculate spawning biomass per recruit for reference points
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
  }
  
  # Calculate reference point spawning biomasses
  SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
  SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
  SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
  SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35
  
  # SPR penalties to enforce reference point constraints
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  # Scale reference points by mean recruitment
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # likelihood/penalties --------------------
  like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  
  # nll ----
  nll = ssqcatch
  nll = nll + like_srv
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec 
  nll = nll + f_regularity
  nll = nll + nll_M
  nll = nll + nll_q
  nll = nll + nll_sigmaR
  nll = nll + sprpen    
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)
  
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

# RTMB model with:
# survey likelihood bias correction
# selectivity scaled to 1
# no input parameter or data changes
slx_scale <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  
  RTMB::getAll(pars, data)
  
  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M = exp(log_M)            # Natural mortality
  a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
  a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
  q = exp(log_q)            # Survey catchability
  F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
  F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
  F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass
  
  # Spawning adjustments
  spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
  spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning
  
  # Index values and dimensions
  A = nrow(age_error)                         # Number of ages in model
  A1 = length(ages)                           # Number of ages in comps
  T = sum(catch_ind)                          # Number of fishery years
  Ts = sum(srv_ind)                           # Number of survey years
  Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
  Tsa = sum(srv_age_ind)                      # Number of survey age comp years
  Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
  L = length(length_bins)                     # Number of length bins
  g = 0.00001                                 # Small number to avoid division by zero
  
  # Containers for model outputs
  Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat = rep(0, A)                                              # Initial numbers-at-age
  catch_pred = rep(0, T)                                           # Predicted catch
  srv_pred = rep(0, Ts)                                            # Predicted survey index
  srv_var = rep(0,Ts)                                              # Survey variance
  fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
  srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
  fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
  spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
  N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit
  
  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
  nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
  nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability
  
  # selectivity ----
  slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
  slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
  slx_block[,3] = to_one(sel_double_logistic(1:A, a50C[2], deltaC[2], adj=0)) # Block 3: double logistic selectivity
  slx_block[,2] = to_one((slx_block[,1] + slx_block[,3]) * 0.5)          # Block 2: average of blocks 1 and 3
  slx_block[,4] = to_one(sel_double_logistic(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
  
  
  for(t in 1:T) {
    slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
  }
  
  slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)
  
  
  # mortality ----
  # Calculate fishing mortality for each year
  Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
  for(t in 1:T){
    Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
    Zat[,t] = Fat[,t] + M           # Total mortality at age and year
  }
  Sat = exp(-Zat)                     # Survivorship at age and year
  
  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)
  
  # Remaining columns: survivors from previous year
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
  }
  
  # Calculate recruits and biomasses
  recruits = Nat[1,]                          # Recruitment time series
  spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
  tot_bio = colSums(Nat * waa)                # Total biomass by year
  
  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj = Sat[,T]^(spawn_fract)
  spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
  
  ## catch ----
  # Calculate predicted catch at age and year
  Cat = Fat / Zat * Nat * (1-Sat)
  catch_pred = colSums(Cat * waa) # Predicted catch biomass
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  
  ## survey biomass - bias corrected ----
  isrv = 1
  srv_like = 0.0
  
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
      # survey likelihood (lognormal)
      CV = srv_sd[isrv] / srv_obs[isrv]
      log_sd = sqrt(log(1 + CV^2))
      mu = log(srv_pred[isrv]) - 0.5 * log_sd^2
      srv_like = srv_like - dnorm(log(srv_obs[isrv]), mu, log_sd, log = TRUE)
      isrv = isrv + 1
    }
  }
  
  like_srv = srv_like * srv_wt # Weighted survey likelihood
  
  ## fishery age comp ----
  fish_age_lk = 0.0
  offset = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      # Offset for multinomial likelihood
      offset = offset - fish_age_iss[icomp] * 
        sum((fish_age_obs[,icomp] + g) * 
              log(fish_age_obs[,icomp] + g))
      # Multinomial likelihood for age composition
      fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * 
                                        (fish_age_obs[,icomp] + g) * 
                                        log(fish_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_age_lk = fish_age_lk - offset
  like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood
  
  ## survey age comp ----
  srv_age_lk = 0.0
  offset_sa = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
      # Offset for multinomial likelihood
      offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  srv_age_lk = srv_age_lk - offset_sa
  like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood
  
  ## fishery size comp ----
  icomp = 1
  fish_size_lk = 0.0
  offset_fs = 0.0
  
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
      # Offset for multinomial likelihood
      offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_size_lk = fish_size_lk - offset_fs
  like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood
  
  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  n_rec = nrow(df)
  yrs_rec = df$year
  pred_rec = mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD
  
  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F50 * slx_fish[A,T])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F40 * slx_fish[A,T])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F35 * slx_fish[A,T])))
  
  # Calculate spawning biomass per recruit for reference points
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
  }
  
  # Calculate reference point spawning biomasses
  SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
  SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
  SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
  SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35
  
  # SPR penalties to enforce reference point constraints
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  # Scale reference points by mean recruitment
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # likelihood/penalties --------------------
  like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  
  # nll ----
  nll = ssqcatch
  nll = nll + like_srv
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec 
  nll = nll + f_regularity
  nll = nll + nll_M
  nll = nll + nll_q
  nll = nll + nll_sigmaR
  nll = nll + sprpen    
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)
  
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}


# research models ----
# RTMB model with:
# survey likelihood bias correction
# selectivity scaled to 1
# change waa to time-vary and saa has 3 blocks
tv_wt <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  
  RTMB::getAll(pars, data)
  
  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M = exp(log_M)            # Natural mortality
  a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
  a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
  q = exp(log_q)            # Survey catchability
  F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
  F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
  F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass
  
  # Spawning adjustments
  spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
  spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning
  
  # Index values and dimensions
  A = nrow(age_error)                         # Number of ages in model
  A1 = length(ages)                           # Number of ages in comps
  T = sum(catch_ind)                          # Number of fishery years
  Ts = sum(srv_ind)                           # Number of survey years
  Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
  Tsa = sum(srv_age_ind)                      # Number of survey age comp years
  Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
  L = length(length_bins)                     # Number of length bins
  g = 0.00001                                 # Small number to avoid division by zero
  
  # Containers for model outputs
  Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat = rep(0, A)                                              # Initial numbers-at-age
  catch_pred = rep(0, T)                                           # Predicted catch
  srv_pred = rep(0, Ts)                                            # Predicted survey index
  srv_var = rep(0,Ts)                                              # Survey variance
  fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
  srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
  fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
  spawn_bio = tot_bio = numeric(T)                                  # Spawning and total biomass
  N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit
  wt_mature = matrix(0, A, length(unique(waa_ind)))
  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
  nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
  nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability
  
  # selectivity ----
  slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
  slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
  slx_block[,3] = to_one(sel_double_logistic(1:A, a50C[2], deltaC[2], adj=0)) # Block 3: double logistic selectivity
  slx_block[,2] = to_one((slx_block[,1] + slx_block[,3]) * 0.5)          # Block 2: average of blocks 1 and 3
  slx_block[,4] = to_one(sel_double_logistic(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
  
  
  for(t in 1:T) {
    slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
  }
  
  slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)
  
  
  # mortality ----
  # Calculate fishing mortality for each year
  Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
  for(t in 1:T){
    Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
    Zat[,t] = Fat[,t] + M           # Total mortality at age and year
  }
  Sat = exp(-Zat)                     # Survivorship at age and year
  
  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)
  
  # Remaining columns: survivors from previous year
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
  }
  
  # Calculate recruits and biomasses
  wt_mature = waa * maa * 0.5
  recruits = Nat[1,]                          # Recruitment time series
  for(t in 1:T) {
    spawn_bio[t] = sum(Nat[,t] * wt_mature[,waa_ind[t]])        # Spawning biomass by year
    tot_bio[t] = sum(Nat[,t] * waa[,waa_ind[t]])                # Total biomass by year
  }
  
  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj = Sat[,T]^(spawn_fract)
  spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature[,waa_ind[T]])
  
  ## catch ----
  # Calculate predicted catch at age and year
  Cat = Fat / Zat * Nat * (1-Sat)
  for(t in 1:T) {
    catch_pred[t] = sum(Cat[,t] * waa[,waa_ind[t]]) # Predicted catch biomass
  }
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  
  ## survey biomass ----
  isrv = 1
  srv_like = 0.0
  
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa[,waa_ind[t]]) * q # Predicted survey index
      # Survey likelihood (lognormal, using observed and predicted survey biomass)
      CV = srv_sd[isrv] / srv_obs[isrv]
      log_sd = sqrt(log(1 + CV^2))
      srv_like = srv_like - dnorm(log(srv_obs[isrv]), log(srv_pred[isrv]), log_sd, log = TRUE)
      isrv = isrv + 1
    }
  }
  
  like_srv = srv_like * srv_wt # Weighted survey likelihood
  
  ## fishery age comp ----
  fish_age_lk = 0.0
  offset = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      # Offset for multinomial likelihood
      offset = offset - fish_age_iss[icomp] * 
        sum((fish_age_obs[,icomp] + g) * 
              log(fish_age_obs[,icomp] + g))
      # Multinomial likelihood for age composition
      fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * 
                                        (fish_age_obs[,icomp] + g) * 
                                        log(fish_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_age_lk = fish_age_lk - offset
  like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood
  
  ## survey age comp ----
  srv_age_lk = 0.0
  offset_sa = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
      # Offset for multinomial likelihood
      offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  srv_age_lk = srv_age_lk - offset_sa
  like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood
  
  ## fishery size comp ----
  icomp = 1
  fish_size_lk = 0.0
  offset_fs = 0.0
  
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
      # Offset for multinomial likelihood
      offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_size_lk = fish_size_lk - offset_fs
  like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood
  
  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  n_rec = nrow(df)
  yrs_rec = df$year
  pred_rec = mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD
  
  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F50 * slx_fish[A,T])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F40 * slx_fish[A,T])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F35 * slx_fish[A,T])))
  
  # Calculate spawning biomass per recruit for reference points
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
  }
  
  # Calculate reference point spawning biomasses
  SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
  SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
  SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
  SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35
  
  # SPR penalties to enforce reference point constraints
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  # Scale reference points by mean recruitment
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # likelihood/penalties --------------------
  like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  
  # nll ----
  nll = ssqcatch
  nll = nll + like_srv
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec 
  nll = nll + f_regularity
  nll = nll + nll_M
  nll = nll + nll_q
  nll = nll + nll_sigmaR
  nll = nll + sprpen    
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)
  
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

# same as m25.2 but double normal fishery selectivity for last 2 time blocks
slx_scale_norm <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  
  RTMB::getAll(pars, data)
  
  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M = exp(log_M)            # Natural mortality
  a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
  a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
  q = exp(log_q)            # Survey catchability
  F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
  F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
  F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass
  
  # Spawning adjustments
  spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
  spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning
  
  # Index values and dimensions
  A = nrow(age_error)                         # Number of ages in model
  A1 = length(ages)                           # Number of ages in comps
  T = sum(catch_ind)                          # Number of fishery years
  Ts = sum(srv_ind)                           # Number of survey years
  Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
  Tsa = sum(srv_age_ind)                      # Number of survey age comp years
  Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
  L = length(length_bins)                     # Number of length bins
  g = 0.00001                                 # Small number to avoid division by zero
  
  # Containers for model outputs
  Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat = rep(0, A)                                              # Initial numbers-at-age
  catch_pred = rep(0, T)                                           # Predicted catch
  srv_pred = rep(0, Ts)                                            # Predicted survey index
  srv_var = rep(0,Ts)                                              # Survey variance
  fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
  srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
  fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
  spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
  N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit
  
  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
  nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
  nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability
  
  # selectivity ----
  slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
  slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
  slx_block[,3] = to_one(sel_double_normal(1:A, a50C[2], deltaC[2], adj=0)) # Block 3: double logistic selectivity
  slx_block[,2] = to_one((slx_block[,1] + slx_block[,3]) * 0.5)          # Block 2: average of blocks 1 and 3
  slx_block[,4] = to_one(sel_double_normal(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
  
  
  for(t in 1:T) {
    slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
  }
  
  slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)
  
  
  # mortality ----
  # Calculate fishing mortality for each year
  Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
  for(t in 1:T){
    Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
    Zat[,t] = Fat[,t] + M           # Total mortality at age and year
  }
  Sat = exp(-Zat)                     # Survivorship at age and year
  
  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)
  
  # Remaining columns: survivors from previous year
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
  }
  
  # Calculate recruits and biomasses
  recruits = Nat[1,]                          # Recruitment time series
  spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
  tot_bio = colSums(Nat * waa)                # Total biomass by year
  
  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj = Sat[,T]^(spawn_fract)
  spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
  
  ## catch ----
  # Calculate predicted catch at age and year
  Cat = Fat / Zat * Nat * (1-Sat)
  catch_pred = colSums(Cat * waa) # Predicted catch biomass
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  
  ## survey biomass ----
  isrv = 1
  srv_like = 0.0
  
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
      # Survey likelihood (lognormal, using observed and predicted survey biomass)
      CV = srv_sd[isrv] / srv_obs[isrv]
      log_sd = sqrt(log(1 + CV^2))
      srv_like = srv_like - dnorm(log(srv_obs[isrv]), log(srv_pred[isrv]), log_sd, log = TRUE)
      isrv = isrv + 1
    }
  }
  
  like_srv = srv_like * srv_wt # Weighted survey likelihood
  
  ## fishery age comp ----
  fish_age_lk = 0.0
  offset = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      # Offset for multinomial likelihood
      offset = offset - fish_age_iss[icomp] * 
        sum((fish_age_obs[,icomp] + g) * 
              log(fish_age_obs[,icomp] + g))
      # Multinomial likelihood for age composition
      fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * 
                                        (fish_age_obs[,icomp] + g) * 
                                        log(fish_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_age_lk = fish_age_lk - offset
  like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood
  
  ## survey age comp ----
  srv_age_lk = 0.0
  offset_sa = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
      # Offset for multinomial likelihood
      offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  srv_age_lk = srv_age_lk - offset_sa
  like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood
  
  ## fishery size comp ----
  icomp = 1
  fish_size_lk = 0.0
  offset_fs = 0.0
  
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
      # Offset for multinomial likelihood
      offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_size_lk = fish_size_lk - offset_fs
  like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood
  
  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  n_rec = nrow(df)
  yrs_rec = df$year
  pred_rec = mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD
  
  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F50 * slx_fish[A,T])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F40 * slx_fish[A,T])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F35 * slx_fish[A,T])))
  
  # Calculate spawning biomass per recruit for reference points
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
  }
  
  # Calculate reference point spawning biomasses
  SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
  SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
  SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
  SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35
  
  # SPR penalties to enforce reference point constraints
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  # Scale reference points by mean recruitment
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # likelihood/penalties --------------------
  like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  
  # nll ----
  nll = ssqcatch
  nll = nll + like_srv
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec 
  nll = nll + f_regularity
  nll = nll + nll_M
  nll = nll + nll_q
  nll = nll + nll_sigmaR
  nll = nll + sprpen    
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)
  
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

# same as m25.2 but logistic fishery selectivity for last 2 time blocks
slx_scale_log <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  
  RTMB::getAll(pars, data)
  
  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M = exp(log_M)            # Natural mortality
  a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
  a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
  q = exp(log_q)            # Survey catchability
  F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
  F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
  F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass
  
  # Spawning adjustments
  spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
  spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning
  
  # Index values and dimensions
  A = nrow(age_error)                         # Number of ages in model
  A1 = length(ages)                           # Number of ages in comps
  T = sum(catch_ind)                          # Number of fishery years
  Ts = sum(srv_ind)                           # Number of survey years
  Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
  Tsa = sum(srv_age_ind)                      # Number of survey age comp years
  Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
  L = length(length_bins)                     # Number of length bins
  g = 0.00001                                 # Small number to avoid division by zero
  
  # Containers for model outputs
  Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat = rep(0, A)                                              # Initial numbers-at-age
  catch_pred = rep(0, T)                                           # Predicted catch
  srv_pred = rep(0, Ts)                                            # Predicted survey index
  srv_var = rep(0,Ts)                                              # Survey variance
  fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
  srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
  fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
  spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
  N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit
  
  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
  nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
  nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability
  
  # selectivity ----
  slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
  slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
  slx_block[,3] = sel_logistic(1:A, a50C[2], deltaC[2], adj=0) # Block 3: double logistic selectivity
  slx_block[,2] = to_one((slx_block[,1] + slx_block[,3]) * 0.5)          # Block 2: average of blocks 1 and 3
  slx_block[,4] = sel_logistic(1:A, a50C[3], deltaC[3], adj=0) # Block 4: double logistic selectivity
  
  
  for(t in 1:T) {
    slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
  }
  
  slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)
  
  
  # mortality ----
  # Calculate fishing mortality for each year
  Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
  for(t in 1:T){
    Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
    Zat[,t] = Fat[,t] + M           # Total mortality at age and year
  }
  Sat = exp(-Zat)                     # Survivorship at age and year
  
  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)
  
  # Remaining columns: survivors from previous year
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
  }
  
  # Calculate recruits and biomasses
  recruits = Nat[1,]                          # Recruitment time series
  spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
  tot_bio = colSums(Nat * waa)                # Total biomass by year
  
  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj = Sat[,T]^(spawn_fract)
  spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
  
  ## catch ----
  # Calculate predicted catch at age and year
  Cat = Fat / Zat * Nat * (1-Sat)
  catch_pred = colSums(Cat * waa) # Predicted catch biomass
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  
  ## survey biomass ----
  isrv = 1
  srv_like = 0.0
  
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
      # Survey likelihood (lognormal, using observed and predicted survey biomass)
      CV = srv_sd[isrv] / srv_obs[isrv]
      log_sd = sqrt(log(1 + CV^2))
      srv_like = srv_like - dnorm(log(srv_obs[isrv]), log(srv_pred[isrv]), log_sd, log = TRUE)
      isrv = isrv + 1
    }
  }
  
  like_srv = srv_like * srv_wt # Weighted survey likelihood
  
  ## fishery age comp ----
  fish_age_lk = 0.0
  offset = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      # Offset for multinomial likelihood
      offset = offset - fish_age_iss[icomp] * 
        sum((fish_age_obs[,icomp] + g) * 
              log(fish_age_obs[,icomp] + g))
      # Multinomial likelihood for age composition
      fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * 
                                        (fish_age_obs[,icomp] + g) * 
                                        log(fish_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_age_lk = fish_age_lk - offset
  like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood
  
  ## survey age comp ----
  srv_age_lk = 0.0
  offset_sa = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
      # Offset for multinomial likelihood
      offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  srv_age_lk = srv_age_lk - offset_sa
  like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood
  
  ## fishery size comp ----
  icomp = 1
  fish_size_lk = 0.0
  offset_fs = 0.0
  
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
      # Offset for multinomial likelihood
      offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_size_lk = fish_size_lk - offset_fs
  like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood
  
  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  n_rec = nrow(df)
  yrs_rec = df$year
  pred_rec = mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD
  
  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F50 * slx_fish[A,T])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F40 * slx_fish[A,T])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F35 * slx_fish[A,T])))
  
  # Calculate spawning biomass per recruit for reference points
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
  }
  
  # Calculate reference point spawning biomasses
  SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
  SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
  SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
  SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35
  
  # SPR penalties to enforce reference point constraints
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  # Scale reference points by mean recruitment
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # likelihood/penalties --------------------
  like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  
  # nll ----
  nll = ssqcatch
  nll = nll + like_srv
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec 
  nll = nll + f_regularity
  nll = nll + nll_M
  nll = nll + nll_q
  nll = nll + nll_sigmaR
  nll = nll + sprpen    
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)
  
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

# slx_scale <- function(pars, data) {
#   require(RTMB)
#   "c" <- RTMB::ADoverload("c")
#   "[<-" <- RTMB::ADoverload("[<-")
#   
#   RTMB::getAll(pars, data)
#   
#   # setup -------------
#   # transform
#   # Exponentiate log-parameters to get values on natural scale
#   M = exp(log_M)            # Natural mortality
#   a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
#   a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
#   q = exp(log_q)            # Survey catchability
#   F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
#   F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
#   F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass
#   
#   # Spawning adjustments
#   spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
#   spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning
#   
#   # Index values and dimensions
#   A = nrow(age_error)                         # Number of ages in model
#   A1 = length(ages)                           # Number of ages in comps
#   T = sum(catch_ind)                          # Number of fishery years
#   Ts = sum(srv_ind)                           # Number of survey years
#   Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
#   Tsa = sum(srv_age_ind)                      # Number of survey age comp years
#   Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
#   L = length(length_bins)                     # Number of length bins
#   g = 0.00001                                 # Small number to avoid division by zero
#   
#   # Containers for model outputs
#   Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
#   initNat = rep(0, A)                                              # Initial numbers-at-age
#   catch_pred = rep(0, T)                                           # Predicted catch
#   srv_pred = rep(0, Ts)                                            # Predicted survey index
#   srv_var = rep(0,Ts)                                              # Survey variance
#   fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
#   srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
#   fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
#   spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
#   N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit
#   
#   # priors -----------------
#   # Priors on key parameters (negative log-likelihood contributions)
#   nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
#   nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
#   nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability
#   
#   # selectivity ----
#   slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
#   slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
#   slx_block[,3] = to_one(sel_double_logistic(1:A, a50C[2], deltaC[2], adj=0)) # Block 3: double logistic selectivity
#   slx_block[,2] = to_one((slx_block[,1] + slx_block[,3]) * 0.5)          # Block 2: average of blocks 1 and 3
#   slx_block[,4] = to_one(sel_double_logistic(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
#   
#   
#   for(t in 1:T) {
#     slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
#   }
#   
#   slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)
#   
#   
#   # mortality ----
#   # Calculate fishing mortality for each year
#   Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
#   for(t in 1:T){
#     Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
#     Zat[,t] = Fat[,t] + M           # Total mortality at age and year
#   }
#   Sat = exp(-Zat)                     # Survivorship at age and year
#   
#   ## Nat ----
#   # Populate numbers-at-age matrix (Nat)
#   # First row: recruitment for each year
#   # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
#   for(t in 1:T) {
#     Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
#   }
#   # First column: initial numbers-at-age for each cohort
#   for(a in 2:(A-1)) {
#     Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
#   }
#   Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)
#   
#   # Remaining columns: survivors from previous year
#   for(t in 2:T) {
#     for(a in 2:A) {
#       Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
#     }
#     Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
#   }
#   
#   # Calculate recruits and biomasses
#   recruits = Nat[1,]                          # Recruitment time series
#   spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
#   tot_bio = colSums(Nat * waa)                # Total biomass by year
#   
#   # Adjust spawning biomass in last year for pre-spawning mortality
#   spawn_adj = Sat[,T]^(spawn_fract)
#   spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
#   
#   ## catch ----
#   # Calculate predicted catch at age and year
#   Cat = Fat / Zat * Nat * (1-Sat)
#   catch_pred = colSums(Cat * waa) # Predicted catch biomass
#   ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
#   
#   ## survey biomass ----
#   isrv = 1
#   srv_like = 0.0
#   
#   for(t in 1:T) {
#     if(srv_ind[t]==1) {
#       srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
#       # Survey likelihood (lognormal, using observed and predicted survey biomass)
#       CV = srv_sd[isrv] / srv_obs[isrv]
#       log_sd = sqrt(log(1 + CV^2))
#       srv_like = srv_like - dnorm(log(srv_obs[isrv]), log(srv_pred[isrv]), log_sd, log = TRUE)
#       isrv = isrv + 1
#     }
#   }
#   
#   like_srv = srv_like * srv_wt # Weighted survey likelihood
#   
#   ## fishery age comp ----
#   fish_age_lk = 0.0
#   offset = 0.0
#   icomp = 1
#   
#   for(t in 1:T) {
#     if(fish_age_ind[t] == 1) {
#       # Predicted age composition (with ageing error)
#       fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
#       # Offset for multinomial likelihood
#       offset = offset - fish_age_iss[icomp] * 
#         sum((fish_age_obs[,icomp] + g) * 
#               log(fish_age_obs[,icomp] + g))
#       # Multinomial likelihood for age composition
#       fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * 
#                                         (fish_age_obs[,icomp] + g) * 
#                                         log(fish_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_age_lk = fish_age_lk - offset
#   like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood
#   
#   ## survey age comp ----
#   srv_age_lk = 0.0
#   offset_sa = 0.0
#   icomp = 1
#   
#   for(t in 1:T) {
#     if(srv_age_ind[t] == 1) {
#       # Predicted survey age composition (with ageing error)
#       srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
#       # Offset for multinomial likelihood
#       offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
#       # Multinomial likelihood for survey age composition
#       srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   srv_age_lk = srv_age_lk - offset_sa
#   like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood
#   
#   ## fishery size comp ----
#   icomp = 1
#   fish_size_lk = 0.0
#   offset_fs = 0.0
#   
#   for(t in 1:T) {
#     if(fish_size_ind[t] == 1) {
#       # Predicted size composition (with size-at-age array)
#       fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
#       # Offset for multinomial likelihood
#       offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
#       # Multinomial likelihood for size composition
#       fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
#       icomp = icomp + 1
#     }
#   }
#   fish_size_lk = fish_size_lk - offset_fs
#   like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood
#   
#   # SPR ------------------------
#   # Prepare recruitment data frame for reference point calculations
#   data.frame(log_Rt = log_Rt,
#              pred_rec = Nat[1,],
#              year = years) -> df
#   # Filter years for recruitment estimation (exclude first and last ages)
#   df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
#   n_rec = nrow(df)
#   yrs_rec = df$year
#   pred_rec = mean(df$pred_rec) # Mean predicted recruitment
#   stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD
#   
#   # Calculate numbers per recruit for reference points (F50, F40, F35)
#   for(a in 2:A) {
#     N_spr[a,1] = N_spr[a-1,1] * exp(-M)
#     N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
#     N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
#     N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
#   }
#   # Plus group for per-recruit calculations
#   N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
#   N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) / 
#     (1 - exp(-(M + F50 * slx_fish[A,T])))
#   N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) / 
#     (1 - exp(-(M + F40 * slx_fish[A,T])))
#   N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) / 
#     (1 - exp(-(M + F35 * slx_fish[A,T])))
#   
#   # Calculate spawning biomass per recruit for reference points
#   for(a in 1:A) {
#     sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
#     sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
#     sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
#     sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
#   }
#   
#   # Calculate reference point spawning biomasses
#   SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
#   SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
#   SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
#   SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35
#   
#   # SPR penalties to enforce reference point constraints
#   sprpen = 100. * (SBF50 / SB0 - 0.5)^2
#   sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
#   sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
#   
#   # Scale reference points by mean recruitment
#   B0 = SB0 * pred_rec
#   B40 = SBF40 * pred_rec
#   B35 = SBF35 * pred_rec
#   
#   # likelihood/penalties --------------------
#   like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
#   f_regularity = wt_fmort_reg * sum(log_Ft^2)
#   
#   # nll ----
#   nll = ssqcatch
#   nll = nll + like_srv
#   nll = nll + like_fish_age
#   nll = nll + like_srv_age
#   nll = nll + like_fish_size
#   nll = nll + like_rec 
#   nll = nll + f_regularity
#   nll = nll + nll_M
#   nll = nll + nll_q
#   nll = nll + nll_sigmaR
#   nll = nll + sprpen    
#   
#   # reports -------------------
#   RTMB::REPORT(ages)
#   RTMB::REPORT(years)
#   RTMB::REPORT(M)
#   RTMB::ADREPORT(M)
#   RTMB::REPORT(a50C)
#   RTMB::REPORT(deltaC)
#   RTMB::REPORT(a50S)
#   RTMB::REPORT(deltaS)
#   RTMB::REPORT(q)
#   RTMB::ADREPORT(q)
#   RTMB::REPORT(sigmaR)
#   RTMB::REPORT(log_mean_R)
#   RTMB::REPORT(log_Rt)
#   RTMB::ADREPORT(log_Rt)
#   RTMB::REPORT(log_mean_F)
#   RTMB::REPORT(log_Ft)
#   RTMB::REPORT(waa)
#   RTMB::REPORT(maa)
#   RTMB::REPORT(wt_mature)
#   RTMB::REPORT(yield_ratio)
#   RTMB::REPORT(Fat)
#   RTMB::REPORT(Zat)
#   RTMB::REPORT(Sat)
#   RTMB::REPORT(Cat)
#   RTMB::REPORT(Nat)
#   RTMB::REPORT(slx_srv)
#   RTMB::REPORT(slx_fish)
#   RTMB::REPORT(slx_block)
#   RTMB::REPORT(Ft)
#   RTMB::REPORT(catch_pred)
#   RTMB::REPORT(srv_pred)
#   
#   RTMB::REPORT(fish_age_pred)
#   RTMB::REPORT(srv_age_pred)
#   RTMB::REPORT(fish_size_pred)
#   
#   RTMB::REPORT(tot_bio)
#   RTMB::REPORT(spawn_bio)
#   RTMB::REPORT(recruits)
#   RTMB::ADREPORT(srv_pred)
#   RTMB::ADREPORT(tot_bio)
#   RTMB::ADREPORT(spawn_bio)
#   RTMB::ADREPORT(recruits)
#   RTMB::REPORT(spawn_fract)
#   RTMB::REPORT(B0)
#   RTMB::REPORT(B40)
#   RTMB::REPORT(B35)
#   RTMB::REPORT(F35)
#   RTMB::REPORT(F40)
#   RTMB::REPORT(F50)
#   RTMB::REPORT(pred_rec)
#   RTMB::REPORT(n_rec)
#   RTMB::REPORT(yrs_rec)
#   RTMB::REPORT(stdev_rec)
#   
#   RTMB::REPORT(ssqcatch)
#   RTMB::REPORT(like_srv)
#   RTMB::REPORT(like_fish_age)
#   RTMB::REPORT(like_srv_age)
#   RTMB::REPORT(like_fish_size)
#   RTMB::REPORT(like_rec)
#   RTMB::REPORT(f_regularity)
#   RTMB::REPORT(sprpen)
#   RTMB::REPORT(nll_q)
#   RTMB::REPORT(nll_M)
#   RTMB::REPORT(nll_sigmaR)
#   RTMB::REPORT(nll)
#   # nll = 0.0
#   return(nll)
# }
# random walk on last 2 fishery selectivity time block
# data changes:
#     change fish_block_ind = 4 <- 3
#     n_blk = sum(fish_block_ind==3)
# input parameter changes:
#     remove 3rd log_a50C and deltaC
#   add the following pars
#    log_a50C_rw = rep(1, n_blk)
#    deltaC_rw = rep(5, n_blk)
#    sigma_rw = 0.1
slx_fish_rw <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  
  RTMB::getAll(pars, data)
  
  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M = exp(log_M)            # Natural mortality
  a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
  a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
  q = exp(log_q)            # Survey catchability
  F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
  F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
  F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass
  
  # Spawning adjustments
  spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
  spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning
  
  # Index values and dimensions
  A = nrow(age_error)                         # Number of ages in model
  A1 = length(ages)                           # Number of ages in comps
  T = sum(catch_ind)                          # Number of fishery years
  Ts = sum(srv_ind)                           # Number of survey years
  Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
  Tsa = sum(srv_age_ind)                      # Number of survey age comp years
  Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
  L = length(length_bins)                     # Number of length bins
  g = 0.00001                                 # Small number to avoid division by zero
  
  # Containers for model outputs
  Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat = rep(0, A)                                              # Initial numbers-at-age
  catch_pred = rep(0, T)                                           # Predicted catch
  srv_pred = rep(0, Ts)                                            # Predicted survey index
  srv_var = rep(0,Ts)                                              # Survey variance
  fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
  srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
  fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
  spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
  N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit
  
  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
  nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
  nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability
  
  # selectivity ----
  slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
  slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
  slx_block[,2] = to_one(sel_double_logistic(1:A, a50C[2], deltaC[2], adj=0)) # Block 3: double logistic selectivity
  # slx_block[,2] = to_one((slx_block[,1] + slx_block[,3]) * 0.5)          # Block 2: average of blocks 1 and 3
  # slx_block[,4] = to_one(sel_double_logistic(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
  
  nll_rw_a50C = 0
  for (t in 2:length(log_a50C_rw)) {
    nll_rw_a50C = nll_rw_a50C + (log_a50C_rw[t] - log_a50C_rw[t-1])^2 / (2 * sigma_rw^2)
  }
  
  rw_index = 1
  for(t in 1:T) {
    if (fish_block_ind[t] == 3) {
      a50 = exp(log_a50C_rw[rw_index])
      slx_fish[,t] = to_one(sel_double_logistic(1:A, a50, deltaC[3], adj=0))
      rw_index = rw_index + 1
    } else {
      slx_fish[,t] = slx_block[,fish_block_ind[t]]
    }
  }
  
  # for(t in 1:T) {
  #   slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
  # }
  
  slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)
  
  
  # mortality ----
  # Calculate fishing mortality for each year
  Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
  for(t in 1:T){
    Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
    Zat[,t] = Fat[,t] + M           # Total mortality at age and year
  }
  Sat = exp(-Zat)                     # Survivorship at age and year
  
  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)
  
  # Remaining columns: survivors from previous year
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
  }
  
  # Calculate recruits and biomasses
  recruits = Nat[1,]                          # Recruitment time series
  spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
  tot_bio = colSums(Nat * waa)                # Total biomass by year
  
  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj = Sat[,T]^(spawn_fract)
  spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
  
  ## catch ----
  # Calculate predicted catch at age and year
  Cat = Fat / Zat * Nat * (1-Sat)
  catch_pred = colSums(Cat * waa) # Predicted catch biomass
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  
  ## survey biomass ----
  isrv = 1
  srv_like = 0.0
  
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
      # Survey likelihood (lognormal, using observed and predicted survey biomass)
      CV = srv_sd[isrv] / srv_obs[isrv]
      log_sd = sqrt(log(1 + CV^2))
      srv_like = srv_like - dnorm(log(srv_obs[isrv]), log(srv_pred[isrv]), log_sd, log = TRUE)
      isrv = isrv + 1
    }
  }
  
  like_srv = srv_like * srv_wt # Weighted survey likelihood
  
  ## fishery age comp ----
  fish_age_lk = 0.0
  offset = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      # Offset for multinomial likelihood
      offset = offset - fish_age_iss[icomp] * 
        sum((fish_age_obs[,icomp] + g) * 
              log(fish_age_obs[,icomp] + g))
      # Multinomial likelihood for age composition
      fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * 
                                        (fish_age_obs[,icomp] + g) * 
                                        log(fish_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_age_lk = fish_age_lk - offset
  like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood
  
  ## survey age comp ----
  srv_age_lk = 0.0
  offset_sa = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
      # Offset for multinomial likelihood
      offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  srv_age_lk = srv_age_lk - offset_sa
  like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood
  
  ## fishery size comp ----
  icomp = 1
  fish_size_lk = 0.0
  offset_fs = 0.0
  
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
      # Offset for multinomial likelihood
      offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_size_lk = fish_size_lk - offset_fs
  like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood
  
  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  n_rec = nrow(df)
  yrs_rec = df$year
  pred_rec = mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD
  
  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F50 * slx_fish[A,T])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F40 * slx_fish[A,T])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F35 * slx_fish[A,T])))
  
  # Calculate spawning biomass per recruit for reference points
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
  }
  
  # Calculate reference point spawning biomasses
  SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
  SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
  SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
  SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35
  
  # SPR penalties to enforce reference point constraints
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  # Scale reference points by mean recruitment
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # likelihood/penalties --------------------
  like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  
  # nll ----
  nll = ssqcatch
  nll = nll + like_srv
  nll = nll + nll_rw_a50C
  # nll = nll + nll_rw_deltaC
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec 
  nll = nll + f_regularity
  nll = nll + nll_M
  nll = nll + nll_q
  nll = nll + nll_sigmaR
  nll = nll + sprpen    
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(log_a50C_rw)
  RTMB::ADREPORT(log_a50C_rw)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)
  
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(sigma_rw)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

# random walk on last 2 fishery selectivity time block
# time varying q
slx_fish_rw_tvq <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  
  RTMB::getAll(pars, data)
  
  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M = exp(log_M)            # Natural mortality
  a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
  a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
  q = exp(log_q)            # Survey catchability
  F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
  F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
  F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass
  
  # Spawning adjustments
  spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
  spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning
  
  # Index values and dimensions
  A = nrow(age_error)                         # Number of ages in model
  A1 = length(ages)                           # Number of ages in comps
  T = sum(catch_ind)                          # Number of fishery years
  Ts = sum(srv_ind)                           # Number of survey years
  Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
  Tsa = sum(srv_age_ind)                      # Number of survey age comp years
  Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
  L = length(length_bins)                     # Number of length bins
  g = 0.00001                                 # Small number to avoid division by zero
  
  # Containers for model outputs
  Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat = rep(0, A)                                              # Initial numbers-at-age
  catch_pred = rep(0, T)                                           # Predicted catch
  srv_pred = rep(0, Ts)                                            # Predicted survey index
  srv_var = rep(0,Ts)                                              # Survey variance
  fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
  srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
  fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
  spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
  N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit
  
  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
  nll_q = 0
  for(t in 2:length(log_q)) {
    nll_q = nll_q + (log_q[t] - log_q[t - 1])^2 / (2 * sigma_q^2)
  }        # Prior on survey catchability
  nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability
  
  # selectivity ----
  slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
  slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
  slx_block[,2] = to_one(sel_double_logistic(1:A, a50C[2], deltaC[2], adj=0)) # Block 3: double logistic selectivity
  # slx_block[,2] = to_one((slx_block[,1] + slx_block[,3]) * 0.5)          # Block 2: average of blocks 1 and 3
  # slx_block[,4] = to_one(sel_double_logistic(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
  
  nll_rw_a50C = 0
  for (t in 2:length(log_a50C_rw)) {
    nll_rw_a50C = nll_rw_a50C + (log_a50C_rw[t] - log_a50C_rw[t-1])^2 / (2 * sigma_rw^2)
  }
  
  rw_index = 1
  for(t in 1:T) {
    if (fish_block_ind[t] == 3) {
      a50 = exp(log_a50C_rw[rw_index])
      # delta = deltaC_rw[rw_index]  
      slx_fish[,t] = to_one(sel_double_logistic(1:A, a50, deltaC[3], adj=0))
      rw_index = rw_index + 1
    } else {
      slx_fish[,t] = slx_block[,fish_block_ind[t]]
    }
  }
  
  # for(t in 1:T) {
  #   slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
  # }
  
  slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)
  
  
  # mortality ----
  # Calculate fishing mortality for each year
  Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
  for(t in 1:T){
    Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
    Zat[,t] = Fat[,t] + M           # Total mortality at age and year
  }
  Sat = exp(-Zat)                     # Survivorship at age and year
  
  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)
  
  # Remaining columns: survivors from previous year
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
  }
  
  # Calculate recruits and biomasses
  recruits = Nat[1,]                          # Recruitment time series
  spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
  tot_bio = colSums(Nat * waa)                # Total biomass by year
  
  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj = Sat[,T]^(spawn_fract)
  spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
  
  ## catch ----
  # Calculate predicted catch at age and year
  Cat = Fat / Zat * Nat * (1-Sat)
  catch_pred = colSums(Cat * waa) # Predicted catch biomass
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  
  ## survey biomass ----
  isrv = 1
  srv_like = 0.0
  
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q[isrv] # Predicted survey index
      # Survey likelihood (lognormal, using observed and predicted survey biomass)
      CV = srv_sd[isrv] / srv_obs[isrv]
      log_sd = sqrt(log(1 + CV^2))
      srv_like = srv_like - dnorm(log(srv_obs[isrv]), log(srv_pred[isrv]), log_sd, log = TRUE)
      isrv = isrv + 1
    }
  }
  
  like_srv = srv_like * srv_wt # Weighted survey likelihood
  
  ## fishery age comp ----
  fish_age_lk = 0.0
  offset = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      # Offset for multinomial likelihood
      offset = offset - fish_age_iss[icomp] * 
        sum((fish_age_obs[,icomp] + g) * 
              log(fish_age_obs[,icomp] + g))
      # Multinomial likelihood for age composition
      fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * 
                                        (fish_age_obs[,icomp] + g) * 
                                        log(fish_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_age_lk = fish_age_lk - offset
  like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood
  
  ## survey age comp ----
  srv_age_lk = 0.0
  offset_sa = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
      # Offset for multinomial likelihood
      offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  srv_age_lk = srv_age_lk - offset_sa
  like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood
  
  ## fishery size comp ----
  icomp = 1
  fish_size_lk = 0.0
  offset_fs = 0.0
  
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
      # Offset for multinomial likelihood
      offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_size_lk = fish_size_lk - offset_fs
  like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood
  
  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  n_rec = nrow(df)
  yrs_rec = df$year
  pred_rec = mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD
  
  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F50 * slx_fish[A,T])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F40 * slx_fish[A,T])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F35 * slx_fish[A,T])))
  
  # Calculate spawning biomass per recruit for reference points
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
  }
  
  # Calculate reference point spawning biomasses
  SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
  SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
  SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
  SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35
  
  # SPR penalties to enforce reference point constraints
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  # Scale reference points by mean recruitment
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # likelihood/penalties --------------------
  like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  
  # nll ----
  nll = ssqcatch
  nll = nll + like_srv
  nll = nll + nll_rw_a50C
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec 
  nll = nll + f_regularity
  nll = nll + nll_M
  nll = nll + nll_q
  nll = nll + nll_sigmaR
  nll = nll + sprpen    
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(log_a50C_rw)
  RTMB::ADREPORT(log_a50C_rw)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)
  
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(sigma_rw)
  RTMB::REPORT(sigma_q)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}


# survey slx random wawlk on a50
slx_srv <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  
  RTMB::getAll(pars, data)
  
  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M = exp(log_M)            # Natural mortality
  a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
  a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
  q = exp(log_q)            # Survey catchability
  F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
  F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
  F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass
  
  # Spawning adjustments
  spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
  spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning
  
  # Index values and dimensions
  A = nrow(age_error)                         # Number of ages in model
  A1 = length(ages)                           # Number of ages in comps
  T = sum(catch_ind)                          # Number of fishery years
  Ts = sum(srv_ind)                           # Number of survey years
  Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
  Tsa = sum(srv_age_ind)                      # Number of survey age comp years
  Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
  L = length(length_bins)                     # Number of length bins
  g = 0.00001                                 # Small number to avoid division by zero
  
  # Containers for model outputs
  Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat = rep(0, A)                                              # Initial numbers-at-age
  catch_pred = rep(0, T)                                           # Predicted catch
  srv_pred = rep(0, Ts)                                            # Predicted survey index
  srv_var = rep(0,Ts)                                              # Survey variance
  fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
  srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
  fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
  spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
  N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit
  
  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
  nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
  nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability
  
  # selectivity ----
  slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
  slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
  slx_block[,3] = to_one(sel_double_logistic(1:A, a50C[2], deltaC[2], adj=0)) # Block 3: double logistic selectivity
  slx_block[,2] = to_one((slx_block[,1] + slx_block[,3]) * 0.5)          # Block 2: average of blocks 1 and 3
  slx_block[,4] = to_one(sel_double_logistic(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
  
  
  for(t in 1:T) {
    slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
  }
  
  slx_srv = matrix(0, A, Ts)
  for (t in 1:Ts) {
    slx_srv[,t] = sel_logistic(1:A, a50S[t], deltaS, adj=0)
  }
  
  # slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)
  
  
  # mortality ----
  # Calculate fishing mortality for each year
  Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
  for(t in 1:T){
    Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
    Zat[,t] = Fat[,t] + M           # Total mortality at age and year
  }
  Sat = exp(-Zat)                     # Survivorship at age and year
  
  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)
  
  # Remaining columns: survivors from previous year
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
  }
  
  # Calculate recruits and biomasses
  recruits = Nat[1,]                          # Recruitment time series
  spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
  tot_bio = colSums(Nat * waa)                # Total biomass by year
  
  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj = Sat[,T]^(spawn_fract)
  spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
  
  ## catch ----
  # Calculate predicted catch at age and year
  Cat = Fat / Zat * Nat * (1-Sat)
  catch_pred = colSums(Cat * waa) # Predicted catch biomass
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  
  ## survey biomass ----
  isrv = 1
  srv_like = 0.0
  
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      srv_pred[isrv] = sum(Nat[,t] * slx_srv[,isrv] * waa) * q # Predicted survey index
      # Survey likelihood (lognormal, using observed and predicted survey biomass)
      CV = srv_sd[isrv] / srv_obs[isrv]
      log_sd = sqrt(log(1 + CV^2))
      srv_like = srv_like - dnorm(log(srv_obs[isrv]), log(srv_pred[isrv]), log_sd, log = TRUE)
      isrv = isrv + 1
    }
  }
  
  like_srv = srv_like * srv_wt # Weighted survey likelihood
  
  ## fishery age comp ----
  fish_age_lk = 0.0
  offset = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      # Offset for multinomial likelihood
      offset = offset - fish_age_iss[icomp] * 
        sum((fish_age_obs[,icomp] + g) * 
              log(fish_age_obs[,icomp] + g))
      # Multinomial likelihood for age composition
      fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * 
                                        (fish_age_obs[,icomp] + g) * 
                                        log(fish_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_age_lk = fish_age_lk - offset
  like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood
  
  ## survey age comp ----
  srv_age_lk = 0.0
  offset_sa = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv[,icomp]) / sum(Nat[,t] * slx_srv[,icomp]) * age_error))
      # Offset for multinomial likelihood
      offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  srv_age_lk = srv_age_lk - offset_sa
  like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood
  
  ## fishery size comp ----
  icomp = 1
  fish_size_lk = 0.0
  offset_fs = 0.0
  
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
      # Offset for multinomial likelihood
      offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_size_lk = fish_size_lk - offset_fs
  like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood
  
  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  n_rec = nrow(df)
  yrs_rec = df$year
  pred_rec = mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD
  
  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F50 * slx_fish[A,T])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F40 * slx_fish[A,T])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F35 * slx_fish[A,T])))
  
  # Calculate spawning biomass per recruit for reference points
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
  }
  
  # Calculate reference point spawning biomasses
  SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
  SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
  SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
  SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35
  
  # SPR penalties to enforce reference point constraints
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  # Scale reference points by mean recruitment
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # likelihood/penalties --------------------
  like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  nll_slx_srv = 0
  for (t in 2:Ts) {
    nll_slx_srv = nll_slx_srv + (log_a50S[t] - log_a50S[t-1])^2 / (2 * sigma_slx_srv^2)
  }
  # nll ----
  nll = 0.0
  nll = nll + ssqcatch
  nll = nll + nll_slx_srv
  nll = nll + like_srv
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec 
  nll = nll + f_regularity
  nll = nll + nll_M
  nll = nll + nll_q
  nll = nll + nll_sigmaR
  nll = nll + sprpen    
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)
  
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

# RTMB model with:
# survey likelihood bias correction
# selectivity scaled to 1
# with cohort smoothing
cohort <- function(pars, data) {
  require(RTMB)
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  
  RTMB::getAll(pars, data)
  
  # setup -------------
  # transform
  # Exponentiate log-parameters to get values on natural scale
  M = exp(log_M)            # Natural mortality
  a50C = exp(log_a50C)      # Age at 50% selectivity (fishery)
  a50S = exp(log_a50S)      # Age at 50% selectivity (survey)
  q = exp(log_q)            # Survey catchability
  F50 = exp(log_F50)        # Fishing mortality at 50% spawning biomass
  F40 = exp(log_F40)        # Fishing mortality at 40% spawning biomass
  F35 = exp(log_F35)        # Fishing mortality at 35% spawning biomass
  
  # Spawning adjustments
  spawn_fract = (spawn_mo - 1) / 12           # Fraction of year before spawning
  spawn_adj = exp(-M)^(spawn_fract)           # Mortality adjustment for spawning
  
  # Index values and dimensions
  A = nrow(age_error)                         # Number of ages in model
  A1 = length(ages)                           # Number of ages in comps
  T = sum(catch_ind)                          # Number of fishery years
  Ts = sum(srv_ind)                           # Number of survey years
  Tfa = sum(fish_age_ind)                     # Number of fishery age comp years
  Tsa = sum(srv_age_ind)                      # Number of survey age comp years
  Tfs = sum(fish_size_ind)                    # Number of fishery size comp years
  L = length(length_bins)                     # Number of length bins
  g = 0.00001                                 # Small number to avoid division by zero
  
  # Containers for model outputs
  Bat = Cat = Nat = Fat = Zat = Sat = slx_fish = matrix(0, A, T)   # Matrices for biomass, catch, numbers, F, Z, S, selectivity
  initNat = rep(0, A)                                              # Initial numbers-at-age
  catch_pred = rep(0, T)                                           # Predicted catch
  srv_pred = rep(0, Ts)                                            # Predicted survey index
  srv_var = rep(0,Ts)                                              # Survey variance
  fish_age_pred = matrix(0, A1, Tfa)                               # Predicted fishery age comps
  srv_age_pred = matrix(0, A1, Tsa)                                # Predicted survey age comps
  fish_size_pred = matrix(0, L, Tfs)                               # Predicted fishery size comps
  spawn_bio = tot_bio = rep(0, T)                                  # Spawning and total biomass
  N_spr = sb_spr = matrix(1, A, 4)                                 # Numbers and spawning biomass per recruit
  log_rec_cohort = c(rev(init_log_Rt), log_Rt) 
  
  # priors -----------------
  # Priors on key parameters (negative log-likelihood contributions)
  nll_M = (log(M) - log(mean_M))^2 / (2 * cv_M^2)                # Prior on natural mortality
  nll_q = (log(q) - log(mean_q))^2 / (2 * cv_q^2)                # Prior on survey catchability
  nll_sigmaR = (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2) # Prior on recruitment variability
  
  # selectivity ----
  slx_block = matrix(0, A, 4)                                    # Selectivity blocks for fishery
  slx_block[,1] = sel_logistic(1:A, a50C[1], deltaC[1], adj=0)   # Block 1: logistic selectivity
  slx_block[,3] = to_one(sel_double_logistic(1:A, a50C[2], deltaC[2], adj=0)) # Block 3: double logistic selectivity
  slx_block[,2] = to_one((slx_block[,1] + slx_block[,3]) * 0.5)          # Block 2: average of blocks 1 and 3
  slx_block[,4] = to_one(sel_double_logistic(1:A, a50C[3], deltaC[3], adj=0)) # Block 4: double logistic selectivity
  
  
  for(t in 1:T) {
    slx_fish[,t] = slx_block[,fish_block_ind[t]]                 # Assign selectivity by year
  }
  
  slx_srv = sel_logistic(1:A, a50S, deltaS, adj=0)               # Survey selectivity (logistic)
  
  
  # mortality ----
  # Calculate fishing mortality for each year
  Ft = exp(log_mean_F + log_Ft)  # Annual fishing mortality on natural scale
  for(t in 1:T){
    Fat[,t] = Ft[t] * slx_fish[,t]  # Fishing mortality at age and year
    Zat[,t] = Fat[,t] + M           # Total mortality at age and year
  }
  Sat = exp(-Zat)                     # Survivorship at age and year
  
  ## Nat ----
  # Populate numbers-at-age matrix (Nat)
  # First row: recruitment for each year
  # Use correct log_Rt to match ADMB model (init_log_Rt are in reverse order)
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])  # Recruitment in year t
  }
  # First column: initial numbers-at-age for each cohort
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])  # Initial numbers for ages 2 to A-1
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))         # Plus group (oldest age class)
  
  # Remaining columns: survivors from previous year
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]                 # Survivors from previous age and year
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]              # Plus group accumulates survivors
  }
  
  # Calculate recruits and biomasses
  recruits = Nat[1,]                          # Recruitment time series
  spawn_bio = colSums(Nat * wt_mature)        # Spawning biomass by year
  tot_bio = colSums(Nat * waa)                # Total biomass by year
  
  # Adjust spawning biomass in last year for pre-spawning mortality
  spawn_adj = Sat[,T]^(spawn_fract)
  spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
  
  ## catch ----
  # Calculate predicted catch at age and year
  Cat = Fat / Zat * Nat * (1-Sat)
  catch_pred = colSums(Cat * waa) # Predicted catch biomass
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2) # Catch likelihood (sum of squared log differences)
  
  ## survey biomass ----
  isrv = 1
  srv_like = 0.0
  
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      srv_pred[isrv] = sum(Nat[,t] * slx_srv * waa) * q # Predicted survey index
      # Survey likelihood (lognormal, using observed and predicted survey biomass)
      CV = srv_sd[isrv] / srv_obs[isrv]
      log_sd = sqrt(log(1 + CV^2))
      srv_like = srv_like - dnorm(log(srv_obs[isrv]), log(srv_pred[isrv]), log_sd, log = TRUE)
      isrv = isrv + 1
    }
  }
  
  like_srv = srv_like * srv_wt # Weighted survey likelihood
  
  ## fishery age comp ----
  fish_age_lk = 0.0
  offset = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      # Predicted age composition (with ageing error)
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      # Offset for multinomial likelihood
      offset = offset - fish_age_iss[icomp] * 
        sum((fish_age_obs[,icomp] + g) * 
              log(fish_age_obs[,icomp] + g))
      # Multinomial likelihood for age composition
      fish_age_lk = fish_age_lk - sum(fish_age_iss[icomp] * 
                                        (fish_age_obs[,icomp] + g) * 
                                        log(fish_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_age_lk = fish_age_lk - offset
  like_fish_age = fish_age_lk * fish_age_wt # Weighted fishery age comp likelihood
  
  ## survey age comp ----
  srv_age_lk = 0.0
  offset_sa = 0.0
  icomp = 1
  
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      # Predicted survey age composition (with ageing error)
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx_srv) / sum(Nat[,t] * slx_srv) * age_error))
      # Offset for multinomial likelihood
      offset_sa = offset_sa - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_obs[,icomp] + g))
      # Multinomial likelihood for survey age composition
      srv_age_lk = srv_age_lk - srv_age_iss[icomp] * sum((srv_age_obs[,icomp] + g) * log(srv_age_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  srv_age_lk = srv_age_lk - offset_sa
  like_srv_age = srv_age_lk * srv_age_wt # Weighted survey age comp likelihood
  
  ## fishery size comp ----
  icomp = 1
  fish_size_lk = 0.0
  offset_fs = 0.0
  
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      # Predicted size composition (with size-at-age array)
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * saa_array[,,fish_saa_ind[t]]))
      # Offset for multinomial likelihood
      offset_fs = offset_fs - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_obs[,icomp] + g))
      # Multinomial likelihood for size composition
      fish_size_lk = fish_size_lk - fish_size_iss[icomp] * sum((fish_size_obs[,icomp] + g) * log(fish_size_pred[,icomp] + g))
      icomp = icomp + 1
    }
  }
  fish_size_lk = fish_size_lk - offset_fs
  like_fish_size = fish_size_lk * fish_size_wt # Weighted fishery size comp likelihood
  
  # SPR ------------------------
  # Prepare recruitment data frame for reference point calculations
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # Filter years for recruitment estimation (exclude first and last ages)
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  n_rec = nrow(df)
  yrs_rec = df$year
  pred_rec = mean(df$pred_rec) # Mean predicted recruitment
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1)) # Recruitment SD
  
  # Calculate numbers per recruit for reference points (F50, F40, F35)
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx_fish[a-1,T]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx_fish[a-1,T]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx_fish[a-1,T]))
  }
  # Plus group for per-recruit calculations
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F50 * slx_fish[A,T])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F40 * slx_fish[A,T])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx_fish[A-1,T])) / 
    (1 - exp(-(M + F35 * slx_fish[A,T])))
  
  # Calculate spawning biomass per recruit for reference points
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx_fish[a,T]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx_fish[a,T]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx_fish[a,T]))
  }
  
  # Calculate reference point spawning biomasses
  SB0 = sum(sb_spr[,1])    # Unfished spawning biomass per recruit
  SBF50 = sum(sb_spr[,2])  # Spawning biomass per recruit at F50
  SBF40 = sum(sb_spr[,3])  # Spawning biomass per recruit at F40
  SBF35 = sum(sb_spr[,4])  # Spawning biomass per recruit at F35
  
  # SPR penalties to enforce reference point constraints
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  # Scale reference points by mean recruitment
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # likelihood/penalties --------------------
  like_rec = (sum(c(log_Rt, init_log_Rt)^2) / (2 * sigmaR^2) + length(c(log_Rt, init_log_Rt)) * log(sigmaR) ) * wt_rec_var
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  cohort_penalty = 0
  n_cohort = length(log_rec_cohort)
  for (i in 2:(n_cohort-1)) {
    d2 = log_rec_cohort[i-1] - 2 * log_rec_cohort[i] + log_rec_cohort[i+1]
    cohort_penalty = cohort_penalty + d2^2
  }
  
  # nll ----
  nll = ssqcatch
  nll = nll + like_srv
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec 
  nll = nll + f_regularity
  nll = nll + cohort_penalty
  nll = nll + nll_M
  nll = nll + nll_q
  nll = nll + nll_sigmaR
  nll = nll + sprpen    
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::ADREPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::ADREPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::ADREPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(maa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx_srv)
  RTMB::REPORT(slx_fish)
  RTMB::REPORT(slx_block)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(recruits)
  RTMB::ADREPORT(srv_pred)
  RTMB::ADREPORT(tot_bio)
  RTMB::ADREPORT(spawn_bio)
  RTMB::ADREPORT(recruits)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(n_rec)
  RTMB::REPORT(yrs_rec)
  RTMB::REPORT(stdev_rec)
  
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_sigmaR)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

