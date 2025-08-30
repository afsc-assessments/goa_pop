# bridge GOA Pacific ocean perch from ADMB to RTMB
# ben.williams@noaa.gov
# 2025-05

# load ----
library(RTMButils)
library(tidyverse)
# library(Matrix)
# library(here)
# library(scico)
theme_set(afscassess::theme_report())

# globals ---- 
year = 2025
ages = 2:25
years = 1961:2023
length_bins = 16:45

# source(here::here(year, 'r', "utils.R"))
source(here::here(year, 'r', "models.R"))

# bridge data and results from 2023 ADMB model
REP = readLines(here::here(year, "base", "model_20_1.rep"))
PAR = readLines(here::here(year, "base", "model_20_1.par"))
STD = read_table(here::here(year, "base", "model_20_1.std"))
DAT = readLines(here::here(year, "base", "goa_pop_2023.dat"))
CTL = readLines(here::here(year, "base", "goa_pop_2023.ctl"))

# data ----
# build data inputs directly from the data used in the ADMB model

# maa = as.numeric(stringr::str_split(REP[grep('Maturity', REP)], " ")[[1]][3:30])
# waa = as.numeric(stringr::str_split(REP[grep('Weight', REP)], " ")[[1]][3:30])
# wt_mature = waa * maa / 2
# fish_block_ind = case_when(
#   years <= 1976 ~ 1,
#   years <= 1995 ~ 2,
#   years <= 2006 ~ 3,
#   TRUE ~ 4
# )
# 
# a50_vals = c(as.numeric(PAR[grep('\\ba502:', PAR)+1]), as.numeric(PAR[grep('\\ba50:', PAR)+1]), as.numeric(PAR[grep('\\ba503:', PAR)+1]))
# delta_vals = c(as.numeric(PAR[grep('\\bdelta2:', PAR)+1]), as.numeric(PAR[grep('\\bdelta:', PAR)+1]),as.numeric(PAR[grep('\\bdelta3:', PAR)+1]))
# 
# catch_obs = c(as.double(stringr::str_split(DAT[grep("Fishery catch", DAT)+3], " ")[[1]]))
# 
# srv_yrs = as.numeric(stringr::str_split(DAT[grep('Trawl survey years', DAT)+1], " ")[[1]])
# srv_ind = ifelse(years %in% srv_yrs, 1, 0)
# srv_obs = as.numeric(stringr::str_split(DAT[grep('Trawl survey years', DAT)+3], " ")[[1]])
# srv_sd = as.numeric(stringr::str_split(DAT[grep('Trawl survey years', DAT)+5], " ")[[1]])
# 
# fish_age_yrs = as.numeric(stringr::str_split(DAT[grep('Fishery Age Composition', DAT)+5], " ")[[1]])
# as.data.frame(t(mapply(rbind, stringr::str_split(REP[grep('Obs_P_fish_age', REP)+1:length(fish_age_yrs)], " ")))) %>%
#   dplyr::select(-V2) -> fish_ac
# names(fish_ac) <- c('year', paste0('age-', ages))
# fish_ac %>%
#   mutate(across(c(year, starts_with('age')), as.numeric)) -> fish_ac
# 
# fish_age_ind = ifelse(years %in% fish_age_yrs, 1, 0)
# fish_age_iss = sqrt(as.numeric(stringr::str_split(DAT[grep('nsamples_fish_age', DAT)+1], " ")[[1]]))
# fish_age_obs = fish_ac %>%
#   select(starts_with('age')) %>%
#   as.matrix() %>%
#   t() %>%
#   unname()
# 
# srv_age_yrs = as.numeric(stringr::str_split(DAT[grep('Trawl Survey Age Composition', DAT)+5], " ")[[1]])
# as.data.frame(t(mapply(rbind, stringr::str_split(REP[grep('Obs_P_srv1_age', REP)+1:length(srv_age_yrs)], " ")))) %>%
#   dplyr::select(-V2) -> srv_ac
# names(srv_ac) <- c('year', paste0('age-', ages))
# srv_ac %>%
#   mutate(across(c(year, starts_with('age')), as.numeric)) -> srv_ac
# 
# srv_age_ind = ifelse(years %in% srv_age_yrs, 1, 0)
# srv_age_iss = sqrt(as.numeric(stringr::str_split(DAT[grep('nsamples_srv1_age', DAT)+1], " ")[[1]]))
# srv_age_obs = srv_ac %>%
#   select(starts_with('age')) %>%
#   as.matrix() %>%
#   t() %>%
#   unname()
# 
# fish_size_yrs = as.numeric(stringr::str_split(DAT[grep('Fishery Size Composition', DAT)+5], " ")[[1]])
# as.data.frame(t(mapply(rbind, stringr::str_split(REP[grep('Obs_P_fish_size', REP)+1:length(fish_size_yrs)], " ")))) %>%
#   dplyr::select(-V2) -> fish_sc
# names(fish_sc) <- c('year', paste0('l-', length_bins))
# fish_sc %>%
#   mutate(across(c(year, starts_with('l')), as.numeric)) -> fish_sc
# 
# fish_size_ind = ifelse(years %in% fish_size_yrs, 1, 0)
# fish_size_iss = as.numeric(stringr::str_split(DAT[grep('Number of samples:', DAT)[3]+1], " ")[[1]])
# fish_size_obs = fish_sc %>%
#   select(starts_with('l')) %>%
#   as.matrix() %>%
#   t() %>%
#   unname()
# 
# as.data.frame(t(mapply(rbind, stringr::str_split(DAT[grep('Size-age transition matrix', DAT)+3:(2+length(maa))], " ")))) %>%
#   mutate(across(everything(), as.numeric)) %>%
#   as.matrix() -> saa
# 
# as.data.frame(t(mapply(rbind, stringr::str_split(DAT[grep('Size-age matrix using mean estimation', DAT)+1:(length(maa))], " ")))) %>%
#   mutate(across(everything(), as.numeric)) %>%
#   as.matrix() -> saa2
# 
# fish_saa_ind = ifelse(years >= 1990, 2, 1)
# saa_array <- array(NA, dim = c(28, 30, 2))
# 
# # Assign the two matrices into the 3D array
# saa_array[,,1] <- saa
# saa_array[,,2] <- saa2
# 
# as.data.frame(t(mapply(rbind, stringr::str_split(DAT[grep('Ageing error matrix', DAT)+2:(1+length(maa))], " ")))) %>%
#   mutate(across(everything(), as.numeric)) %>%
#   as.matrix() -> ae
# 
# yield_ratio = as.numeric(stringr::str_split(CTL[grep("yield", CTL)], " ")[[1]][1])
# 
# data = list(ages = as.integer(ages),
#             years = as.integer(years),
#             length_bins = as.integer(length_bins),
#             waa = waa,
#             maa = maa,
#             wt_mature = wt_mature,
#             spawn_mo = 5,
#             fish_block_ind = as.integer(fish_block_ind),
#             catch_obs = catch_obs,
#             catch_ind = as.integer(c(rep(1,18), rep(1, 45))),
#             catch_wt = rep(50, length(years)),
#             # fish_block_ind = rep(1, length(years)),
#             # fish_sel_type = rep(0, length(years)),
#             srv_yrs = as.integer(srv_yrs),
#             srv_ind = as.integer(srv_ind),
#             srv_obs = srv_obs,
#             srv_sd = srv_sd,
#             srv_wt = 1.0,
#             fish_age_yrs = as.integer(fish_age_yrs),
#             fish_age_ind = as.integer(fish_age_ind),
#             fish_age_iss = fish_age_iss,
#             fish_age_obs = fish_age_obs,
#             fish_age_wt = 1.0,
#             srv_age_yrs = as.integer(srv_age_yrs),
#             srv_age_ind = as.integer(srv_age_ind),
#             srv_age_iss = srv_age_iss,
#             srv_age_obs = srv_age_obs,
#             srv_age_wt = 1.0,
#             fish_size_yrs = as.integer(fish_size_yrs),
#             fish_size_ind = as.integer(fish_size_ind),
#             fish_size_iss = fish_size_iss,
#             fish_size_obs = fish_size_obs,
#             fish_size_wt = 1.0,
#             age_error = unname(ae),
#             size_age = unname(saa),
#             fish_saa_ind = as.integer(fish_saa_ind),
#             saa_array = saa_array,
#             wt_fmort_reg = 0.1,
#             wt_rec_var = 1,
#             mean_M = 0.0614,
#             cv_M = 0.1,
#             mean_q = 1.15,
#             cv_q = 0.447213595,
#             mean_sigmaR = 1.7,
#             cv_sigmaR = 0.2,
#             yield_ratio = yield_ratio
# )
# str(data)
# saveRDS(data, here::here(year, 'rtmb_bridge', "data.rds"))
# 
# # parameters ----
# build parameter inputs directly from the ADMB model 

# A = nrow(ae)
# pars = list(log_M = as.numeric(PAR[grep('logm', PAR)+1]),
#             log_a50C = log(a50_vals),
#             deltaC = delta_vals,
#             log_a50S = log(as.numeric(PAR[grep('\\ba50_srv1:', PAR)+1])),
#             deltaS = as.numeric(PAR[grep('\\bdelta_srv1:', PAR)+1]),
#             log_q = as.numeric(PAR[grep('\\blog_q_srv1:', PAR)+1]),
#             log_mean_R = as.numeric(PAR[grep('\\blog_mean_rec:', PAR)+1]),
#             init_log_Rt = rev(as.numeric(stringr::str_split(PAR[grep('\\log_rec_dev:', PAR)+1], " ")[[1]])[2:(A-1)]),
#             log_Rt = as.numeric(stringr::str_split(PAR[grep('\\log_rec_dev:', PAR)+1], " ")[[1]])[-(1:(A-1))],
#             log_mean_F = as.numeric(PAR[grep('\\blog_avg_F:', PAR)+1]),
#             log_Ft = as.numeric(stringr::str_split(PAR[grep('\\blog_F_devs:', PAR)+1], " ")[[1]])[2:64],
#             log_F35 = log(as.numeric(PAR[grep('\\bmF35:', PAR)+1])),
#             log_F40 = log(as.numeric(PAR[grep('\\bmF40:', PAR)+1])),
#             log_F50 = log(as.numeric(PAR[grep('\\bmF50:', PAR)+1])),
#             sigmaR = as.numeric(PAR[grep('\\bsigr:', PAR)+1])
# )
# 
# str(pars)
# saveRDS(pars, here::here(year, 'rtmb_bridge', "pars.rds"))

data = readRDS(here::here(year, 'rtmb_bridge', "data.rds"))
pars = readRDS(here::here(year, 'rtmb_bridge', "pars.rds"))

# original ADMB starting pars
# pars = list(log_M = log(0.0614),
#             log_a50C = log(c(6, 2.5, 2.5)),
#             deltaC = c(1.5,4.5, 4.5),
#             log_a50S = log(7.3),
#             deltaS = 3.8,
#             log_q = log(1.15),
#             log_mean_R = 3.0,
#             init_log_Rt = rep(0, 26),
#             log_Rt = rep(0, sum(data$catch_ind)),
#             log_mean_F = 0,
#             log_Ft = rep(0, sum(data$catch_ind)),
#             log_F35 = 0,
#             log_F40 = 0,
#             log_F50 = 0,
#             sigmaR = 1.7)

# limits ----
# not currently used 
# lower = list(
#   pars$log_M = log(0.03),           # Minimum natural mortality
#   log_a50C = rep(log(0.5), 3), # Minimum age at 50% selectivity
#   deltaC = rep(0.1, 3),        # Minimum selectivity slope
#   log_a50S = log(0.5),         # Minimum survey selectivity age
#   deltaS = 0.1,                # Minimum survey selectivity slope
#   log_q = log(0.01),           # Minimum catchability
#   log_mean_R = 2,              # Minimum mean recruitment
#   init_log_Rt = rep(-10, 26),  # Minimum initial recruitment devs
#   log_Rt = rep(-10, sum(data$catch_ind)), # Minimum recruitment devs
#   log_mean_F = -5,             # Minimum mean fishing mortality
#   log_Ft = rep(-10, sum(data$catch_ind)), # Minimum F devs
#   log_F35 = -5,               # Minimum F35
#   log_F40 = -5,               # Minimum F40
#   log_F50 = -5,               # Minimum F50
#   sigmaR = 0.1                # Minimum recruitment SD
# )

# upper = list(
#   pars$log_M = log(0.27),           # Maximum natural mortality
#   log_a50C = rep(log(20), 3), # Maximum age at 50% selectivity
#   deltaC = rep(20, 3),        # Maximum selectivity slope
#   log_a50S = log(20),         # Maximum survey selectivity age
#   deltaS = 20,                # Maximum survey selectivity slope
#   log_q = log(10),            # Maximum catchability
#   log_mean_R = 10,            # Maximum mean recruitment
#   init_log_Rt = rep(10, 26),  # Maximum initial recruitment devs
#   log_Rt = rep(10, sum(data$catch_ind)), # Maximum recruitment devs
#   log_mean_F = 1,             # Maximum mean fishing mortality
#   log_Ft = rep(5, sum(data$catch_ind)),  # Maximum F devs
#   log_F35 = 1,                # Maximum F35
#   log_F40 = 1,                # Maximum F40
#   log_F50 = 1,                # Maximum F50
#   sigmaR = 1.5                  # Maximum recruitment SD
# )

# model runs ----
# match the base model w/o optimizing
obj <- RTMB::MakeADFun(cmb(base, data),
                       pars)
rpt <- obj$report(obj$env$last.par.best)
proj_bio(rpt)

# match the base model w/optimizing
fit <- nlminb(start = obj$par,
              objective = obj$fn,
              gradient = obj$gr,
              control = list(iter.max=100000,
                             eval.max=20000))
rpt <- obj$report(obj$env$last.par.best)
saveRDS(obj, here::here(2025, 'rtmb_bridge', 'results', 'obj.RDS'))
saveRDS(fit, here::here(2025, 'rtmb_bridge', 'results', 'fit.RDS'))
saveRDS(rpt, here::here(2025, 'rtmb_bridge', 'results', 'report.RDS'))

prj = proj_bio(rpt)[1,]
sd_fit = sdreport(obj, getJointPrecision = TRUE)

# compare ----
## SD
as.data.frame(summary(sd_fit)) %>% 
  tibble::rownames_to_column("name") %>% 
  dplyr::rename(value = Estimate, std.dev = `Std. Error`)-> sds

bind_rows(filter(STD, name=='log_F_devs') %>% 
            mutate(name = "log_Ft",
                   year = rpt$years,
                   ll = value - std.dev,
                   ul = value + std.dev,
                   id = 'ADMB')) %>% 
  bind_rows(sds %>%
              filter(stringr::str_detect(name, "Ft")) %>%
              mutate(name = gsub("\\.\\d+", "", name),
                     year = rpt$years,
                     ll = value -std.dev,
                     ul = value + std.dev,
                     id = "RTMB")) %>% 
  bind_rows(filter(STD, name=='pred_rec') %>% 
  mutate(name = "recruits",
         year = rpt$years,
         ll = value - std.dev,
         ul = value + std.dev,
         id = 'ADMB')) %>% 
  bind_rows(sds %>%
              filter(stringr::str_detect(name, "recruits")) %>%
              mutate(name = gsub("\\.\\d+", "", name),
                     year = rpt$years,
                     ll = value -std.dev,
                     ul = value + std.dev,
                     id = "RTMB")) %>% 
bind_rows(filter(STD, name=='spawn_biom') %>% 
  mutate(name = "spawn_bio",
         year = rpt$years,
         ll = value - std.dev,
         ul = value + std.dev,
         id = 'ADMB')) %>% 
  bind_rows(sds %>%
              filter(stringr::str_detect(name, "spawn_bio")) %>%
              mutate(name = gsub("\\.\\d+", "", name),
                     year = rpt$years,
                     ll = value -std.dev,
                     ul = value + std.dev,
                     id = "RTMB")) %>% 
  bind_rows(filter(STD, name=='tot_biom') %>% 
  mutate(name = "tot_bio",
         year = rpt$years,
         ll = value - std.dev,
         ul = value + std.dev,
         id = 'ADMB')) %>% 
  bind_rows(sds %>%
              filter(stringr::str_detect(name, "tot_bio")) %>%
              mutate(name = gsub("\\.\\d+", "", name),
                     year = rpt$years,
                     ll = value -std.dev,
                     ul = value + std.dev,
                     id = "RTMB")) %>% 
  select(name, value, year, ll, ul, Model = id) %>% 
  mutate(name = case_when(name=='spawn_bio' ~ "Spawning Biomass",
                          name=='tot_bio' ~ "Total Biomass",
                          name=='log_Ft' ~ "Log Fishing Mortality Devs",
                          name=='recruits' ~ "Recruitment")) %>% 
  ggplot(aes(year, value, color = Model, fill = Model, linetype = Model)) + 
  geom_ribbon(aes(ymin=ll, ymax=ul), alpha = 0.2, color = NA) +
  geom_line() +
  facet_wrap(~name, scales = "free") +
  scico::scale_color_scico_d(end = 0.7) +
  scico::scale_fill_scico_d(end = 0.7) +
  xlab("Year") +
  ylab("")
ggsave(here::here(2025, "sep_pt", "figs", "biomass.png"))


filter(STD, name=='log_mean_rec') %>% 
  bind_rows(filter(sds, name=="log_mean_R"))

filter(STD, name=='sigr') %>% 
  bind_rows(filter(sds, name=="sigmaR"))

filter(STD, name=='nattymort') %>% 
  bind_rows(filter(sds, name=="log_M")) %>% 
  mutate(value = case_when(name=='log_M'~ exp(value), 
                            TRUE ~ value),
         std.dev = case_when(name=='log_M'~ value * std.dev, # delta method
                           TRUE ~ std.dev))
filter(STD, name=='q_srv1') %>% 
  bind_rows(filter(sds, name=="log_q")) %>% 
  mutate(value = case_when(name=='log_q'~ exp(value), 
                           TRUE ~ value),
         std.dev = case_when(name=='log_q'~ value * std.dev, # delta method
                             TRUE ~ std.dev))
filter(STD, name=='q_srv1') %>% 
  bind_rows(filter(sds, name=="log_q")) %>% 
  mutate(value = case_when(name=='log_q'~ exp(value), 
                           TRUE ~ value),
         std.dev = case_when(name=='log_q'~ value * std.dev, # delta method
                             TRUE ~ std.dev))

# get ADMB likelihood values
ssqcatch = round(as.numeric(stringr::str_split(REP[grep("SSQ Catch Likelihood", REP)], " ")[[1]][2]), 4)
srv_like = round(as.numeric(stringr::str_split(REP[grep("Bottom Trawl Survey Likelihood", REP)], " ")[[1]][2]), 4)
fa_like = round(as.numeric(stringr::str_split(REP[grep("Fishery Age Composition Likelihood", REP)], " ")[[1]][2]), 4)
sa_like = round(as.numeric(stringr::str_split(REP[grep("Bottom Trawl Survey Age Composition Likelihood", REP)], " ")[[1]][2]), 4)
ss_like = round(as.numeric(stringr::str_split(REP[grep("Fishery Size Composition Likelihood", REP)], " ")[[1]][2]), 4)
recdev = round(as.numeric(stringr::str_split(REP[grep("Recruitment Deviations Likelihood", REP)], " ")[[1]][2]), 4)
fmort = round(as.numeric(stringr::str_split(REP[grep("Fishing Mortality Deviations Penalty", REP)], " ")[[1]][2]), 4)
pM = round(as.numeric(stringr::str_split(REP[grep("Priors M", REP)], " ")[[1]][2]), 4)
pq = round(as.numeric(stringr::str_split(REP[grep("Priors q", REP)], " ")[[1]][2]), 4)
pS = round(as.numeric(stringr::str_split(REP[grep("Priors SigmaR", REP)], " ")[[1]][2]), 4)
spr = round(as.numeric(stringr::str_split(REP[grep("Spawner-recruit penalty", REP)], " ")[[1]][2]), 4)
datalike = round(as.numeric(stringr::str_split(REP[grep("Data Likelihood", REP)], " ")[[1]][2]), 4)

# likelihood table ----
data.frame(Likelihood = c("Catch", "Survey", "Fish age", "Survey age", "Fish size", "Recruitment", "F regularity", "SPR penalty", "M prior", "q prior", "Sigma R prior", "Sub total"),
           ADMB = c(ssqcatch, srv_like, fa_like, sa_like, ss_like, recdev, fmort, spr, pM, pq, pS, 
                    sum(c(ssqcatch, srv_like, fa_like, sa_like, ss_like, recdev, fmort, spr, pM, pq, pS))),
           RTMB = round(c(rpt$ssqcatch, rpt$like_srv, rpt$like_fish_age, rpt$like_srv_age, 
                          rpt$like_fish_size, rpt$like_rec, rpt$f_regularity,
                          rpt$sprpen, rpt$nll_M, rpt$nll_q, rpt$nll_sigmaR, 
                          sum(c(rpt$ssqcatch, rpt$like_srv, rpt$like_fish_age, rpt$like_srv_age, 
                                rpt$like_fish_size, rpt$like_rec, rpt$f_regularity,
                                rpt$sprpen, rpt$nll_M, rpt$nll_q, rpt$nll_sigmaR))), 4)) %>% 
  mutate(Difference =  ADMB - RTMB) %>% 
  vroom::vroom_write(here::here(2025, "sep_pt", "tables", "like_tbl.csv"), delim=",")


# key parameters ----
tot = as.numeric(REP[grep("TotalBiomass for 2024", REP)+1])
sp = as.numeric(REP[grep("Female_Spawning Biomass for 2024", REP)+1])
ofl = as.numeric(REP[grep("\\bOFL for 2024", REP)+1])
fofl = as.numeric(REP[grep("F_OFL for 2024", REP)+1])
abc = as.numeric(REP[grep("\\bABC for 2024", REP)+1])
b40 = as.numeric(REP[grep("B_40", REP)+1])
fabc = as.numeric(REP[grep("F_ABC for 2024", REP)+1])
prj = proj_bio(rpt)[1,]

# data = readRDS(here::here(year, 'rtmb_bridge', "data.rds"))
pars = readRDS(here::here(year, 'rtmb_bridge', "pars.rds"))

data.frame(Item = c("M", "q", "Log mean recruitment", "Log mean F", 'a50_1', 'delta_1', 'a50_3', 'delta_3', 'a50_4', 'delta_4', 'a50_survey', 'delta_survey',
                    "2024 Total biomass", "2024 Spawning biomass", "B40", "2024 OFL", "2024 F OFL", " 2024 ABC", "2024 F ABC"),
           ADMB = round(c(exp(pars$log_M), exp(pars$log_q), pars$log_mean_R, pars$log_mean_F, exp(pars$log_a50C)[1], 
                          pars$deltaC[1], exp(pars$log_a50C)[2], pars$deltaC[2], exp(pars$log_a50C)[3], pars$deltaC[3], 
                          exp(pars$log_a50S), pars$deltaS, tot, sp, b40, ofl, fofl, abc, fabc), 4),
           RTMB = round(c(rpt$M, rpt$q, rpt$log_mean_R, rpt$log_mean_F, rpt$a50C[1], rpt$deltaC[1], 
                          rpt$a50C[2], rpt$deltaC[2], rpt$a50C[3], rpt$deltaC[3], rpt$a50S, rpt$deltaS,
                          prj$tot_bio, prj$spawn_bio, rpt$B40, prj$catch_ofl, prj$F35, prj$catch_abc, prj$F40),4)) %>% 
  mutate(Difference = ADMB - RTMB) %>% 
  # flextable::flextable() %>% 
  # flextable::colformat_double(
  #   i = c(13:15,17),  j = 2:4,
  #   big.mark = ",", 
  #   digits = 2, 
  #   na_str = "N/A"
  # ) %>% 
  vroom::vroom_write(here::here(2025, "sep_pt", "tables", "par_tbl.csv"), delim=",")

# data likelihood
rpt$ssqcatch + rpt$like_srv + rpt$like_srv_age + rpt$like_fish_age + rpt$like_fish_size # data objective function


# basic comparison plots ----
# admb slx
data.frame(age = 2:29,
           X1 = as.numeric(stringr::str_split(REP[grep("Selectivity", REP)[[1]]], " ")[[1]][-1]),
           X2 = as.numeric(stringr::str_split(REP[grep("Selectivity", REP)[[2]]], " ")[[1]][-1]),
           X3 = as.numeric(stringr::str_split(REP[grep("Selectivity", REP)[[3]]], " ")[[1]][-1]),
           X4 = as.numeric(stringr::str_split(REP[grep("Selectivity", REP)[[4]]], " ")[[1]][-1]),
           srv_slx = as.numeric(stringr::str_split(REP[grep("Bottom_Trawl_Survey_Selectivity", REP)], " ")[[1]][-c(1:2)]),
           id = "ADMB") %>% 
  bind_rows(
    data.frame(age = 2:29,
           rpt$slx_block,
           srv_slx = rpt$slx_srv,
           id = "RTMB")) %>% 
  rename('Fishery 1' = X1, 'Fishery 2' = X2, 'Fishery 3' = X3, 'Fishery 4' = X4, "Survey" = srv_slx, Model = id) %>% 
  tidyr::pivot_longer(-c(age, Model)) %>% 
  ggplot(aes(age, value, color = Model, linetype = Model)) + 
  geom_line() + 
  facet_wrap(~name, scales = "free_y") +
  xlab("Age") +
  ylab("Selectivity") +
  scico::scale_colour_scico_d(end=0.7) +
  theme(legend.position = c(0.8, 0.25))

ggsave(here::here(2025, "sep_pt", "figs", "slx.png"))

