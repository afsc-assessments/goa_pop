# September 2025 plan team examinations

# load ----
library(RTMButils)
# devtools::unload('RTMButils')
library(tidyverse)
# library(Matrix)
# library(here)
# library(scico)
library(patchwork)
theme_set(afscassess::theme_report())

# source(here::here(2025, 'r', "utils.R"))
source(here::here(2025, 'r', "models.R"))

# data ----
# setup data and pars for all models 
# 2023 data & pars from ADMB model
data = readRDS(here::here(2025, 'rtmb_bridge', "data.rds"))
# pars = readRDS(here::here(2025, 'rtmb_bridge', "pars.rds"))

# original ADMB starting pars
pars = list(log_M = log(0.0614),
            log_a50C = log(c(6, 2.5, 2.5)),
            deltaC = c(1.5,4.5, 4.5),
            log_a50S = log(7.3),
            deltaS = 3.8,
            log_q = log(1.15),
            log_mean_R = 3.0,
            init_log_Rt = rep(0, 26),
            log_Rt = rep(0, sum(data$catch_ind)),
            log_mean_F = 0,
            log_Ft = rep(0, sum(data$catch_ind)),
            log_F35 = 0,
            log_F40 = 0,
            log_F50 = 0,
            sigmaR = 1.7)

# og management model -----
m25 = run_model(base, data, pars) # base model run

# improve gradient
cmb = function(f, d) function(p) f(p, d)
obj =  RTMB::MakeADFun(cmb(base, data),
                       pars)

optim = nlminb(start = obj$par,
               objective = obj$fn,
               gradient = obj$gr,
               control = list(iter.max=1e+07,
                              eval.max=1e+07))

newton_loops = 3

try_improve <- tryCatch(expr = for (i in 1:newton_loops) {
  g = as.numeric(obj$gr(optim$par))
  h = optimHess(optim$par, fn = obj$fn, gr = obj$gr)
  optim$par = optim$par - solve(h, g)
  optim$objective = obj$fn(optim$par)
}, error = function(e) {
  e
}, warning = function(w) {
  w
})

rpt = obj$report(obj$env$last.par.best)
proj = proj_bio(rpt) # function to project the next 2 years
sd = sdreport(obj)
m25_test = list(obj = obj, fit = optim, rpt = rpt, proj = proj, sd = sd)


fit_check(m25)
fit_check(m25_test)


