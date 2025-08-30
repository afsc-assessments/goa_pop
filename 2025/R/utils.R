# note that the 2nd time period in the selectivity blocks is not scaled to one
# so have this add on to account for that 
to_one <- function(x) {
  x / max(x)
}


#' standard MESA logistic selectivity
#'
#' @param age age to examine
#' @param a50 inflection point
#' @param delta steepness of slope
#' @param adj based upon index that starts at 1 - unless adjusted to start at recruitment age
#' doesn't really do anything other than scale
sel_logistic <- function(age, a50, delta, adj) {
  sel = 1 / (1 + exp(-log(19) * ((age + adj) - a50) / delta))
  sel #= sel / max(sel)
}

#' standard MESA double logistic selectivity
#'
#' @param age age to examine
#' @param a50 inflection point
#' @param delta steepness of slope
#' @param adj based upon index that starts at 1 - unless adjusted to start at recruitment age
#' doesn't really do anything other than scale
sel_double_logistic <- function(age, a50, delta, adj) {
  expa50 = exp(a50)  # log scale for dome
  denom = 0.5 * (sqrt(expa50^2 + 4 * delta^2) - expa50)
  sel = (((age + adj) / expa50)^(expa50 / denom)) * exp((expa50 - (age + adj)) / denom)
  sel #= sel / max(sel)
}

sel_double_normal <- function(age, a50, delta, adj = 0) {
  sel <- exp(-((age + adj - a50)^2) / (2 * delta^2))
  sel
}

# Helper function for Dirichlet-multinomial likelihood
dirmultinom <- function(x, prob, theta, log = FALSE) {
  # x: observed proportions
  # prob: predicted proportions  
  # theta: effective sample size
  n = sum(x)
  alpha = prob * theta
  ll = lgamma(n + 1) + lgamma(theta) - lgamma(n + theta) +
    sum(lgamma(x + alpha) - lgamma(alpha) - lgamma(x + 1))
  if(log) return(ll) else return(exp(ll))
}

#' Projected biomass for 2 years
#'
#' @param report RTMB report object
#' @param Tproj number of years to project, default: 2
proj_bio <- function(report, Tproj=2) {
  
  # values
  F40 = report$F40
  F35 = report$F35
  B40 = report$B40
  Nat = report$Nat
  Sat = report$Sat
  ages = report$ages
  years = report$years
  waa = report$waa
  wt_mature = report$wt_mature
  spawn_frac = report$spawn_fract
  yield_ratio = report$yield_ratio
  M = report$M
  pred_rec = report$pred_rec
  stdev_rec = report$stdev_rec 
  A = nrow(Nat) # number of ages
  T = ncol(Nat) # number of years
  slx = report$slx_fish[,T]
  
  # storage
  N = Cat = Cat_ofl= Zabc = Zofl = S = matrix(0, A, Tproj)
  tot_bio = spawn_bio = F40_proj = F35_proj= rep(0, Tproj)
  # setup
  F40_proj[1] = F40
  F35_proj[1] = F35
  
  # total F
  Fabc_tot = slx * F40_proj[1]
  Fofl_tot = slx * F35_proj[1]
  
  # first projection year
  N[1,] = pred_rec
  for(a in 1:(A-1)) {
    N[a+1,1] = Nat[a,T] * Sat[a,T]
  }
  N[A,1] = Nat[A-1,T] * Sat[A-1,T] + Nat[A,T] * Sat[A,T]
  spawn_bio[1] = sum(N[,1] * exp(-yield_ratio * Fabc_tot - M)^spawn_frac * wt_mature)
  
  for(t in 1:Tproj) {
    # tier check
    if((spawn_bio[t] / B40) > 1) {
      F40_proj[t] = F40
      F35_proj[t] = F35
    } else {
      F40_proj[t] = F40_proj[t] * (spawn_bio[t] / B40 - 0.05) / 0.95
      F35_proj[t] = F35_proj[t] * (spawn_bio[t] / B40 - 0.05) / 0.95
    }
    # update Fs
    Fabc_tot = slx * F40_proj[t]
    Fofl_tot = slx * F35_proj[t]
    Z = Fabc_tot + M
    Zofl = Fofl_tot + M
    S = exp(-Z)
    
    # catch
    Cat[,t] = yield_ratio * N[,t] * Fabc_tot / Z * (1 - S)
    Cat_ofl[,t] = yield_ratio * N[,t] * Fofl_tot / Zofl * (1 - exp(-Zofl))
    
    if(t<Tproj) {
      for(a in 1:(A-1)){
        N[a+1,t+1] = N[a,t] * exp(-yield_ratio * Fabc_tot[a] - M)
      }
      N[A,t+1] = N[A-1,t] * exp(-yield_ratio * Fabc_tot[A-1] - M) +
        N[A,t] * exp(-yield_ratio * Fabc_tot[A] - M)
      
      tot_bio[t+1] = sum(N[,t+1] * waa)
      spawn_bio[t+1] = sum(N[,t+1] * exp(-yield_ratio * Fabc_tot - M)^spawn_frac * wt_mature)
    }
  }
  catch = colSums(Cat * waa / yield_ratio)
  catch_ofl = colSums(Cat_ofl * waa / yield_ratio)
  tot_bio = colSums(N * waa)
  
  data.frame(year = max(years)+1:Tproj,
             spawn_bio = spawn_bio,
             tot_bio = tot_bio,
             catch_abc = catch,
             catch_ofl = catch_ofl,
             F40 = F40_proj,
             F35 = F35_proj)
}

#' pull likelihoods from report 
#'
#' @param report model report
#' @param model name for column
#' @param addl any additional parameters to pull
#' @param exclude any parameters to exclude
#'
get_likes <- function(report, model = "Model", addl=NULL, exclude=NULL) {
  
  items = paste(c("like", "nll", "spr", "regularity", "ssqcatch", addl), collapse = "|")
  selected = report[grep(items, names(report))]
  if (!is.null(exclude)) {
    exclude = paste(exclude, collapse = "|")
    selected = selected[!grepl(exclude, names(selected))]
  }
  
  df <- data.frame(
    item = names(selected),
    value = round(unlist(selected),4),
    row.names = NULL
  )
  
  names(df)[names(df) == "value"] <- model
  df
}

#' pull parameters from report and projection
#'
#' @param report model report
#' @param proj model projection
#' @param model name for column
#' @param addl any additional parameters to pull
#' @param exclude any parameters to exclude
#'
get_pars <- function(report, proj, model = "Model", addl=NULL, exclude=NULL) {
  prj = proj[1,]
  items = paste0("^", c("M", "q", "log_mean_R", "log_mean_F", "a50C", "deltaC", "a50S", "deltaS", addl, "$"), collapse = "|")
  selected = report[grep(items, names(report))]
  
  items2 = data.frame(item = c("tot_bio", "spawn_bio", "catch_ofl", "F35", "catch_abc", "F40"),
                      value = round(c(prj$tot_bio, prj$spawn_bio, prj$catch_ofl, prj$F35, prj$catch_abc, prj$F40), 4))
  
  # Flatten values and preserve indices
  flat = lapply(seq_along(selected), function(i) {
    value = selected[[i]]
    base_name = names(selected)[i]
    
    # If value is a vector, give indexed names
    if (length(value) > 1) {
      data.frame(
        item = paste0(base_name, seq_along(value)),
        value = value
      )
    } else {
      data.frame(
        item = base_name,
        value = round(value, 4)
      )
    }
  })
  
  df <- do.call(rbind, flat)
  
  df = dplyr::bind_rows(df, items2)
  
  if (!is.null(exclude)) {
    exclude = paste(exclude, collapse = "|")
    df = df[!grepl(exclude, df$item),]
  }
  
  
  names(df)[names(df) == "value"] <- model
  df
  
}


#' Run RTMB model
#'
#' @param model the model function
#' @param data named data list
#' @param pars named parameter list
#' @param map named mapping list, default: NULL
#' @param lower unnamed vector of lower parameter limits, default: NULL
#' @param upper unnamed vector of upper parameter limits, default: NULL
#' @param random vector of parameter(s) to be random effects, default: NULL
run_model <- function(model, data, pars, map=NULL, lower=NULL, upper=NULL, random = NULL) {
    obj =  RTMB::MakeADFun(cmb(model, data), 
                           pars,
                           map = map,
                           random = random)
 
  if(!is.null(lower) & !is.null(upper)) {
    fit = nlminb(start = obj$par,
                 objective = obj$fn,
                 gradient = obj$gr,
                 control = list(iter.max=100000,
                                eval.max=20000),
                 lower = lower,
                 upper=upper)
  }  
  fit = nlminb(start = obj$par,
               objective = obj$fn,
               gradient = obj$gr,
               control = list(iter.max=100000,
                              eval.max=20000))
  rpt = obj$report(obj$env$last.par.best)
  proj = proj_bio(rpt) # function to project the next 2 years
  list(obj=obj, fit=fit, rpt=rpt, proj=proj) 
}


# use Grant's test for model letter or number
model_test <- function(m1, m2) {
  if(sqrt(sum(((m1$rpt$spawn_bio / m2$rpt$spawn_bio - 1)^2) / length(m2$rpt$years))) < 0.1) {
    return('letter')
  } else {
    return ("number")
  }
}

filt <- function(item, data) {
  newd = data
  newd[[item]] = 0
  newd
}


#' Profile a parameter using RTMB with re-optimization
#'
#' @param par_name parameter to profile (e.g., "log_M")
#' @param par_values values to profile the parameter over
#' @param derived character vector of derived quantity likelihoods to extract from report
#' @param data input data to RTMB model
#' @param map any other parameters that are fixed
#' @param obj RTMB model object function (from your original model_fn)
#' @param fit Original fit from rtmb::fit()
#' @param data The data input used in the model
#' @param model The model to fit the data to
#'
#' @return A data frame with parameter value, log-likelihood, and derived values
#' @export
profiles <- function(par_name = "log_M", par_values = log(seq(0.03, 0.1, 0.005)), 
                     derived = c("like_fish_age", "like_srv_age"), 
                     data, map = NULL, obj, fit, model) {
  
  results <- data.frame(value = par_values,
                        log_like = NA_real_,
                        matrix(NA_real_, nrow = length(par_values), ncol = length(derived)))
  colnames(results)[3:ncol(results)] <- derived
  
  nms = unique(names(fit$par))
  split_list <- split(fit$par, names(fit$par))
  pars = lapply(split_list, unname)
  pars = pars[nms]
  
  # put any mapped items back into the pars
  if(!is.null(map)) {
    pars[[names(map)]] = obj$report()[[names(map)]]
  }
  # create map if none provided
  if(is.null(map)) {
    map = list()
  }
  map[[par_name]] = factor(NA)
  
  for (i in seq_along(par_values)) {
    cat("Profiling", par_name, "at", par_values[i], "\n")
    
    # Copy parameter list
    pars_i <- pars
    pars_i[[par_name]] <- par_values[i]
    
    # Rebuild obj with fixed parameter
    obj_i  = RTMB::MakeADFun(cmb(model, data),,
                             parameters = pars_i,
                             map = map)
    
    # Optimize the rest
    fit_i <- nlminb(start = obj_i$par,
                    objective = obj_i$fn,
                    gradient = obj_i$gr,
                    control = list(iter.max = 100000, eval.max = 20000))
    
    results$log_like[i] <- fit_i$objective 
    
    # store raw derived values
    rep <- obj_i$report(fit_i$par)
    for (j in seq_along(derived)) {
      key <- derived[j]
      val <- rep[[key]]
      results[i, key] <- if (is.atomic(val) && length(val) > 1) tail(val, 1) else val
    }
  }
  
  results 

}



#' plot parameter profile
#'
#' @param data profile data.frame with columns value, log_like, ...
#' @param exp  exponentiate the parameter value?
#' @export
plot_profile <- function(data, exp = FALSE) {
  if(isTRUE(exp)) {
    data %>%
      dplyr::mutate(value = exp(value)) -> data
  }
  data %>% 
    dplyr::mutate(dplyr::across(-value, ~ .x - min(.x))) %>%
    tidyr::pivot_longer(-value, values_to = "nll") %>%
    ggplot2::ggplot(ggplot2::aes(value, nll, color = name)) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 0, lty=3) +
    # ggplot2::coord_cartesian(y = c(0,3)) +
    scico::scale_color_scico_d(palette = 'roma') +
    ylab("Change in -log likelihood")
}

check_model_fits <- function(fit, pars) {
  # Basic convergence info
  cat("=== Convergence Information ===\n")
  cat("Convergence code:", fit$convergence, "\n")
  cat("Message:", fit$message, "\n")
  cat("Iterations:", fit$iterations, "\n")
  cat("Objective value:", fit$objective, "\n\n")
  
  # Parameter bounds check
  cat("=== Parameter Bounds Check ===\n")
  at_bounds <- which(abs(fit$par - unlist(lower)) < 1e-5 | 
                       abs(fit$par - unlist(upper)) < 1e-5)
  if(length(at_bounds) > 0) {
    cat("Parameters at bounds:", names(fit$par)[at_bounds], "\n\n")
  } else {
    cat("No parameters at bounds\n\n")
  }
  
  # Gradient check
  cat("=== Gradient Information ===\n")
  final_grad <- obj$gr(fit$par)
  max_grad <- max(abs(final_grad))
  cat("Maximum absolute gradient:", max_grad, "\n")
  if(max_grad > 0.001) {
    cat("WARNING: Large gradients present\n\n")
  }
  
  # Compare initial vs final values
  cat("=== Parameter Changes ===\n")
  initial <- unlist(pars)
  final <- fit$par
  rel_change <- abs((final - initial)/initial)
  changed <- which(rel_change > 0.5)
  if(length(changed) > 0) {
    changes_df <- data.frame(
      Parameter = names(initial)[changed],
      Initial = initial[changed],
      Final = final[changed],
      Rel_Change = rel_change[changed]
    )
    print(changes_df[order(-changes_df$Rel_Change), ])
  } else {
    cat("No parameters changed by more than 50%\n")
  }
}


# jitter ----
jitter_pars <- function(pars, jitter_amt = 0.1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  jitter_value = function(x, jitter_amt) {
    if (length(x) > 1) {
      return(x + rnorm(length(x), mean = 0, sd = jitter_amt))
    } else {
      return(x + rnorm(1, mean = 0, sd = jitter_amt))
    }
  }
  # Apply jitter to each parameter
  jits <- lapply(pars, function(p) jitter_value(p, jitter_amt))
  
  # Preserve names
  names(jits) <- names(pars)
  return(jits)
  
}



# jittered_pars_list <- list() 
# results_df <- data.frame() 
# 
# for (i in 1:10) {
#   message("Running replicate ", i)
#   
#   jittered <- jitter_pars(pars, seed = 100 + i)
#   
#   # Save jittered parameters
#   jittered_pars_list[[i]] <- jittered
#   
#   obj <- RTMB::MakeADFun(cmb(bridge, data), 
#                          jittered)
#   
#   fit <- try(nlminb(start = obj$par,
#                     objective = obj$fn,
#                     gradient = obj$gr,
#                     control = list(iter.max = 100000, eval.max = 20000)), silent = TRUE)
#   
#   if (inherits(fit, "try-error")) {
#     warning(paste("Optimization failed at replicate", i))
#     next
#   }
#   
#   rpt <- obj$report(obj$env$last.par.best)
#   prj <- proj_bio(rpt)[1, ]
#   
#   like_vals <- c(rpt$ssqcatch, rpt$like_srv, rpt$like_fish_age, rpt$like_srv_age, 
#                  rpt$like_fish_size, rpt$like_rec, rpt$f_regularity,
#                  rpt$sprpen, rpt$nll_M, rpt$nll_q, rpt$nll_sigmaR)
#   
#   total_nll <- sum(like_vals)
#   
#   results_df <- rbind(
#     results_df,
#     data.frame(
#       replicate = i,
#       ssqcatch = round(rpt$ssqcatch, 4),
#       like_srv = round(rpt$like_srv, 4),
#       like_fish_age = round(rpt$like_fish_age, 4),
#       like_srv_age = round(rpt$like_srv_age, 4),
#       like_fish_size = round(rpt$like_fish_size, 4),
#       like_rec = round(rpt$like_rec, 4),
#       f_regularity = round(rpt$f_regularity, 4),
#       sprpen = round(rpt$sprpen, 4),
#       nll_M = round(rpt$nll_M, 4),
#       nll_q = round(rpt$nll_q, 4),
#       nll_sigmaR = round(rpt$nll_sigmaR, 4),
#       total_nll = round(total_nll, 4)
#     )
#   )
#   results_df  %>% 
#     flextable::flextable()
# }

effn <- function(obs, pred) {
  
  colSums((1 - pred) * pred) / colSums((obs - pred)^2) 
}

sdnr <- function(obs, pred, iss) {
  n = ncol(obs) 
  sdnr = vector(length = n)
  for(i in 1:n) {
    
    sdnr[i] = sd((obs[,i] - pred[,i]) / sqrt(pred[,i] * (1 - pred[,i]) / iss[i]))
  }
  sdnr
}


osa <- function(obs, pred, iss, yrs, ind, label = 'Age', outlier=3) {
  
  obs = round(iss * obs / colSums(obs))
  pred = pred / colSums(pred)
  res = compResidual::resMulti(obs, pred)
  mat = matrix(res, nrow=nrow(res), ncol=ncol(res))
  df = as.data.frame(mat)
  names(df) <- yrs
  
  df %>%
    mutate(ind = head(ind, -1)) %>%
    tidyr::pivot_longer( -ind, names_to = 'year') %>%
    mutate(year = as.numeric(year),
           Year = factor(year),
           id = ifelse(value<0 , 'a', 'b'),
           Outlier = ifelse(abs(value) >= outlier, "Yes", "No"),
           Outlier = factor(Outlier, levels = c('No', 'Yes'))) -> df
  
  df %>%
    ggplot(aes(year, ind, color = value, size = value, shape = Outlier) ) +
    geom_point(show.legend=TRUE) +
    scale_size_area(guide=F, max_size = 3) +
    scico::scale_color_scico(limits = c(-4, 4), palette = 'vik') +
    afscassess::scale_y_tickr(data = df, var = ind, start=0) +
    afscassess::scale_x_tickr(data = df, var = year) +
    scale_shape_manual(values = c(19,8), drop = FALSE) +
    ylab(label) +
    xlab('Year') +
    ggtitle('OSA') -> osa
  
  df %>% 
    ggplot() + 
    stat_qq(aes(sample=value)) +
    geom_abline(slope=1, intercept = 0, lty=3) -> qq
  list(osa=osa, qq=qq)
  
}
pearson <- function(obs, pred, iss, yrs, ind, label = 'Age', outlier=3) {
  
  as.data.frame(iss * (obs -pred) / sqrt(iss * pred)) %>%
    mutate(ind = ind) -> df
  names(df) <- c(yrs, 'ind')
  
  df %>%
    tidyr::pivot_longer( -ind, names_to = 'year') %>%
    mutate(year = as.numeric(year),
           Year = factor(year),
           id = ifelse(value<0 , 'a', 'b'),
           Outlier = ifelse(abs(value) >= outlier, "Yes", "No"),
           Outlier = factor(Outlier, levels = c("No", "Yes"))) -> df
  
  
  df %>%
    ggplot(aes(year, ind, color = value, size = value, shape = Outlier) ) +
    geom_point(show.legend=TRUE) +
    scale_size_area(guide=F, max_size = 3) +
    scico::scale_color_scico(limits = c(-4, 4), palette = 'vik') +
    afscassess::scale_y_tickr(data = df, var = ind, start=0) +
    afscassess::scale_x_tickr(data = df, var = year) +
    scale_shape_manual(values = c(19,8), drop = FALSE) +
    ylab(label) +
    xlab('Year') +
    ggtitle('Pearson')
  
}

#' aggregate residual plot for RTMB model
#'
#' @param obs observed comp
#' @param pred predictedcomp
#' @param ind vector or ages of lengths
#' @param label axis label
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obs = fish_age_obs
#' pred = report1$fish_age_pred
#' ind = ages
#' agg <- function(obs, pred, ind, label = 'Age')
#' }
agg <- function(obs, pred, ind, label = 'Age') {
  
  df = data.frame(obs = rowSums(obs)/sum(obs),
                  pred = rowSums(pred)/sum(pred),
                  ind = ind)
  
  df %>%
    ggplot(aes(ind, pred)) +
    geom_bar(aes(y=obs), stat = 'identity', alpha=0.4) +
    geom_point() +
    geom_line() +
    afscassess::scale_x_tickr(data=df, var=ind, start=0) +
    xlab(label) +
    ylab('') +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
}

annual <- function(obs, pred, ind, yrs, label = 'Age') {
  obs = as.data.frame(obs)
  pred = as.data.frame(pred)
  names(obs) <- names(pred) <- yrs
  obs %>% 
    dplyr::mutate(type = 'obs',
                  ind = ind) %>% 
    bind_rows(
      pred %>% 
        dplyr::mutate(type = 'pred',
                      ind = ind)
    ) %>% 
    tidyr::pivot_longer(-c(ind, type)) %>% 
    mutate(year = as.numeric(name),
           cohort = year - ind) -> df
  
  df %>% 
    filter(type=='obs') %>% 
    ggplot(aes(ind, value)) + 
    geom_col(alpha = 0.7, aes(fill = factor(cohort)), show.legend = FALSE) + 
    geom_line(data = filter(df, type=='pred')) +
    facet_wrap(~year, ncol = 2, dir='v') + 
    afscassess::scale_x_tickr(data=df, var=ind, start=0) +
    scico::scale_fill_scico_d(palette = 'roma') +
    xlab(label) +
    ylab("") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
}
qq <- function() {
  
}

sample_size <- function(obs, pred, iss, yrs){
  data.frame(year = yrs,
             ISS = iss,
             effN = effn(obs, pred),
             sdnr = sdnr(obs, pred, iss)) %>% 
    tidyr::pivot_longer(-year) %>% 
    mutate(grp = ifelse(name=='sdnr', 'SDNR', 'Sample size')) -> df
  
  df %>% 
    ggplot(aes(year, value, color = name)) +
    geom_point() +
    facet_wrap(~grp, scales = 'free_y', dir = 'h') +
    scale_color_manual("", breaks = c('effN', 'ISS'), values = c("#7E1700","#5DC0D2",1)) +
    expand_limits(y = 0) +
    theme(legend.position=c(0.1,0.8)) +
    afscassess::scale_x_tickr(data=df, var=year, to=10, start = 1960) +
    xlab('Year') +
    ylab('Value')
}
#' get all residual plots for RTMB model
#'
#' @param obs observed comp
#' @param pred predictedcomp
#' @param iss input sample size
#' @param yrs comp years
#' @param ind vector or ages of lengths
#' @param label axis label
#'
#' @export
#'
#' @examples
#' \dontrun{
#' iss = data$fish_age_iss
#' obs = fish_age_obs
#' pred = report1$fish_age_pred
#' yrs = fish_age_yrs
#' label = "Age"
#' ind = ages
#' fish_age_resids <- resids(obs = fish_age_obs, pred = report1$fish_age_pred, iss = data$fish_age_iss,
#'                          yrs = fish_age_yrs, ind = ages)
#'  fish_age_resids$osa + ggtitle('osa fishery age comp residuals')
#' }
resids <- function(obs, pred, iss, yrs, ind, label = 'Age', outlier=3) {
  list(osa = osa(obs, pred, iss, yrs, ind, label, outlier)$osa,
       qq = osa(obs, pred, iss, yrs, ind, label, outlier)$qq,
       pearson = pearson(obs, pred, iss, yrs, ind, label, outlier),
       agg = agg(obs, pred, ind, label),
       annual = annual(obs, pred, ind, yrs, label),
       ss = sample_size(obs, pred, iss, yrs) )
}
