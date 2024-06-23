#### load packages ####

library(afex)
library(dplyr)
library(ggplot2)
library(magrittr)
library(stringr)

#### z scores and outliers #####

z <- function(x){(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)}

z_percent <- function(x, z){(which(abs(x) >= z) %>% length()) / length(!is.na(x)) * 100}

load_z_table <- function(){readr::read_csv("../data/z_table_long.csv")}

# z_table <- readr::read_csv("../data/z_table.csv") %>%
#   dplyr::filter(!is.na(Z)) %>%
#   tidyr::pivot_longer(., names_to = "dec", -Z) %>%
#   dplyr::transmute(
#     z = Z - as.numeric(dec),
#     prob = value
#   ) %>%
#   dplyr::arrange(prob)

#write.csv(z_table, "../data/z_table_long.csv", row.names = FALSE)

#### working with models ####

get_obj_name <- function(object) {
  deparse(substitute(object))
}

get_coef <- function(model, label, model_name
                     ) {
  summary(model)$coefficients %>%
    tibble::as_tibble(., rownames = "coefficient") %>%
    dplyr::transmute(
      model_name = model_name,
      model_type = label,
      coefficient = coefficient,
      b = round(Estimate, digits = 5),
      se = `Std. Error`,
      t = `t value`,
      p = `Pr(>|t|)`,
      ci_lower = confint(model)[1:nrow(summary(model)$coefficients),1],
      ci_upper = confint(model)[1:nrow(summary(model)$coefficients),2]
    )
}

get_coef_se <- function(model, label = "HC4", model_name){

  parameters <- parameters::model_parameters(model, robust = TRUE, vcov_type = "HC4") %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      model_name = model_name,
      model_type = label,
      coefficient = Parameter,
      b = Coefficient,
      se = SE,
      t = t,
      p = p,
      ci_lower = CI_low,
      ci_upper = CI_high
    )

  return(parameters)

}

get_mode_number <- function(vector){

  binned_resid <- vector %>%
    round(.,digits = 2) %>%
    tibble::as_tibble() %>%
    dplyr::group_by(value) %>%
    dplyr::summarise(
      n = dplyr::n()
    )

  mode_n <- binned_resid %>% dplyr::filter(n == max(n))

  binned_resid %<>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      rel_freq = n / mode_n$n[1] * 100,
      sd = sd(value),
      sd_2_3rds = sd/3*2,
      dist_from_mode = (value - mode_n$value[1]) %>% abs(),
      alt_mode = dplyr::if_else(rel_freq >= 80 & dist_from_mode > sd_2_3rds, 1, 0)
    )

  alt_mode_tib <- binned_resid %>% dplyr::filter(alt_mode == 1)

  return(nrow(alt_mode_tib) + 1) # represents the number of modes
}



get_het <- function(z_model, categorical_only = FALSE){

  var_names = variable.names(z_model)[2:length(variable.names(z_model))]
  var_names = stringr::str_remove_all(var_names, pattern = "`")
  var_names_in_data = var_names %in% names(z_model$model)

  if(any(var_names_in_data) == FALSE){

    data = z_model$model %>%
      dplyr::mutate(
        fitted = z_model$fitted.values,
        filter_var = factor(round(fitted, digits = 4)),
        residuals_z_sqrt = sqrt(abs(z(z_model$residuals)))
      )

    levels = levels(data$filter_var)
    level_combinations = gtools::combinations(n = length(levels), r = 2, v = levels)

    slopes <- c()

    for(i in 1:nrow(level_combinations)){

      step_mod_i = lm(residuals_z_sqrt ~ fitted,
                      data = data %>% dplyr::filter(filter_var %in% level_combinations[i, ]))

      cc_i = step_mod_i$coefficients[2]

      slopes = c(slopes, cc_i)

    }

    return(
      data.frame(
        slope = abs(max(slopes)),
        avg_slope = mean(abs(slopes), na.rm = TRUE),
        predictors = "categorical only"
        )
      )

  } else {

    resid_tib <- tibble::tibble(
      residuals = z_model$residuals,
      residuals_z_sqrt = sqrt(abs(z(residuals))),
      fitted = z_model$fitted.values,
    )

    step_mod <- lm(residuals_z_sqrt ~ fitted * I(fitted >= 0), data = resid_tib)

    slope_1 = abs(summary(step_mod)$coefficients[2])
    slope_2 = abs(summary(step_mod)$coefficients[2] + summary(step_mod)$coefficients[4])

    if(slope_1 >= slope_2){
      return(
        data.frame(
          slope = slope_1,
          avg_slope = mean(c(abs(slope_1), abs(slope_1)), na.rm = TRUE),
          predictors = "mixed or numeric only"
        )
      )
    } else {
      return(
        data.frame(
          slope = slope_2,
          avg_slope = mean(c(abs(slope_1), abs(slope_1)), na.rm = TRUE),
          predictors = "mixed or numeric only"
        )
      )
    }
  }
}


fit_and_save <- function(data, model, vars_to_z, paper_id = 9999, model_numer = 99){

  data %<>%
    dplyr::mutate(
      dplyr::across(
        .cols = vars_to_z,
        .fns = z,
        .names = "{.col}"
      )
    )

  model_formula = model$call %>% as.character()

  z_model = lm(
    formula = as.formula(model_formula[2]),
    data = data
  )

  saveRDS(object = z_model,
          file = paste0("../models/ols/", paper_id, "_mod_z_", model_numer, ".rds"))

  robust_model = robustbase::lmrob(
    formula = as.formula(model_formula[2]),
    data = data,
    #setting = "KS2011"
  )

  model_list = list(z_model, robust_model)


  saveRDS(object = robust_model,
          file = paste0("../models/lmrob/", paper_id, "_mod_rob_", model_numer, ".rds"))

  return(model_list)

}

fit_z <- function(data, model, vars_to_z){
  
  data %<>%
    dplyr::mutate(
      dplyr::across(
        .cols = vars_to_z,
        .fns = z,
        .names = "{.col}"
      )
    )
  
  model_formula = model$call %>% as.character()
  
  z_model = lm(
    formula = as.formula(model_formula[2]),
    data = data
  )
  
  return(z_model)
  
}


summarise_and_save <- function(data, var_to_summarise, group_by = NULL, paper_id, summary_number){

  data %<>%
    dplyr::mutate(
      dplyr::across(
        .cols = var_to_summarise,
        .fns = z,
        .names = "{.col}"
      )
    )

  if(is.null(group_by)){

    var_summary = data %>%
      dplyr::group_by() %>%
      dplyr::summarise(
        paper     = paste0(paper_id),
        outcome   = paste0(var_to_summarise),
        n         = dplyr::n(),
        skew      = skewness(get(var_to_summarise), na.rm = TRUE),
        kurt      = kurtosis(get(var_to_summarise), na.rm = TRUE),
        z_1.96    = z_percent(get(var_to_summarise), 1.96),
        z_2.58    = z_percent(get(var_to_summarise), 2.58),
        z_3.29    = z_percent(get(var_to_summarise), 3.29),
        mode_n    = get_mode_number(get(var_to_summarise)),
        var       = NA
      ) %>%
      dplyr::ungroup() %>% #
      dplyr::mutate(
        var_ratio = max(var) / min(var)
      ) %>%
      dplyr::select(variables)

  } else if (length(group_by) == 2){

    var_summary = data %>%
      dplyr::group_by(get(group_by[1]), get(group_by[2])) %>%
      dplyr::summarise(
        paper     = paste0(paper_id),
        outcome   = paste0(var_to_summarise),
        n         = dplyr::n(),
        skew      = skewness(get(var_to_summarise), na.rm = TRUE),
        kurt      = kurtosis(get(var_to_summarise), na.rm = TRUE),
        z_1.96    = z_percent(get(var_to_summarise), 1.96),
        z_2.58    = z_percent(get(var_to_summarise), 2.58),
        z_3.29    = z_percent(get(var_to_summarise), 3.29),
        mode_n    = get_mode_number(get(var_to_summarise)),
        var       = sd(get(var_to_summarise), na.rm = TRUE)^2
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        var_ratio = max(var) / min(var)
      ) %>%
      dplyr::select(variables)

  } else if (length(group_by) == 3){

    var_summary = data %>%
      dplyr::group_by(get(group_by[1]), get(group_by[2]), get(group_by[3])) %>%
      dplyr::summarise(
        paper     = paste0(paper_id),
        outcome   = paste0(var_to_summarise),
        n         = dplyr::n(),
        skew      = skewness(get(var_to_summarise), na.rm = TRUE),
        kurt      = kurtosis(get(var_to_summarise), na.rm = TRUE),
        z_1.96    = z_percent(get(var_to_summarise), 1.96),
        z_2.58    = z_percent(get(var_to_summarise), 2.58),
        z_3.29    = z_percent(get(var_to_summarise), 3.29),
        mode_n    = get_mode_number(get(var_to_summarise)),
        var       = sd(get(var_to_summarise), na.rm = TRUE)^2
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        var_ratio = max(var) / min(var)
      ) %>%
      dplyr::select(variables)

  } else {

    var_summary = data %>%
      dplyr::group_by(get(group_by[1])) %>%
      dplyr::summarise(
        paper     = paste0(paper_id),
        outcome   = paste0(var_to_summarise),
        n         = dplyr::n(),
        skew      = skewness(get(var_to_summarise), na.rm = TRUE),
        kurt      = kurtosis(get(var_to_summarise), na.rm = TRUE),
        z_1.96    = z_percent(get(var_to_summarise), 1.96),
        z_2.58    = z_percent(get(var_to_summarise), 2.58),
        z_3.29    = z_percent(get(var_to_summarise), 3.29),
        mode_n    = get_mode_number(get(var_to_summarise)),
        var       = sd(get(var_to_summarise), na.rm = TRUE)^2
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        var_ratio = max(var) / min(var)
      ) %>%
      dplyr::select(variables)

  }

  write.csv(x = var_summary,
            file = paste0("../summaries/", paper_id, "_z_", var_to_summarise, "_summary_", summary_number, ".csv"),
            row.names = FALSE)

  return(var_summary)

}

get_model_info <- function(z_model, model_name){

  het_vals <- get_het(z_model)

  tibble::tibble(
    model_name = model_name,
    outcome   = as.character(z_model$call$formula)[2],
    n         = z_model$residuals %>% length(),
    skew      = skewness(z_model$residuals, na.rm = TRUE),
    kurt      = kurtosis(z_model$residuals, na.rm = TRUE),
    z_1.96    = z_percent(z_model$residuals, 1.96),
    z_2.58    = z_percent(z_model$residuals, 2.58),
    z_3.29    = z_percent(z_model$residuals, 3.29),
    mode_n    = get_mode_number(z_model$residuals),
    het_slope = het_vals$slope,
    het_avg_slope = het_vals$avg_slope,
    predictor_class = het_vals$predictors
  )

}

#### other ####

loop_progress <- function(){

  Sys.sleep(0.001)
  print(i)
  flush.console()

}

t2 <- function(x){
  x %>% t() %>% t()
}

reverse_score <- function(var, max_score){
  return((var - (max_score + 1)) %>% abs())
}

#### Other small helpers ####

extract_num <- function(x){
  as.numeric(gsub("[^0-9.\\-]+", "", as.character(x)))
}

weight_by_n <- function(x, w){
  return(x*w)
}


export_ols <- function(mod_z, key_effects) {
  
  mod_z_list <- list(
    model = mod_z, 
    key_effects = key_effects
  )
  
  file = paste0("../objects/lm/", 
                deparse(substitute(mod_z)), 
                "_list.rds")
  
  saveRDS(object = mod_z_list, 
          file = file)
  
}


reset_df <- function(){
  
  set.seed(1234)
  
  n = 200
  h_level = 0.2
  exp = 2
  
  b0 = 20
  b1 = 1.5
  x_sd = 20
  
  x = rnorm(n = n, mean = 50, sd = x_sd)
  
  h = function(x) x_sd + (h_level * x)^exp
  sd = h(x)
  
  error = rnorm(n = n, mean = 0, sd = sd)
  
  y = b0 + b1*x + error
  
  df = data.frame(
    x = x,
    y = y,
    error = error
  )
  
  # p1 <- df %>%
  #   ggplot2::ggplot(., aes(x = x, y = y)) +
  #   stat_smooth(formula = "y~x", method = "loess") +
  #   geom_point() +
  #   theme_light()
  
  mod <- lm(y ~ x, df)
  # plot(mod, which = c(1,3))
  
  df %<>% 
    dplyr::mutate(
      resid = mod$residuals,
      z_resid = z(mod$residuals), 
      z_predict = z(mod$fitted.values)
    )
  
  df <<- df
}


progress <- function(){
  
  mods <- list.files(here::here("objects/lm")) %>% 
    as.data.frame()
  
  mods %<>% 
    tidyr::separate(., col = ., into = c("id", "."), sep = "_mod_") %>% 
    dplyr::filter(!str_detect(., "lmer")) %>% 
    dplyr::mutate(
      design = if_else(str_detect(., "rm"), "within", "between"), 
      mod_num = readr::parse_number(.), 
      pub_status = if_else(str_detect(id, "u"), "unpublished", "published")
    ) %>% 
    dplyr::transmute(id, mod_num, design, pub_status)
  
  n_models <- nrow(mods)
  n_papers <- length(unique(mods$id))
  
  `Counts by publication` <- mods %>% 
    dplyr::group_by(id) %>% 
    dplyr::filter(row_number() %in% 1) %>% 
    group_by(pub_status) %>% 
    dplyr::summarise(n = n()) 
  
  `Counts by design` <- mods %>% 
    group_by(design) %>% 
    dplyr::summarise(n = n()) 
  
  `Counts by publication and design` <- mods %>% 
    group_by(design, pub_status) %>% 
    dplyr::summarise(n = n()) %>%
    dplyr::arrange(n)
  
  `Counts by model count` <- mods %>% 
    group_by(id) %>% 
    dplyr::summarise(n = n()) %>%
    dplyr::arrange(n)
  
  
  
  cat("\n", "Papers re-analysed: ", n_papers, "\n", 
      "Models total: ", n_models, "\n")
  
  list(
    n_models = n_models, 
    n_papers = n_papers,
    mods = mods,
    `Counts by publication` = `Counts by publication`, 
    `Counts by design` = `Counts by design`, 
    `Counts by publication and design` = `Counts by publication and design`, 
    `Counts by model count` = `Counts by model count`
  )
  
  
}


mods_info_export <- function(mod_id)
{
  mods <- progress() 
  id_mod <- mods$mods %>%
    dplyr::filter(id == mod_id) %>% 
    dplyr::arrange(mod_num)
  
  write.csv(id_mod, paste0(mod_id, "_models.csv"))
}

mod_coeffs <- function(mod){
  rownames(summary(mod)$coefficients)[-1]
}

gmc <- function(x) {x - mean(x)}

#### authors' functions ####

zscore <- function (v) 
{
  (v - mean(v, na.rm = T))/sqrt(var(v, na.rm = T))
}


























