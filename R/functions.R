##### This script contains all the functions for the targets pipeline

# Preparing data
prepare_isometric_data <- function(data) {
  isometric_data <- data |>
  select(1:15, contains("max")) |>
  pivot_longer(16:33, 
               names_to = "name",
               values_to = "value") |>
  separate(name,
           into = c("exercise", "x", "day", "trial")) |>
  select(-x) |> 
  mutate_if(is.character,as.factor) |>
  group_by(participant, exercise) |>
  mutate(mean_i = mean(value),
         mean_diff = value - mean_i)
  
  contrasts(isometric_data$day) <- contr.helmert(3)
  
  isometric_data
}

prepare_isokinetic_data <- function(data) {
  isokinetic_data <- data |>
  select(1:15, contains("best")) |>
  pivot_longer(16:27, 
               names_to = "name",
               values_to = "value") |>
  separate(name,
           into = c("day", "exercise", "x", "y", "z", "a", "con_ecc")) |>
  select(-x,-y,-z,-a) |>
  mutate(day = case_when(
    day == "t1" ~ "d1",
    day == "t2" ~ "d2"
  )) |> 
  mutate_if(is.character,as.factor) |>
    pivot_wider(names_from = c("day", "exercise", "con_ecc"),
                values_from = "value",
                id_cols = c("participant")) |>
    mutate(delta_cp_con = d2_cp_con - d1_cp_con,
           mean_cp_con = (d1_cp_con + d2_cp_con)/1,
           delta_lp_con = d2_lp_con - d1_lp_con,
           mean_lp_con = (d1_lp_con + d2_lp_con)/1,
           delta_row_con = d2_row_con - d1_row_con,
           mean_row_con = (d1_row_con + d2_row_con)/1,
           delta_cp_ecc = d2_cp_ecc - d1_cp_ecc,
           mean_cp_ecc = (d1_cp_ecc + d2_cp_ecc)/1,
           delta_lp_ecc = d2_lp_ecc - d1_lp_ecc,
           mean_lp_ecc = (d1_lp_ecc + d2_lp_ecc)/1,
           delta_row_ecc = d2_row_ecc - d1_row_ecc,
           mean_row_ecc = (d1_row_ecc + d2_row_ecc)/1) 
  
  isokinetic_data
}

# Misc functions
make_plot_tiff <- function(plot, width, height, path) {
  ggsave(
    path,
    plot,
    width = width,
    height = height,
    device = "tiff",
    dpi = 300
  )
  
}

# Isometric models and plots
fit_isometric_model <- function(isometric_data) {
  isometric_data_wide <- isometric_data |>
    select(participant, exercise, day, trial, value) |>
    pivot_wider(names_from = "exercise",
                values_from = "value",
                id_cols = c("participant", "day", "trial"))
  
  isometric_cp_bf <- bf(cp ~ 1 + day + (1|a|participant) + (1|b|participant:day))
  isometric_lp_bf <- bf(lp ~ 1 + day + (1|a|participant) + (1|b|participant:day))
  isometric_row_bf <- bf(row ~ 1 + day + (1|a|participant) + (1|b|participant:day))
  
  brm_model_isometric <- brm(isometric_cp_bf + isometric_lp_bf + isometric_row_bf + set_rescor(TRUE),
                             data = isometric_data_wide,
                             chains = 4,
                             cores = 4,
                             seed = 1988,
                             warmup = 2000,
                             iter = 8000,
                             control = list(adapt_delta = 0.99),
                             save_pars = save_pars(all = TRUE))
}

get_isometric_draws <- function(isometric_model, isometric_data) {
  draws_isometric <- spread_draws(isometric_model,
                                  b_cp_day1, b_cp_day2,
                                  `sd_participant:day__cp_Intercept`,
                                  sigma_cp,
                                  b_lp_day1, b_lp_day2,
                                  `sd_participant:day__lp_Intercept`,
                                  sigma_lp,
                                  b_row_day1, b_row_day2,
                                  `sd_participant:day__row_Intercept`,
                                  sigma_row
  ) |>
    mutate(
      cp_loam_bw = 1.96 * sqrt( ((length(unique(isometric_data$day))-1)/length(unique(isometric_data$day)))*`sd_participant:day__cp_Intercept`^2 + 
                                  ((length(unique(isometric_data$day))*length(unique(isometric_data$trial))-1)/(length(unique(isometric_data$day))*length(unique(isometric_data$trial))))*(sigma_cp^2)),
      
      cp_loam_w = 1.96 * sqrt( ((length(unique(isometric_data$trial))*length(unique(isometric_data$day))-1)/(length(unique(isometric_data$trial))*length(unique(isometric_data$day))))*sigma_cp^2),
      
      lp_loam_bw = 1.96 * sqrt( ((length(unique(isometric_data$day))-1)/length(unique(isometric_data$day)))*`sd_participant:day__lp_Intercept`^2 + 
                                  ((length(unique(isometric_data$day))*length(unique(isometric_data$trial))-1)/(length(unique(isometric_data$day))*length(unique(isometric_data$trial))))*(sigma_lp^2)),
      
      lp_loam_w = 1.96 * sqrt( ((length(unique(isometric_data$trial))*length(unique(isometric_data$day))-1)/(length(unique(isometric_data$trial))*length(unique(isometric_data$day))))*sigma_lp^2),
      
      row_loam_bw = 1.96 * sqrt( ((length(unique(isometric_data$day))-1)/length(unique(isometric_data$day)))*`sd_participant:day__row_Intercept`^2 + 
                                   ((length(unique(isometric_data$day))*length(unique(isometric_data$trial))-1)/(length(unique(isometric_data$day))*length(unique(isometric_data$trial))))*(sigma_row^2)),
      
      row_loam_w = 1.96 * sqrt( ((length(unique(isometric_data$trial))*length(unique(isometric_data$day))-1)/(length(unique(isometric_data$trial))*length(unique(isometric_data$day))))*sigma_row^2)
    ) 
}

get_isometric_bias_summary <- function(isometric_draws) {
  isometric_bias_summary <- isometric_draws |>
    select(b_cp_day1, b_cp_day2, b_lp_day1, b_lp_day2, b_row_day1, b_row_day2) |>
    pivot_longer(1:6, 
                 names_to = "coef",
                 values_to = "mean")  |>
    separate(coef, into = c("coef", "exercise", "day")) |>
    unite("coef", c("coef","day"), sep = "_") |>
    mutate(coef = case_when(
      coef == "b_day1" ~ "Day 1 to Day 2",
      coef == "b_day2" ~ "Day 2 to Day 3"
    ),
    exercise = case_when(
      exercise == "cp" ~ "Chest Press",
      exercise == "lp" ~ "Leg Press",
      exercise == "row" ~ "Row"
    )) |>
    group_by(coef, exercise) |>
    mean_qi(mean)
}

get_isometric_loamr_summary <- function(isometric_data, isometric_draws) {
  isometric_loamr_summary <- isometric_draws |>
    select(cp_loam_bw, cp_loam_w, lp_loam_bw, lp_loam_w, row_loam_bw, row_loam_w) |>
    pivot_longer(1:6, 
                 names_to = "coef",
                 values_to = "mean") |> 
    group_by(coef) |>
    mean_qi() |>
    separate(coef, into = c("exercise", "loam", "bw_w")) |>
    left_join(isometric_data |> select(exercise, mean_i), by = "exercise") |>
    mutate(exercise = case_when(
      exercise == "cp" ~ "Chest Press",
      exercise == "lp" ~ "Leg Press",
      exercise == "row" ~ "Row"
    ))
}

plot_isometric_bias <- function(isometric_draws, isometric_bias_summary) {
  bias_plot <- isometric_draws |>
    select(b_cp_day1, b_cp_day2, b_lp_day1, b_lp_day2, b_row_day1, b_row_day2) |>
    pivot_longer(1:6, 
                 names_to = "coef",
                 values_to = "draw")  |>
    separate(coef, into = c("coef", "exercise", "day")) |>
    unite("coef", c("coef","day"), sep = "_") |>
    mutate(coef = case_when(
      coef == "b_day1" ~ "Day 1 to Day 2",
      coef == "b_day2" ~ "Day 2 to Day 3"
    ),
    exercise = case_when(
      exercise == "cp" ~ "Chest Press",
      exercise == "lp" ~ "Leg Press",
      exercise == "row" ~ "Row"
    )) |>
    ggplot(aes(x = coef, y = draw)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    stat_slabinterval(alpha = 0.75, .width = .95) +
    geom_text(data = isometric_bias_summary,
              aes(x=coef, y=mean, label=round(mean,2)), nudge_x=-0.2, size = 2.5) +
    geom_text(data = isometric_bias_summary,
              aes(x=coef, y=.lower, label=round(.lower,2)), nudge_x=-0.2, size = 2.5) +
    geom_text(data = isometric_bias_summary,
              aes(x=coef, y=.upper, label=round(.upper,2)), nudge_x=-0.2, size = 2.5) +
    facet_grid(.~exercise) +
    labs(
      title = "Between-day mean bias",
      x = "Contrast Coefficients for Day (Helmert Coding)",
      y = expression(beta[j]-beta[j-1]),
      color = "Day"
    ) +
    theme_bw()
}

plot_isometric_loamr <- function(isometric_data, isometric_loamr_summary) {
  loam_plot <- isometric_data |>
    mutate(exercise = case_when(
      exercise == "cp" ~ "Chest Press",
      exercise == "lp" ~ "Leg Press",
      exercise == "row" ~ "Row"
    )) |>
    ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    geom_ribbon(data = isometric_loamr_summary |> filter(bw_w == "w"),
                aes(x=mean_i, ymin = .lower, ymax = .upper, group=exercise),
                alpha = 0.75, fill = "#E69F00") +
    geom_ribbon(data = isometric_loamr_summary |> filter(bw_w == "w"),
                aes(x=mean_i, ymin = -.lower, ymax = -.upper, group=exercise),
                alpha = 0.75, fill = "#E69F00") +
    geom_line(data = isometric_loamr_summary |> filter(bw_w == "w"),
              aes(x=mean_i, y = mean, group=exercise),
              alpha = 0.75, linetype = "dotted") +
    geom_line(data = isometric_loamr_summary |> filter(bw_w == "w"),
              aes(x=mean_i, y = -mean, group=exercise), linetype = "dotted") +
    
    
    geom_ribbon(data = isometric_loamr_summary |> filter(bw_w == "bw"),
                aes(x=mean_i, ymin = .lower, ymax = .upper, group=exercise),
                alpha = 0.75, fill = "#56B4E9") +
    geom_ribbon(data = isometric_loamr_summary |> filter(bw_w == "bw"),
                aes(x=mean_i, ymin = -.lower, ymax = -.upper, group=exercise),
                alpha = 0.75, fill = "#56B4E9") +
    geom_line(data = isometric_loamr_summary |> filter(bw_w == "bw"),
              aes(x=mean_i, y = mean, group=exercise)) +
    geom_line(data = isometric_loamr_summary |> filter(bw_w == "bw"),
              aes(x=mean_i, y = -mean, group=exercise)) +
    
    geom_text(data = isometric_loamr_summary |> filter(bw_w == "w") |> group_by(exercise) |> 
                mutate(x = (max(mean_i) - min(mean_i))/2) |> slice_min(mean_i, n=1) |> slice_head(n=1),
              aes(x=mean_i + x, y=-75, group=exercise, 
                  label=glue::glue("Within-day: +/- {round(mean,2)} [95% QI: {round(.lower,2)},{round(.upper,2)}]")),
              size = 2.5) +
    geom_text(data = isometric_loamr_summary |> filter(bw_w == "bw") |> group_by(exercise) |> 
                mutate(x = (max(mean_i) - min(mean_i))/2) |> slice_min(mean_i, n=1) |> slice_head(n=1),
              aes(x=mean_i + x, y=-85, group=exercise, 
                  label=glue::glue("Between-day: +/- {round(mean,2)} [95% QI: {round(.lower,2)},{round(.upper,2)}]")),
              size = 2.5) +
    
    geom_point(aes(x=mean_i, y=mean_diff), size = 2, color = "black", fill = NA, shape = 21, alpha = 0.75) +
    geom_point(aes(x=mean_i, y=mean_diff, color=day), alpha = 0.75) +
    scale_color_viridis_d() +
    scale_y_continuous(breaks = seq(-100, 100, by=25)) +
    ggh4x::facet_grid2(.~exercise, scales = "free", independent = "x") +
    labs(
      title = "Limits of Agreement with the Mean",
      x = expression(bar(y)[i..]),
      y = expression(y[ijk]-bar(y)[i..]),
      color = "Day"
    ) +
    theme_bw()
    
}

# Isokinetic models and plots
fit_isokinetic_model <- function(isokinetic_data) {
  
  cp_con_bf <- bf(delta_cp_con ~ 1)
  
  lp_con_bf <- bf(delta_lp_con ~ 1)
  
  row_con_bf <- bf(delta_row_con ~ 1)
  
  cp_ecc_bf <- bf(delta_cp_ecc ~ 1)
  
  lp_ecc_bf <- bf(delta_lp_ecc ~ 1)
  
  row_ecc_bf <- bf(delta_row_ecc ~ 1)
  
  
  brm_model_isokinetic <- brm(cp_con_bf + lp_con_bf + row_con_bf + 
                                cp_ecc_bf + lp_ecc_bf + row_ecc_bf + 
                                set_rescor(rescor = TRUE),
                              data = isokinetic_data,
                              chains = 4,
                              cores = 4,
                              seed = 1988,
                              warmup = 2000,
                              iter = 8000,
                              control = list(adapt_delta = 0.99),
                              save_pars = save_pars(all = TRUE))
  
}

get_isokinetic_draws <- function(isokinetic_model) {
  draws_isokinetic <- gather_draws(isokinetic_model,
                                   b_deltacpcon_Intercept, b_deltalpcon_Intercept, b_deltarowcon_Intercept,
                                   b_deltacpecc_Intercept, b_deltalpecc_Intercept, b_deltarowecc_Intercept,
                                   sigma_deltacpcon, sigma_deltalpcon, sigma_deltarowcon,
                                   sigma_deltacpecc, sigma_deltalpecc, sigma_deltarowecc,
  ) |>
    mutate(exercise = case_when(
      str_detect(.variable, pattern = "cp") ~ "cp",
      str_detect(.variable, pattern = "lp") ~ "lp", 
      str_detect(.variable, pattern = "row") ~ "row"
    ),
    con_ecc = case_when(
      str_detect(.variable, pattern = "con") ~ "con", 
      str_detect(.variable, pattern = "ecc") ~ "ecc"
    ),
    .variable = case_when(
      str_detect(.variable, pattern = "b_delta") ~ "bias", 
      str_detect(.variable, pattern = "sigma") ~ "sigma", 
      str_detect(.variable, pattern = "sigma") ~ "sigma"
    )) |>
    pivot_wider(names_from = ".variable",
                values_from = ".value",
                id_cols = c(".chain", ".iteration", ".draw", "exercise", "con_ecc")) |>
    mutate(
      loa = 1.96 * sigma,
    ) 
}

get_isokinetic_bias_loa_summary <- function(isokinetic_draws, isokinetic_data) {
  isokinetic_bias_loa_summary <- isokinetic_draws |>
    group_by(exercise, con_ecc) |>
    mean_qi(bias, loa) |>
    left_join(isokinetic_data |>  select(contains("mean")) |>
                pivot_longer(1:6,
                             names_to = "what",
                             values_to = "mean") |>
                mutate(con_ecc = case_when(
                  str_detect(what, pattern = "con") ~ "con", 
                  str_detect(what, pattern = "ecc") ~ "ecc", 
                ), 
                exercise = case_when(
                  str_detect(what, pattern = "cp") ~ "cp", 
                  str_detect(what, pattern = "lp") ~ "lp", 
                  str_detect(what, pattern = "row") ~ "row"
                )),
              by = c("exercise", "con_ecc")) |>
    mutate(exercise = case_when(
      exercise == "cp" ~ "Chest Press",
      exercise == "lp" ~ "Leg Press",
      exercise == "row" ~ "Row"
    ),
    con_ecc = case_when(
      con_ecc == "con" ~ "Concentric",
      con_ecc == "ecc" ~ "Eccentric"
    ))
}

plot_isokinetic_BA <- function(isokinetic_bias_loa_summary, isokinetic_data) {
  loa_plot_isokinetic <- isokinetic_data |> select(participant, contains("mean"), contains("delta")) |>
    pivot_longer(2:13,
                 names_to = c("coef", "exercise", "con_ecc"),
                 values_to =  "value",
                 names_sep = "_") |>
    pivot_wider(names_from = "coef",
                values_from = "value",
                id_cols = c("participant", "exercise", "con_ecc")) |>
    mutate(exercise = case_when(
      exercise == "cp" ~ "Chest Press",
      exercise == "lp" ~ "Leg Press",
      exercise == "row" ~ "Row"
    ),
    con_ecc = case_when(
      con_ecc == "con" ~ "Concentric",
      con_ecc == "ecc" ~ "Eccentric"
    )) |>
    ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    geom_ribbon(data = isokinetic_bias_loa_summary,
                aes(x=mean, ymin = bias.lower, ymax = bias.upper),
                alpha = 0.5) +
    geom_ribbon(data = isokinetic_bias_loa_summary,
                aes(x=mean, ymin = bias-loa.upper, ymax = bias-loa.lower),
                alpha = 0.5) +
    geom_ribbon(data = isokinetic_bias_loa_summary,
                aes(x=mean, ymin = bias+loa.upper, ymax = bias+loa.lower),
                alpha = 0.5) +
    geom_line(data = isokinetic_bias_loa_summary,
              aes(x=mean, y = bias)) +
    geom_line(data = isokinetic_bias_loa_summary,
              aes(x=mean, y = bias-loa)) +
    geom_line(data = isokinetic_bias_loa_summary,
              aes(x=mean, y = bias+loa)) +
    geom_point(aes(x=mean, y=delta), alpha = 0.75) +
    
    geom_text(data = isokinetic_bias_loa_summary |> group_by(exercise, con_ecc) |> 
                mutate(x = (max(mean) - min(mean))/2, y = (-loa.upper*1.25)) |> slice_min(mean, n=1) |> slice_head(n=1),
              aes(x=mean + x, y=y, group=exercise, 
                  label=glue::glue("Mean bias: {round(bias,2)} [95% QI: {round(bias.lower,2)},{round(bias.upper,2)}]\nLimits of Agreement: +/- {round(loa,2)} [95% QI: {round(loa.lower,2)},{round(loa.upper,2)}]")),
              size = 2.5, vjust=-0.05) +
    scale_y_continuous(breaks = seq(-200, 200, by=25)) +
    ggh4x::facet_grid2(con_ecc~exercise, scales = "free", independent = "all") +
    
    labs(
      title = "Mean Bias and Limits of Agreement",
      caption = "Note, the upper, middle, and lower solid horizontal lines (point estimate) with grey ribbons (95% quantile intervals) show the between day upper 95% limit of agreement, the mean bias, and lower 95% limit of agreement respectively",
      x = expression(bar(y)[i.]),
      y = expression(y[i2]-y[i1])
    ) +
    theme_bw() +
    theme(plot.caption = element_text(size=6))
}