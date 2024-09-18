##### This script is just a sandbox for building and checking functions before adding to pipeline

library(tidyverse)
library(here)
library(janitor)
library(loamr)
library(patchwork)
library(rstan)
library(brms)
library(tidybayes)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

data <- read_csv(here("data","data.csv")) |>
  clean_names()

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
  mutate_if(is.character,as.factor)

contrasts(isokinetic_data$day) <- contr.helmert(3)



# Look at descriptives
isometric_data |>
  group_by(exercise, day, trial) |>
  summarise(mean = mean(value),
            sd = sd(value)) |>
  ggplot(aes(x=day, y = mean)) +
  geom_pointrange(aes(y=mean, ymin = mean-sd, ymax = mean+sd)) +
  facet_grid(exercise~trial, scales = "free_y")

isokinetic_data |>
  group_by(exercise, day, con_ecc) |>
  summarise(mean = mean(value),
            sd = sd(value)) |>
  ggplot(aes(x=day, y = mean)) +
  geom_pointrange(aes(y=mean, ymin = mean-sd, ymax = mean+sd)) +
  facet_grid(exercise~con_ecc, scales = "free_y")

# Bayesian LoA with the mean approach 

# Try a multivariate model
isometric_data_wide <- isometric_data |>
  select(participant, exercise, day, trial, value) |>
  pivot_wider(names_from = "exercise",
              values_from = "value",
              id_cols = c("participant", "day", "trial"))

isometric_cp_bf <- bf(cp ~ 1 + day + (1|a|participant) + (1|b|participant:day))
isometric_lp_bf <- bf(lp ~ 1 + day + (1|a|participant) + (1|b|participant:day))
isometric_row_bf <- bf(row ~ 1 + day + (1|a|participant) + (1|b|participant:day))

isometric_model <- brm(isometric_cp_bf + isometric_lp_bf + isometric_row_bf,
                              data = isometric_data_wide,
                              chains = 4,
                              cores = 4,
                              seed = 1988,
                              warmup = 2000,
                              iter = 8000,
                              control = list(adapt_delta = 0.99),
                              save_pars = save_pars(all = TRUE))


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
  

bias_summary_isometric <- draws_isometric |>
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

bias_plot <- draws_isometric |>
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
  geom_text(data = bias_summary_isometric,
            aes(x=coef, y=mean, label=round(mean,2)), nudge_x=-0.2, size = 2.5) +
  geom_text(data = bias_summary_isometric,
            aes(x=coef, y=.lower, label=round(.lower,2)), nudge_x=-0.2, size = 2.5) +
  geom_text(data = bias_summary_isometric,
            aes(x=coef, y=.upper, label=round(.upper,2)), nudge_x=-0.2, size = 2.5) +
  facet_grid(.~exercise) +
  labs(
    title = "Between-day mean bias",
    x = "Contrast Coefficients for Day (Helmert Coding)",
    y = expression(beta[j]-beta[j-1]),
    color = "Day"
  ) +
  theme_bw()


loamr_summary_isometric <- draws_isometric |>
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
  

loam_plot <- isometric_data |>
  mutate(exercise = case_when(
    exercise == "cp" ~ "Chest Press",
    exercise == "lp" ~ "Leg Press",
    exercise == "row" ~ "Row"
  )) |>
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  geom_ribbon(data = loamr_summary_isometric |> filter(bw_w == "w"),
              aes(x=mean_i, ymin = .lower, ymax = .upper, group=exercise),
              alpha = 0.75, fill = "#E69F00") +
  geom_ribbon(data = loamr_summary_isometric |> filter(bw_w == "w"),
              aes(x=mean_i, ymin = -.lower, ymax = -.upper, group=exercise),
              alpha = 0.75, fill = "#E69F00") +
  geom_line(data = loamr_summary_isometric |> filter(bw_w == "w"),
              aes(x=mean_i, y = mean, group=exercise),
              alpha = 0.75, linetype = "dotted") +
  geom_line(data = loamr_summary_isometric |> filter(bw_w == "w"),
            aes(x=mean_i, y = -mean, group=exercise), linetype = "dotted") +
  
  
  geom_ribbon(data = loamr_summary_isometric |> filter(bw_w == "bw"),
              aes(x=mean_i, ymin = .lower, ymax = .upper, group=exercise),
              alpha = 0.75, fill = "#56B4E9") +
  geom_ribbon(data = loamr_summary_isometric |> filter(bw_w == "bw"),
              aes(x=mean_i, ymin = -.lower, ymax = -.upper, group=exercise),
              alpha = 0.75, fill = "#56B4E9") +
  geom_line(data = loamr_summary_isometric |> filter(bw_w == "bw"),
            aes(x=mean_i, y = mean, group=exercise)) +
  geom_line(data = loamr_summary_isometric |> filter(bw_w == "bw"),
            aes(x=mean_i, y = -mean, group=exercise)) +
  
  # geom_label(data = isometric_data |> mutate(exercise = case_when(
  #   exercise == "cp" ~ "Chest Press",
  #   exercise == "lp" ~ "Leg Press",
  #   exercise == "row" ~ "Row"
  # )) |> group_by(participant, exercise) |> slice_max(mean_diff, n=1),
  #            aes(x = mean_i, y = mean_diff + 10, label = participant)) +
  
  geom_text(data = loamr_summary_isometric |> filter(bw_w == "w") |> group_by(exercise) |> 
              mutate(x = (max(mean_i) - min(mean_i))/2) |> slice_min(mean_i, n=1) |> slice_head(n=1),
            aes(x=mean_i + x, y=-75, group=exercise, 
                label=glue::glue("Within-day: +/- {round(mean,2)} [95% QI: {round(.lower,2)},{round(.upper,2)}]")),
            size = 2.5) +
  geom_text(data = loamr_summary_isometric |> filter(bw_w == "bw") |> group_by(exercise) |> 
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
    caption = "Note, the solid horizontal lines with pale blue ribbons show the between day limits of agreement,\nthe dotted horizontal lines with pale orange ribbons show the within day limits of agreement",
    x = expression(bar(y)[i..]),
    y = expression(y[ijk]-bar(y)[i..]),
    color = "Day"
  ) +
  theme_bw()

bias_plot / loam_plot










### Isometric
isometric_data_cp <- isometric_data |>
  select(participant, exercise, day, trial, value) |>
  filter(exercise == "cp") |>
  group_by(participant) |>
  mutate(mean_i = mean(value),
         mean_diff = value - mean_i)


isometric_model_cp <- brm(value ~ 1 + day + (1|participant/day),
                         data = isometric_data_cp,
                         chains = 4,
                         cores = 4,
                         seed = 1988,
                         warmup = 2000,
                         iter = 8000,
                         control = list(adapt_delta = 0.99),
                         save_pars = save_pars(all = TRUE))

draws_isometric_cp <- spread_draws(isometric_model_cp,
                                 b_Intercept,
                                 b_day1, b_day2,
                                 sd_participant__Intercept,
                                 `sd_participant:day__Intercept`,
                                 sigma) |>
  mutate(
    loam_bw_days = 1.96 * sqrt( ((length(unique(isometric_data_cp$day))-1)/length(unique(isometric_data_cp$day)))*`sd_participant:day__Intercept`^2 + 
                                  ((length(unique(isometric_data_cp$day))*length(unique(isometric_data_cp$trial))-1)/(length(unique(isometric_data_cp$day))*length(unique(isometric_data_cp$trial))))*(sigma^2)),
    
    loam_w_days = 1.96 * sqrt( ((length(unique(isometric_data_cp$trial))*length(unique(isometric_data_cp$day))-1)/(length(unique(isometric_data_cp$trial))*length(unique(isometric_data_cp$day))))*sigma^2)
         )

bias_summary_isometric_cp <- draws_isometric_cp |>
  select(b_day1, b_day2) |>
  pivot_longer(1:2, 
               names_to = "coef",
               values_to = "mean")  |>
  mutate(coef = case_when(
    coef == "b_day1" ~ "Day 1 to Day 2",
    coef == "b_day2" ~ "Day 2 to Day 3"
  )) |>
  group_by(coef) |>
  mean_qi()

bias_plot_cp <- draws_isometric_cp |>
  select(b_day1, b_day2) |>
  pivot_longer(1:2, 
               names_to = "coef",
               values_to = "draw") |>
  mutate(coef = case_when(
    coef == "b_day1" ~ "Day 1 to Day 2",
    coef == "b_day2" ~ "Day 2 to Day 3"
  )) |>
  ggplot(aes(x = coef, y = draw)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_slabinterval(alpha = 0.75, .width = .95) +
  geom_text(data = bias_summary_isometric_cp,
            aes(x=coef, y=mean, label=round(mean,2)), nudge_x=-0.15) +
  geom_text(data = bias_summary_isometric_cp,
            aes(x=coef, y=.lower, label=round(.lower,2)), nudge_x=-0.15) +
  geom_text(data = bias_summary_isometric_cp,
            aes(x=coef, y=.upper, label=round(.upper,2)), nudge_x=-0.15) +
  labs(
    title = "Chest Press Exercise",
    subtitle = "Between-day mean bias",
    x = "Contrast Coefficients for Day (Helmert Coding)",
    y = expression(beta[j]-beta[j-1]),
    color = "Day"
  ) +
  theme_bw()

  
loamr_summary_isometric_cp <- draws_isometric_cp |>
  mean_qi(loam_bw_days, loam_w_days)

loam_plot_cp <- isometric_data_cp |>
  ggplot(aes(x=mean_i, y=mean_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = loamr_summary_isometric_cp$loam_w_days.lower, ymax = loamr_summary_isometric_cp$loam_w_days.upper,
           alpha = 0.75, fill = "#E69F00") +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = -loamr_summary_isometric_cp$loam_w_days.lower, ymax = -loamr_summary_isometric_cp$loam_w_days.upper,
           alpha = 0.75, fill = "#E69F00") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = loamr_summary_isometric_cp$loam_bw_days.lower, ymax = loamr_summary_isometric_cp$loam_bw_days.upper,
           alpha = 0.5, fill = "#56B4E9") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = -loamr_summary_isometric_cp$loam_bw_days.lower, ymax = -loamr_summary_isometric_cp$loam_bw_days.upper,
           alpha = 0.5, fill = "#56B4E9") +

  geom_hline(yintercept = c(-loamr_summary_isometric_cp$loam_w_days,loamr_summary_isometric_cp$loam_w_days), linetype = "dotted") +
  geom_hline(yintercept = c(-loamr_summary_isometric_cp$loam_bw_days,loamr_summary_isometric_cp$loam_bw_days)) +
  geom_label(data = isometric_data_cp |> group_by(participant) |> slice_max(mean_diff, n=1),
             aes(x = mean_i, y = mean_diff + 10, label = participant)) +
  geom_point(size = 2, color = "black", fill = NA, shape = 21, alpha = 0.75) +
  geom_point(aes(color=day), alpha = 0.75) +
  scale_color_viridis_d() +
  scale_y_continuous(breaks = seq(-100, 100, by=25)) +
  labs(
    title = "Chest Press Exercise",
    subtitle = glue::glue("Within-day trial to trial limits of agreement: +/- {round(loamr_summary_isometric_cp$loam_w_days,2)} [95% QI: {round(loamr_summary_isometric_cp$loam_w_days.lower,2)},{round(loamr_summary_isometric_cp$loam_w_days.upper,2)}]\nBetween-day trial to trial limits of agreement: +/- {round(loamr_summary_isometric_cp$loam_bw_days,2)} [95% QI: {round(loamr_summary_isometric_cp$loam_bw_days.lower,2)},{round(loamr_summary_isometric_cp$loam_bw_days.upper,2)}]"),
    x = expression(bar(y)[i..]),
    y = expression(y[ijk]-bar(y)[i..]),
    color = "Day"
  ) +
  theme_bw()
  

isometric_data_lp <- isometric_data |>
  select(participant, exercise, day, trial, value) |>
  filter(exercise == "lp") |>
  group_by(participant) |>
  mutate(mean_i = mean(value),
         mean_diff = value - mean_i)


isometric_model_lp <- brm(value ~ 1 + day + (1|participant/day),
                              data = isometric_data_lp,
                              chains = 4,
                              cores = 4,
                              seed = 1988,
                              warmup = 2000,
                              iter = 8000,
                              control = list(adapt_delta = 0.99),
                              save_pars = save_pars(all = TRUE))

draws_isometric_lp <- spread_draws(isometric_model_lp,
                                   b_Intercept,
                                   b_day1, b_day2,
                                   sd_participant__Intercept,
                                   `sd_participant:day__Intercept`,
                                   sigma) |>
  mutate(
    loam_bw_days = 1.96 * sqrt( ((length(unique(isometric_data_lp$day))-1)/length(unique(isometric_data_lp$day)))*`sd_participant:day__Intercept`^2 + 
                                  ((length(unique(isometric_data_lp$day))*length(unique(isometric_data_lp$trial))-1)/(length(unique(isometric_data_lp$day))*length(unique(isometric_data_lp$trial))))*(sigma^2)),
    
    loam_w_days = 1.96 * sqrt( ((length(unique(isometric_data_lp$trial))*length(unique(isometric_data_lp$day))-1)/(length(unique(isometric_data_lp$trial))*length(unique(isometric_data_lp$day))))*sigma^2)
  )

bias_summary_isometric_lp <- draws_isometric_lp |>
  select(b_day1, b_day2) |>
  pivot_longer(1:2, 
               names_to = "coef",
               values_to = "mean")  |>
  mutate(coef = case_when(
    coef == "b_day1" ~ "Day 1 to Day 2",
    coef == "b_day2" ~ "Day 2 to Day 3"
  )) |>
  group_by(coef) |>
  mean_qi()

bias_plot_lp <- draws_isometric_lp |>
  select(b_day1, b_day2) |>
  pivot_longer(1:2, 
               names_to = "coef",
               values_to = "draw") |>
  mutate(coef = case_when(
    coef == "b_day1" ~ "Day 1 to Day 2",
    coef == "b_day2" ~ "Day 2 to Day 3"
  )) |>
  ggplot(aes(x = coef, y = draw)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_slabinterval(alpha = 0.75, .width = .95) +
  geom_text(data = bias_summary_isometric_lp,
            aes(x=coef, y=mean, label=round(mean,2)), nudge_x=-0.15) +
  geom_text(data = bias_summary_isometric_lp,
            aes(x=coef, y=.lower, label=round(.lower,2)), nudge_x=-0.15) +
  geom_text(data = bias_summary_isometric_lp,
            aes(x=coef, y=.upper, label=round(.upper,2)), nudge_x=-0.15) +
  labs(
    title = "Leg Press Exercise",
    subtitle = "Between-day mean bias",
    x = "Contrast Coefficients for Day (Helmert Coding)",
    y = expression(beta[j]-beta[j-1]),
    color = "Day"
  ) +
  theme_bw()

loamr_summary_isometric_lp <- draws_isometric_lp |>
  mean_qi(loam_bw_days, loam_w_days)

loam_plot_lp <- isometric_data_lp |>
  ggplot(aes(x=mean_i, y=mean_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = loamr_summary_isometric_lp$loam_w_days.lower, ymax = loamr_summary_isometric_lp$loam_w_days.upper,
           alpha = 0.75, fill = "#E69F00") +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = -loamr_summary_isometric_lp$loam_w_days.lower, ymax = -loamr_summary_isometric_lp$loam_w_days.upper,
           alpha = 0.75, fill = "#E69F00") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = loamr_summary_isometric_lp$loam_bw_days.lower, ymax = loamr_summary_isometric_lp$loam_bw_days.upper,
           alpha = 0.5, fill = "#56B4E9") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = -loamr_summary_isometric_lp$loam_bw_days.lower, ymax = -loamr_summary_isometric_lp$loam_bw_days.upper,
           alpha = 0.5, fill = "#56B4E9") +
  
  geom_hline(yintercept = c(-loamr_summary_isometric_lp$loam_w_days,loamr_summary_isometric_lp$loam_w_days), linetype = "dotted") +
  geom_hline(yintercept = c(-loamr_summary_isometric_lp$loam_bw_days,loamr_summary_isometric_lp$loam_bw_days)) +
  geom_label(data = isometric_data_lp |> group_by(participant) |> slice_max(mean_diff, n=1),
             aes(x = mean_i, y = mean_diff + 10, label = participant)) +
  geom_point(size = 2, color = "black", fill = NA, shape = 21, alpha = 0.75) +
  geom_point(aes(color=day), alpha = 0.75) +
  scale_color_viridis_d() +
  scale_y_continuous(breaks = seq(-100, 100, by=25)) +
  labs(
    title = "Leg Press Exercise",
    subtitle = glue::glue("Within-day trial to trial limits of agreement: +/- {round(loamr_summary_isometric_lp$loam_w_days,2)} [95% QI: {round(loamr_summary_isometric_lp$loam_w_days.lower,2)},{round(loamr_summary_isometric_lp$loam_w_days.upper,2)}]\nBetween-day trial to trial limits of agreement: +/- {round(loamr_summary_isometric_lp$loam_bw_days,2)} [95% QI: {round(loamr_summary_isometric_lp$loam_bw_days.lower,2)},{round(loamr_summary_isometric_lp$loam_bw_days.upper,2)}]"),
    x = expression(bar(y)[i..]),
    y = expression(y[ijk]-bar(y)[i..]),
    color = "Day"
  ) +
  theme_bw()

isometric_data_row <- isometric_data |>
  select(participant, exercise, day, trial, value) |>
  filter(exercise == "row") |>
  group_by(participant) |>
  mutate(mean_i = mean(value),
         mean_diff = value - mean_i)


isometric_model_row <- brm(value ~ 1 + day + (1|participant/day),
                              data = isometric_data_row,
                              chains = 4,
                              cores = 4,
                              seed = 1988,
                              warmup = 2000,
                              iter = 8000,
                              control = list(adapt_delta = 0.99),
                              save_pars = save_pars(all = TRUE))

draws_isometric_row <- spread_draws(isometric_model_row,
                                   b_Intercept,
                                   b_day1, b_day2,
                                   sd_participant__Intercept,
                                   `sd_participant:day__Intercept`,
                                   sigma) |>
  mutate(
    loam_bw_days = 1.96 * sqrt( ((length(unique(isometric_data_row$day))-1)/length(unique(isometric_data_row$day)))*`sd_participant:day__Intercept`^2 + 
                                  ((length(unique(isometric_data_row$day))*length(unique(isometric_data_row$trial))-1)/(length(unique(isometric_data_row$day))*length(unique(isometric_data_row$trial))))*(sigma^2)),
    
    loam_w_days = 1.96 * sqrt( ((length(unique(isometric_data_row$trial))*length(unique(isometric_data_row$day))-1)/(length(unique(isometric_data_row$trial))*length(unique(isometric_data_row$day))))*sigma^2)
  )

bias_summary_isometric_row <- draws_isometric_row |>
  select(b_day1, b_day2) |>
  pivot_longer(1:2, 
               names_to = "coef",
               values_to = "mean")  |>
  mutate(coef = case_when(
    coef == "b_day1" ~ "Day 1 to Day 2",
    coef == "b_day2" ~ "Day 2 to Day 3"
  )) |>
  group_by(coef) |>
  mean_qi()

bias_plot_row <- draws_isometric_row |>
  select(b_day1, b_day2) |>
  pivot_longer(1:2, 
               names_to = "coef",
               values_to = "draw") |>
  mutate(coef = case_when(
    coef == "b_day1" ~ "Day 1 to Day 2",
    coef == "b_day2" ~ "Day 2 to Day 3"
  )) |>
  ggplot(aes(x = coef, y = draw)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_slabinterval(alpha = 0.75, .width = .95) +
  geom_text(data = bias_summary_isometric_row,
            aes(x=coef, y=mean, label=round(mean,2)), nudge_x=-0.15) +
  geom_text(data = bias_summary_isometric_row,
            aes(x=coef, y=.lower, label=round(.lower,2)), nudge_x=-0.15) +
  geom_text(data = bias_summary_isometric_row,
            aes(x=coef, y=.upper, label=round(.upper,2)), nudge_x=-0.15) +
  labs(
    title = "Row Exercise",
    subtitle = "Between-day mean bias",
    x = "Contrast Coefficients for Day (Helmert Coding)",
    y = expression(beta[j]-beta[j-1]),
    color = "Day"
  ) +
  theme_bw()

loamr_summary_isometric_row <- draws_isometric_row |>
  mean_qi(loam_bw_days, loam_w_days)

loam_plot_row <- isometric_data_row |>
  ggplot(aes(x=mean_i, y=mean_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = loamr_summary_isometric_row$loam_w_days.lower, ymax = loamr_summary_isometric_row$loam_w_days.upper,
           alpha = 0.75, fill = "#E69F00") +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = -loamr_summary_isometric_row$loam_w_days.lower, ymax = -loamr_summary_isometric_row$loam_w_days.upper,
           alpha = 0.75, fill = "#E69F00") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = loamr_summary_isometric_row$loam_bw_days.lower, ymax = loamr_summary_isometric_row$loam_bw_days.upper,
           alpha = 0.5, fill = "#56B4E9") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = -loamr_summary_isometric_row$loam_bw_days.lower, ymax = -loamr_summary_isometric_row$loam_bw_days.upper,
           alpha = 0.5, fill = "#56B4E9") +
  
  geom_hline(yintercept = c(-loamr_summary_isometric_row$loam_w_days,loamr_summary_isometric_row$loam_w_days), linetype = "dotted") +
  geom_hline(yintercept = c(-loamr_summary_isometric_row$loam_bw_days,loamr_summary_isometric_row$loam_bw_days)) +
  geom_label(data = isometric_data_row |> group_by(participant) |> slice_max(mean_diff, n=1),
             aes(x = mean_i, y = mean_diff + 10, label = participant)) +
  geom_point(size = 2, color = "black", fill = NA, shape = 21, alpha = 0.75) +
  geom_point(aes(color=day), alpha = 0.75) +
  scale_color_viridis_d() +
  scale_y_continuous(breaks = seq(-100, 100, by=25)) +
  labs(
    title = "Row Exercise",
    subtitle = glue::glue("Within-day trial to trial limits of agreement: +/- {round(loamr_summary_isometric_row$loam_w_days,2)} [95% QI: {round(loamr_summary_isometric_row$loam_w_days.lower,2)},{round(loamr_summary_isometric_row$loam_w_days.upper,2)}]\nBetween-day trial to trial limits of agreement: +/- {round(loamr_summary_isometric_row$loam_bw_days,2)} [95% QI: {round(loamr_summary_isometric_row$loam_bw_days.lower,2)},{round(loamr_summary_isometric_row$loam_bw_days.upper,2)}]"),
    x = expression(bar(y)[i..]),
    y = expression(y[ijk]-bar(y)[i..]),
    color = "Day"
  ) +
  theme_bw()

((bias_plot_cp / bias_plot_lp / bias_plot_row) |
(loam_plot_cp / loam_plot_lp / loam_plot_row)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Isometric Tests",
                  caption = "Note, the solid horizontal lines with pale blue ribbons show the between day limits of agreement,\nthe dotted horizontal lines with pale orange ribbons show the within day limits of agreement,\nnumeric labels are participant numbers")


### Isokinetic
isokinetic_data_cp <- isokinetic_data |>
  select(participant, exercise, day, con_ecc, value) |>
  filter(exercise == "cp") |>
  group_by(participant, con_ecc) |>
  mutate(mean_i = mean(value),
         mean_diff = value - mean_i)

# isokinetic_data_cp_wide <- isokinetic_data_cp |>
#   pivot_wider(names_from = c("day", "con_ecc"),
#               values_from = "value",
#               id_cols = c("participant", "exercise")) |>
  

con_bf <- bf(con ~ 1 + day + (1|i|participant))
ecc_bf <- bf(ecc ~ 1 + day + (1|i|participant))

brm_model_isokinetic_cp <- brm(con_bf + ecc_bf,
                              data = isokinetic_data_cp_wide,
                              chains = 4,
                              cores = 4,
                              seed = 1988,
                              warmup = 2000,
                              iter = 8000,
                              control = list(adapt_delta = 0.99),
                              save_pars = save_pars(all = TRUE))

draws_isokinetic_cp <- spread_draws(brm_model_isokinetic_cp,
                                    sigma_con,
                                    sigma_ecc
                                   ) |>
  mutate(
    loam_bw_days_con = 1.96 * sigma_con,
    loam_bw_days_ecc = 1.96 * sigma_ecc
  ) 

loamr_summary_isokinetic_cp <- draws_isokinetic_cp |>
  mean_qi(loam_bw_days_con, loam_bw_days_ecc) 

loam_plot_cp_con <- isokinetic_data_cp |>
  filter(con_ecc == "con") |>
  ggplot(aes(x=mean_i, y=mean_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = loamr_summary_isokinetic_cp$loam_bw_days_con.lower, ymax = loamr_summary_isokinetic_cp$loam_bw_days_con.upper,
           alpha = 0.5, fill = "#56B4E9") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = -loamr_summary_isokinetic_cp$loam_bw_days_con.lower, ymax = -loamr_summary_isokinetic_cp$loam_bw_days_con.upper,
           alpha = 0.5, fill = "#56B4E9") +
  
  geom_hline(yintercept = c(-loamr_summary_isokinetic_cp$loam_bw_days_con,loamr_summary_isokinetic_cp$loam_bw_days_con)) +
  
  geom_label(data = isokinetic_data_cp |> filter(con_ecc == "con") |> group_by(participant) |> slice_max(mean_diff, n=1),
             aes(x = mean_i, y = mean_diff + 10, label = participant)) +
  geom_point(size = 2, color = "black", fill = NA, shape = 21, alpha = 0.75) +
  geom_point(aes(color=day), alpha = 0.75) +
  scale_color_viridis_d() +
  scale_y_continuous(breaks = seq(-100, 100, by=10)) +
  labs(
    title = "Chest Press Exercise - Concentric",
    subtitle = glue::glue("Between-day trial to trial limits of agreement: +/- {round(loamr_summary_isokinetic_cp$loam_bw_days_con,2)} [95% QI: {round(loamr_summary_isokinetic_cp$loam_bw_days_con.lower,2)},{round(loamr_summary_isokinetic_cp$loam_bw_days_con.upper,2)}]"),
    x = "Mean Force (N)",
    y = "Difference from the Participant Mean (N)",
    color = "Day"
  ) +
  theme_bw()

loam_plot_cp_ecc <- isokinetic_data_cp |>
  filter(con_ecc == "ecc") |>
  ggplot(aes(x=mean_i, y=mean_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = loamr_summary_isokinetic_cp$loam_bw_days_ecc.lower, ymax = loamr_summary_isokinetic_cp$loam_bw_days_ecc.upper,
           alpha = 0.5, fill = "#56B4E9") +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = -loamr_summary_isokinetic_cp$loam_bw_days_ecc.lower, ymax = -loamr_summary_isokinetic_cp$loam_bw_days_ecc.upper,
           alpha = 0.5, fill = "#56B4E9") +
  
  geom_hline(yintercept = c(-loamr_summary_isokinetic_cp$loam_bw_days_ecc,loamr_summary_isokinetic_cp$loam_bw_days_ecc)) +
  
  geom_label(data = isokinetic_data_cp |> filter(con_ecc == "ecc") |> group_by(participant) |> slice_max(mean_diff, n=1),
             aes(x = mean_i, y = mean_diff + 10, label = participant)) +
  geom_point(size = 2, color = "black", fill = NA, shape = 21, alpha = 0.75) +
  geom_point(aes(color=day), alpha = 0.75) +
  scale_color_viridis_d() +
  scale_y_continuous(breaks = seq(-100, 100, by=10)) +
  labs(
    title = "Chest Press Exercise - Eccentric",
    subtitle = glue::glue("Between-day trial to trial limits of agreement: +/- {round(loamr_summary_isokinetic_cp$loam_bw_days_ecc,2)} [95% QI: {round(loamr_summary_isokinetic_cp$loam_bw_days_ecc.lower,2)},{round(loamr_summary_isokinetic_cp$loam_bw_days_ecc.upper,2)}]"),
    x = "Mean Force (N)",
    y = "Difference from the Participant Mean (N)",
    color = "Day"
  ) +
  theme_bw()

(loam_plot_cp_con | loam_plot_cp_ecc ) +
  # (loam_plot_lp_con | loam_plot_lp_ecc ) +
  # (loam_plot_row_con | loam_plot_row_ecc ) 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Limits of Agreement with the Mean",
                  caption = "Note, the solid horizontal lines with pale blue ribbons show the between day limits of agreement,\nnumeric labels are participant numbers")

