# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) 

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse",
    "here",
    "janitor",
    "patchwork",
    "rstan",
    "brms",
    "tidybayes",
    "bayesplot",
    "grateful",
    "quarto",
    "cli"
  ) # Packages that your targets need for their tasks.
  
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  # 
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("R/functions.R")
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  
  # Read in and prepare data
  tar_target(file, here("data","data.csv"), format = "file"),
  tar_target(data, read_csv(file) |> clean_names()),
  tar_target(isometric_data, prepare_isometric_data(data)),
  tar_target(isokinetic_data, prepare_isokinetic_data(data)),
  
  # Fit and plot isometric model
  tar_target(isometric_model, fit_isometric_model(isometric_data)),
  tar_target(isometric_model_rhat, make_rhat_plot(isometric_model)),
  tar_target(isometric_model_trace, make_trace_plot(isometric_model)),
  tar_target(isometric_model_pp_check_cp, make_pp_check(isometric_model, "cp")),
  tar_target(isometric_model_pp_check_lp, make_pp_check(isometric_model, "lp")),
  tar_target(isometric_model_pp_check_row, make_pp_check(isometric_model, "row")),
  tar_target(isometric_draws, get_isometric_draws(isometric_model, isometric_data)),
  tar_target(isometric_bias_summary, get_isometric_bias_summary(isometric_draws)),
  tar_target(isometric_loamr_summary, get_isometric_loamr_summary(isometric_data,isometric_draws)),
  tar_target(isometric_bias_plot, plot_isometric_bias(isometric_draws, isometric_bias_summary)),
  tar_target(isometric_loamr_plot, plot_isometric_loamr(isometric_data, isometric_loamr_summary)),
  tar_target(combined_isometric_plot, (isometric_bias_plot / isometric_loamr_plot) +
               plot_annotation(caption = "Note, the solid horizontal lines (point estimates) with pale blue ribbons (95% quantile intervals) show the between day limits of agreement,\nthe dotted horizontal lines (point estimates) with pale orange ribbons (95% quantile intervals) show the within day limits of agreement") &
               theme(plot.caption = element_text(size=6))),
  tar_target(isometric_plot_tiff, make_plot_tiff(combined_isometric_plot, 10, 6.66, "plots/isometric_plot.tiff")),
  
  # Fit and plot isokinetic model
  tar_target(isokinetic_model, fit_isokinetic_model(isokinetic_data)),
  tar_target(isokinetic_model_rhat, make_rhat_plot(isokinetic_model)),
  tar_target(isokinetic_model_trace, make_trace_plot(isokinetic_model)),
  tar_target(isokinetic_model_pp_check_cp_con, make_pp_check(isokinetic_model, "deltacpcon")),
  tar_target(isokinetic_model_pp_check_lp_con, make_pp_check(isokinetic_model, "deltalpcon")),
  tar_target(isokinetic_model_pp_check_row_con, make_pp_check(isokinetic_model, "deltarowcon")),
  tar_target(isokinetic_model_pp_check_cp_ecc, make_pp_check(isokinetic_model, "deltacpecc")),
  tar_target(isokinetic_model_pp_check_lp_ecc, make_pp_check(isokinetic_model, "deltalpecc")),
  tar_target(isokinetic_model_pp_check_row_ecc, make_pp_check(isokinetic_model, "deltarowecc")),
  tar_target(isokinetic_draws, get_isokinetic_draws(isokinetic_model)),
  tar_target(isokinetic_bias_loa_summary, get_isokinetic_bias_loa_summary(isokinetic_draws, isokinetic_data)),
  tar_target(isokinetic_BA_plot, plot_isokinetic_BA(isokinetic_bias_loa_summary, isokinetic_data)),
  tar_target(isokinetic_BA_plot_tiff, make_plot_tiff(isokinetic_BA_plot, 10, 6.66, "plots/isokinetic_BA_plot.tiff")),
  
  # Reporting
  tar_target(grateful_report, cite_packages(out.dir = ".", cite.tidyverse = TRUE, out.format = "pdf"))
  # tar_quarto(isometric_model_diagnostics, path = "plots/isometric_model_diagnostics.qmd"),
  # tar_quarto(isokinetic_model_diagnostics, path = "plots/isokinetic_model_diagnostics.qmd")
  # tar_quarto(analysis_results, path = "analysis_results.qmd")
  
  )








