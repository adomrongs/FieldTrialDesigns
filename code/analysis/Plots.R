# Script to plot the results of the simulation from the combined output of each
# scenario (outputs/<Scenario>_output.csv). Figures are saved in outputs/plots

# ------- Prepare Data -------
# Load packages
rm(list = ls())
library(tidyverse)

# Define scenarios, colours and line types (the same for each scenario in every plot).
# Every scenario has the same budget: S2 50,000 EUR + S3 and S4 240,000 EUR per year
#   controls (red): PS_RCBD dashed (phenotypic selection), GS_RCBD solid
#   extremes (blues, dark to light): ES1, ES2, ES2_ST
#   intermediates (greens, light to dark): I1 to I4
#   ES2_ST_G (purple): ES2 entries in sparse testing + the savings to genotype more S2
scenarios <- c(
  "PS_RCBD",
  "GS_RCBD",
  "ES1",
  "ES2",
  "ES2_ST",
  "I1",
  "I2",
  "I3",
  "I4",
  "ES2_ST_G"
)
colors_scenario <- c(
  PS_RCBD = "#DD5129",
  GS_RCBD = "#DD5129",
  ES1 = "#0B2E59",
  ES2 = "#4A86D6",
  ES2_ST = "#8AB6F2",
  I1 = "#A8DDA0",
  I2 = "#6BBF63",
  I3 = "#2E8B3E",
  I4 = "#0B3D1E",
  ES2_ST_G = "#7B3294"
)
linetypes_scenario <- c(
  PS_RCBD = "dashed",
  GS_RCBD = "solid",
  ES1 = "solid",
  ES2 = "solid",
  ES2_ST = "solid",
  I1 = "solid",
  I2 = "solid",
  I3 = "solid",
  I4 = "solid",
  ES2_ST_G = "solid"
)

# Read outputs
burnin <- read.csv("outputs/BurnIn_output.csv")
future_files <- file.path("outputs", paste0(scenarios, "_output.csv"))
future <- bind_rows(
  lapply(future_files[file.exists(future_files)], read.csv) # scenarios already run
)
last_burnin <- max(burnin$year)

# Genetic gain in S1 and S5 w.r.t. the end of the burn-in of each replicate
gv_start <- burnin |>
  filter(year == last_burnin) |>
  select(rep, start_S1 = mean_gv_S1, start_S5 = mean_gv_S5)

results <- future |>
  left_join(gv_start, by = "rep") |>
  mutate(
    gain_S1 = mean_gv_S1 - start_S1,
    gain_S5 = mean_gv_S5 - start_S5,
    scenario = factor(scenario, levels = scenarios)
  ) |>
  select(year, rep, scenario, gain_S1, gain_S5, acc_S2, acc_S3, acc_S4)

# Mean and SD of every variable across the replicates, per scenario and year.
# The SD of the genetic gain measures its stability (as in Werner et al. 2025)
summary_results <- results |>
  pivot_longer(-c(year, rep, scenario), names_to = "variable") |>
  filter(!is.na(value)) |>
  summarise(mean = mean(value), sd = sd(value), .by = c(scenario, year, variable))

# Names of the variables shown in the facet strips
labels_variable <- c(
  gain_S1 = "Genetic gain in S1 (TPE)",
  gain_S5 = "Genetic gain in S5 (TPE)",
  acc_S2 = "Selection accuracy in S2 (TPE)",
  acc_S3 = "Selection accuracy in S3 (TPE)",
  acc_S4 = "Selection accuracy in S4 (TPE)"
)

dir.create("outputs/plots", showWarnings = FALSE, recursive = TRUE)

# Groups of scenarios plotted together (the controls in every group)
controls <- c("PS_RCBD", "GS_RCBD")
groups <- list(
  extremes = c(controls, "ES1", "ES2", "ES2_ST", "ES2_ST_G"),
  intermediates = c(controls, "I1", "I2", "I3", "I4")
)

for (group in names(groups)) {
  scenarios_group <- groups[[group]]

  # ------- Genetic gain -------
  # Future years, starting from 0 at the end of the burn-in
  p1 <- summary_results |>
    filter(variable %in% c("gain_S1", "gain_S5"), scenario %in% scenarios_group) |>
    bind_rows(expand.grid(
      scenario = factor(scenarios_group, levels = scenarios),
      variable = c("gain_S1", "gain_S5"),
      year = last_burnin,
      mean = 0,
      stringsAsFactors = FALSE
    )) |>
    ggplot(aes(year, mean, colour = scenario, linetype = scenario)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = colors_scenario) +
    scale_linetype_manual(values = linetypes_scenario) +
    facet_wrap(~variable, scales = "free_y", labeller = as_labeller(labels_variable)) +
    labs(x = "Year", y = NULL, colour = NULL, linetype = NULL) +
    theme_bw() +
    theme(legend.position = "top")

  ggsave(
    file.path("outputs/plots", paste0("genetic_gain_", group, ".png")),
    p1,
    width = 10,
    height = 5,
    dpi = 300,
    bg = "white"
  )

  # ------- Selection accuracy -------
  p2 <- summary_results |>
    filter(
      variable %in% c("acc_S2", "acc_S3", "acc_S4"),
      scenario %in% scenarios_group
    ) |>
    ggplot(aes(year, mean, colour = scenario, linetype = scenario)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = colors_scenario) +
    scale_linetype_manual(values = linetypes_scenario) +
    facet_wrap(~variable, labeller = as_labeller(labels_variable)) +
    labs(x = "Year", y = NULL, colour = NULL, linetype = NULL) +
    theme_bw() +
    theme(legend.position = "top")

  ggsave(
    file.path("outputs/plots", paste0("selection_accuracy_", group, ".png")),
    p2,
    width = 15,
    height = 5,
    dpi = 300,
    bg = "white"
  )

  # ------- SD of the genetic gain -------
  # SD across the replicates (stability of the gain), starting from 0 at the
  # end of the burn-in
  p3 <- summary_results |>
    filter(variable %in% c("gain_S1", "gain_S5"), scenario %in% scenarios_group) |>
    bind_rows(expand.grid(
      scenario = factor(scenarios_group, levels = scenarios),
      variable = c("gain_S1", "gain_S5"),
      year = last_burnin,
      sd = 0,
      stringsAsFactors = FALSE
    )) |>
    ggplot(aes(year, sd, colour = scenario, linetype = scenario)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = colors_scenario) +
    scale_linetype_manual(values = linetypes_scenario) +
    facet_wrap(~variable, scales = "free_y", labeller = as_labeller(labels_variable)) +
    labs(x = "Year", y = "SD across replicates", colour = NULL, linetype = NULL) +
    theme_bw() +
    theme(legend.position = "top")

  ggsave(
    file.path("outputs/plots", paste0("sd_genetic_gain_", group, ".png")),
    p3,
    width = 10,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}

# SD of the genetic gain in the last year, per scenario
summary_results |>
  filter(variable %in% c("gain_S1", "gain_S5"), year == max(year)) |>
  select(scenario, variable, mean, sd) |>
  pivot_wider(names_from = variable, values_from = c(mean, sd)) |>
  arrange(scenario) |>
  print()
