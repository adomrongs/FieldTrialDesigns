# Script to plot the results of the simulation from the combined output of each
# scenario (outputs/<Scenario>_output.csv). Figures are saved in outputs/plots

# ------- Prepare Data -------
# Load packages
rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)

# Define scenarios and colours (the same colour for each scenario in every plot)
scenarios <- c("PS_RCBD", "GS_RCBD", "GS_Augmented", "GS_pREP_odw")
colors_scenario <- c(
  PS_RCBD = "#335c67",
  GS_RCBD = "#2d6a4f",
  GS_Augmented = "#fff3b0",
  GS_pREP_odw = "#540b0e"
)

# colors_scenario <- c(
#   BurnIn = "#335c67",
#   PS_RCBD = "#fff3b0",
#   GS_RCBD = "#e09f3e",
#   GS_Augmented = "#9e2a2b",
#   GS_pREP_odw = "#540b0e"
# )

# Cost per plot (phenotyping) and per sample (genotyping), EUR
cost_plot_S2 <- 100
cost_plot_S3 <- 250
cost_geno <- 40

# Annual cost of what differs between the programs: S2 (phenotyped in PS,
# genotyped in GS) + S3 plots (design). S1, S4 and S5 are the same in all of them
plots_S3 <- c(PS_RCBD = 600, GS_RCBD = 600, GS_Augmented = 240, GS_pREP_odw = 250)
cost_S2 <- c(
  PS_RCBD = 1000 * cost_plot_S2,
  GS_RCBD = 1000 * cost_geno,
  GS_Augmented = 1000 * cost_geno,
  GS_pREP_odw = 1000 * cost_geno
)
cost_scenario <- cost_S2 + plots_S3 * cost_plot_S3
cost_scenario

# Read outputs
burnin <- read.csv("outputs/BurnIn_output.csv")
future <- bind_rows(
  lapply(file.path("outputs", paste0(scenarios, "_output.csv")), read.csv)
)
last_burnin <- max(burnin$year)

# Genetic gain w.r.t. the end of the burn-in of each replicate, and gain per
# EUR of the annual cost of the program
gv_start <- burnin |>
  filter(year == last_burnin) |>
  select(rep, gv_start = mean_gv_tpe)

results <- bind_rows(burnin, future) |>
  left_join(gv_start, by = "rep") |>
  mutate(
    gain = mean_gv_tpe - gv_start,
    gain_cost = gain / cost_scenario[scenario],
    scenario = factor(scenario, levels = names(colors_scenario))
  ) |>
  select(-gv_start)

# Mean and 95 % CI (t distribution, as t.test) of every variable across the
# replicates, per scenario and year
summary_results <- results |>
  pivot_longer(-c(year, rep, scenario), names_to = "variable") |>
  filter(!is.na(value)) |>
  summarise(
    n = n(),
    mean = mean(value),
    se = sd(value) / sqrt(n),
    .by = c(scenario, year, variable)
  ) |>
  mutate(
    lower = mean - qt(0.975, n - 1) * se,
    upper = mean + qt(0.975, n - 1) * se
  )

# Names of the variables shown in the facet strips
labels_variable <- c(
  gain = "Genetic gain (TPE)",
  gain_cost = "Genetic gain per EUR of annual cost (TPE)",
  acc_S3_met = "Selection accuracy in S3 (trial environment)",
  acc_S3_tpe = "Selection accuracy in S3 (TPE)",
  acc_gs = "Selection accuracy in S2 (TPE)"
)

dir.create("outputs/plots", showWarnings = FALSE, recursive = TRUE)

# ------- Genetic gain -------
# Future years, starting from 0 at the end of the burn-in
plot_gain <- summary_results |>
  filter(variable == "gain", scenario %in% scenarios) |>
  bind_rows(data.frame(
    scenario = factor(scenarios, levels = names(colors_scenario)),
    variable = "gain",
    year = last_burnin,
    mean = 0,
    lower = 0,
    upper = 0
  )) |>
  ggplot(aes(year, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = colors_scenario) +
  scale_fill_manual(values = colors_scenario) +
  facet_wrap(~variable, labeller = as_labeller(labels_variable)) +
  labs(x = "Year", y = NULL, colour = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")
plot_gain
ggsave(
  "outputs/plots/genetic_gain.png",
  plot_gain,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ------- Genetic gain per cost -------
# Gain since the end of the burn-in divided by the annual cost of the program
plot_gain_cost <- summary_results |>
  filter(variable == "gain_cost") |>
  ggplot(aes(year, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = colors_scenario) +
  scale_fill_manual(values = colors_scenario) +
  facet_wrap(~variable, labeller = as_labeller(labels_variable)) +
  labs(x = "Year", y = NULL, colour = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")
plot_gain_cost
ggsave(
  "outputs/plots/genetic_gain_cost.png",
  plot_gain_cost,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ------- S3 selection accuracy -------
# Accuracy of the values used to select S4 from S3 (raw means in PS, GEBVs in
# GS) w.r.t. the genetic value in the S3 trial environment and in the TPE
plot_acc_S3 <- summary_results |>
  filter(variable %in% c("acc_S3_met", "acc_S3_tpe"), scenario %in% scenarios) |>
  ggplot(aes(year, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~variable, labeller = as_labeller(labels_variable)) +
  scale_colour_manual(values = colors_scenario) +
  scale_fill_manual(values = colors_scenario) +
  labs(x = "Year", y = NULL, colour = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")
plot_acc_S3
ggsave(
  "outputs/plots/accuracy_S3.png",
  plot_acc_S3,
  width = 10,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ------- S2 selection accuracy -------
# Accuracy of the GEBVs used to select S3 from S2 (GS scenarios only; the
# training population includes the S3 trials of each design)
plot_acc_S2 <- summary_results |>
  filter(variable == "acc_gs", scenario %in% scenarios) |>
  ggplot(aes(year, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~variable, labeller = as_labeller(labels_variable)) +
  scale_colour_manual(values = colors_scenario) +
  scale_fill_manual(values = colors_scenario) +
  labs(x = "Year", y = NULL, colour = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")
plot_acc_S2
ggsave(
  "outputs/plots/accuracy_S2.png",
  plot_acc_S2,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)
