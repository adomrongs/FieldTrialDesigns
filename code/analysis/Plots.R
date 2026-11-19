# Script to plot the results of the simulation from the combined output of each
# scenario (outputs/<Scenario>_output.csv). Figures are saved in outputs/plots

# ------- Prepare Data -------
# Load packages
rm(list = ls())
library(patchwork)
library(tidyverse)

# Define scenarios and colours (the same colour for each scenario in every plot).
# Every scenario has the same budget: S2 50,000 EUR + S3 and S4 240,000 EUR per year
scenarios <- c(
  "PS_RCBD",
  "GS_RCBD",
  "ES2",
  "ES2_ST",
  "I1",
  "I2",
  "I3",
  "I4",
  "ES1",
  "I3_FixedRandom",
  "I3_FixedGRM"
)
colors_scenario <- c(
  PS_RCBD = "#335c67",
  GS_RCBD = "#2d6a4f",
  ES2 = "#7f2704",
  ES2_ST = "#a63603",
  I1 = "#d94801",
  I2 = "#f16913",
  I3 = "#fd8d3c",
  I4 = "#6a51a3",
  ES1 = "#3f007d",
  I3_FixedRandom = "#bdbdbd",
  I3_FixedGRM = "#636363"
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

results <- bind_rows(burnin, future) |>
  left_join(gv_start, by = "rep") |>
  mutate(
    gain_S1 = mean_gv_S1 - start_S1,
    gain_S5 = mean_gv_S5 - start_S5,
    scenario = factor(scenario, levels = c("BurnIn", scenarios))
  ) |>
  select(year, rep, scenario, gain_S1, gain_S5, acc_S2, acc_S3, acc_S4)

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
  gain_S1 = "Genetic gain in S1 (TPE)",
  gain_S5 = "Genetic gain in S5 (TPE)",
  acc_S2 = "Selection accuracy in S2 (TPE)",
  acc_S3 = "Selection accuracy in S3 (TPE)",
  acc_S4 = "Selection accuracy in S4 (TPE)"
)

dir.create("outputs/plots", showWarnings = FALSE, recursive = TRUE)

# ------- Genetic gain -------
# Future years, starting from 0 at the end of the burn-in
p1 <- summary_results |>
  filter(variable %in% c("gain_S1", "gain_S5"), scenario %in% scenarios) |>
  bind_rows(expand.grid(
    scenario = factor(scenarios, levels = c("BurnIn", scenarios)),
    variable = c("gain_S1", "gain_S5"),
    year = last_burnin,
    mean = 0,
    lower = 0,
    upper = 0,
    stringsAsFactors = FALSE
  )) |>
  ggplot(aes(year, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = colors_scenario) +
  scale_fill_manual(values = colors_scenario) +
  facet_wrap(~variable, scales = "free_y", labeller = as_labeller(labels_variable)) +
  labs(x = "Year", y = NULL, colour = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  "outputs/plots/genetic_gain.png",
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
    scenario %in% scenarios
  ) |>
  ggplot(aes(year, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~variable, labeller = as_labeller(labels_variable)) +
  scale_colour_manual(values = colors_scenario) +
  scale_fill_manual(values = colors_scenario) +
  labs(x = "Year", y = NULL, colour = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  "outputs/plots/selection_accuracy.png",
  p2,
  width = 15,
  height = 5,
  dpi = 300,
  bg = "white"
)
