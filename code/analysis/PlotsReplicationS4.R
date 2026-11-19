# Script to plot the experiment on the replication of the S4 p-rep (I3): replicated
# entries rotating between sites (I3), the same random ones (I3_FixedRandom) or the
# same ones chosen with the GRM (I3_FixedGRM). Figures are saved in outputs/plots

# ------- Prepare Data -------
# Load packages
rm(list = ls())
library(tidyverse)

# Define scenarios and colours
scenarios <- c("I3", "I3_FixedRandom", "I3_FixedGRM")
colors_scenario <- c(
  I3 = "#c0392b", # rotating: red
  I3_FixedRandom = "#2e86c1", # fixed random: blue
  I3_FixedGRM = "#27ae60" # fixed GRM: green
)

# Read outputs
burnin <- read.csv("outputs/BurnIn_output.csv")
future <- bind_rows(
  lapply(file.path("outputs", paste0(scenarios, "_output.csv")), read.csv)
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

# Mean and 95 % CI (t distribution) of every variable across the replicates,
# per scenario and year
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

# ------- Selection accuracy -------
p1 <- summary_results |>
  filter(variable %in% c("acc_S2", "acc_S3", "acc_S4")) |>
  ggplot(aes(year, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~variable, labeller = as_labeller(labels_variable)) +
  scale_colour_manual(values = colors_scenario) +
  scale_fill_manual(values = colors_scenario) +
  labs(x = "Year", y = NULL, colour = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  "outputs/plots/replicationS4_accuracy.png",
  p1,
  width = 12,
  height = 4.5,
  dpi = 300,
  bg = "white"
)

# ------- Genetic gain -------
# Future years, starting from 0 at the end of the burn-in
p2 <- summary_results |>
  filter(variable %in% c("gain_S1", "gain_S5")) |>
  bind_rows(expand.grid(
    scenario = factor(scenarios, levels = scenarios),
    variable = c("gain_S1", "gain_S5"),
    year = last_burnin,
    mean = 0,
    lower = 0,
    upper = 0,
    stringsAsFactors = FALSE
  )) |>
  ggplot(aes(year, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~variable, scales = "free_y", labeller = as_labeller(labels_variable)) +
  scale_colour_manual(values = colors_scenario) +
  scale_fill_manual(values = colors_scenario) +
  labs(x = "Year", y = NULL, colour = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  "outputs/plots/replicationS4_gain.png",
  p2,
  width = 9,
  height = 4.5,
  dpi = 300,
  bg = "white"
)
