# =============================================================================
# CostTable.R
# Annual cost of each simulated scenario (EUR, 1 USD = 1 EUR), from the plot and
# genotyping costs in GlobalParameters.R. Writes outputs/cost_table.csv, used by the
# plotting scripts to express genetic gain per unit cost.
# Assumptions:
#   - The number of individuals per stage is the same in every scenario; only the
#     S3 design (plots), S2 phenotyping (PS) and S2 genotyping (GS) change.
#   - GS: the nS2 S2 candidates are genotyped every year; S3-S5 come from them, so
#     they are not genotyped again. The 4 augmented checks are genotyped once (ignored).
#   - Not included yet: cloning per plant, cost per location, S1 seedlings (TODO).
# =============================================================================
library(dplyr)

source("code/processes/GlobalParameters.R")

# Plots of each stage per year (S4 and S5: individuals x reps x locations)
plots_S1 <- nS1 * nRepS1 * nEnvS1
plots_S4 <- nS4 * nRepS4 * nEnvS4
plots_S5 <- nS5 * nRepS5 * nEnvS5

cost_table <- data.frame(
  scenario = c("PS_RCBD", "GS_RCBD", "GS_Augmented", "GS_pREP_R", "GS_pREP_O"),
  selection = c("PS", "GS", "GS", "GS", "GS"),
  plots_S1 = plots_S1,
  plots_S2 = c(nS2 * nRepS2 * nEnvS2, 0, 0, 0, 0), # S2 only phenotyped in PS
  plots_S3 = c(
    nS3 * nRepS3, # RCBD
    nS3 * nRepS3, # RCBD
    nS3 + nChecksS3 * nBlocksAugS3, # augmented: entries once + checks in each block
    nS3 + nRepS3pREP, # p-rep: replicated entries twice
    nS3 + nRepS3pREP
  ),
  plots_S4 = plots_S4,
  plots_S5 = plots_S5,
  n_geno = c(0, nS2, nS2, nS2, nS2) # S2 candidates genotyped per year
) |>
  mutate(
    cost_S1 = plots_S1 * costPlotS1,
    cost_S2 = plots_S2 * costPlotS2,
    cost_S3 = plots_S3 * costPlotS3,
    cost_S4 = plots_S4 * costPlotS4,
    cost_S5 = plots_S5 * costPlotS5,
    cost_geno = n_geno * costGeno,
    cost_total = cost_S1 + cost_S2 + cost_S3 + cost_S4 + cost_S5 + cost_geno,
    saving_vs_PS_RCBD = cost_total[scenario == "PS_RCBD"] - cost_total
  )

cost_table

dir.create("outputs", showWarnings = FALSE)
write.csv(cost_table, file.path("outputs", "cost_table.csv"), row.names = FALSE)
