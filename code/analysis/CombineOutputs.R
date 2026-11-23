# =============================================================================
# CombineOutputs.R
# Joins the per-replicate results of each scenario (outputs/<Scenario>/
# <Scenario>_repXX.csv, one file per cluster job) into outputs/<Scenario>_output.csv
# =============================================================================
library(dplyr)

scenarios <- c(
  "BurnIn",
  "PS_RCBD",
  "GS_RCBD",
  "ES2",
  "ES2_ST",
  "I1",
  "I2",
  "I3",
  "I4",
  "ES1",
  "ES2_ST_G"
)

for (scenario in scenarios) {
  files <- list.files(
    file.path("outputs", scenario),
    pattern = "_rep[0-9]+\\.csv$",
    full.names = TRUE
  )
  if (length(files) == 0) next
  combined <- bind_rows(lapply(files, read.csv))
  write.csv(
    combined,
    file.path("outputs", paste0(scenario, "_output.csv")),
    row.names = FALSE
  )
  cat(scenario, ":", length(files), "replicates\n")
}
