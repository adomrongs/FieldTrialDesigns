# Script to run the PS scenario (reference): S3 RCBD (1 site x 3 blocks, 600
# plots) and S4 RCBD (1 site x 3 blocks, 60 plots), selection on phenotype means

# ------- Prepare Data -------
# Load packages
library(AlphaSimR)
library(FieldSimR)
library(dplyr)
source("code/processes/GlobalParameters.R")

# Define variables
args <- commandArgs(trailingOnly = TRUE)
Rep <- if (length(args) > 0) as.integer(args[1]) else 1
cat("\nWorking on replicate:", Rep, "\n")

load(file.path("outputs", "burnin", sprintf("BurnIn_rep%02d.RData", Rep)))
source("code/processes/GlobalParameters.R")
source("code/processes/HelperFunctions.R")
Scenario <- "PS_RCBD"
S3design <- "RCBD" # 200 entries, 1 site x 3 blocks (60 x 10)
S4design <- "RCBD" # 20 entries, 1 site x 3 blocks (12 x 5)

# Define output
output_years <- (nCycles + 1):(nCycles + nFuture)
source("code/processes/StoreVariables.R")

# ------- Run Future scenario -------
for (year in output_years) {
  cat("Working on year: ", year)
  simYear <- nStages + year

  source("code/processes/UpdateParents.R")

  source("code/processes/AdvanceYear.R")
}

# Save the results of this replicate (combined with code/analysis/CombineOutputs.R)
dir.create(
  file.path("outputs", Scenario),
  showWarnings = FALSE,
  recursive = TRUE
)
write.csv(
  output,
  file.path("outputs", Scenario, sprintf("%s_rep%02d.csv", Scenario, Rep)),
  row.names = FALSE
)
