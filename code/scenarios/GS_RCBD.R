# Script to run the scenario GS_RCBD: Actual program with GS: S3 RCBD (1 site x 3 blocks, 600 plots) and S4 RCBD (1 site x 3 blocks, 60 plots).
# Budget of S3 + S4 = 240,000 EUR/year (37.5 % to S4)

# ------- Prepare Data -------
# Load packages
library(AlphaSimR)
library(FieldSimR)
library(sommer)
library(odw)
library(dplyr)
source("code/processes/GlobalParameters.R")

# Define variables
args <- commandArgs(trailingOnly = TRUE)
Rep <- if (length(args) > 0) as.integer(args[1]) else 1
cat("\nWorking on replicate:", Rep, "\n")

load(file.path("outputs", "burnin", sprintf("BurnIn_rep%02d.RData", Rep)))
source("code/processes/GlobalParameters.R")
source("code/processes/HelperFunctions.R")
Scenario <- "GS_RCBD"

# Designs of S3 and S4
S3design <- "RCBD" # 200 entries, 1 site x 3 blocks (60 x 10)
S4design <- "RCBD" # 20 entries, 1 site x 3 blocks (12 x 5)

# GS genotypes nS2GS S2 (same budget as the 500 phenotyped in PS): select them
# again from the S1 of the last burn-in year
S2 <- selectInd(pop = S1, nInd = nS2GS)

# Define output
output_years <- (nCycles + 1):(nCycles + nFuture)
source("code/processes/StoreVariables.R")

# ------- Run Future scenario -------
for (year in output_years) {
  cat("Working on year: ", year)
  simYear <- nStages + year

  source("code/processes/FitModel.R")

  source("code/processes/UpdateParents.R")

  source("code/processes/AdvanceYearGS.R")
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
