# Script to run the scenario I3: intermediate 3: S3 sparse testing in 6 sites (240 entries, 360 plots); S4 p-rep in 4 sites (100 plots).
# Budget of S3 + S4 = 240,000 EUR/year (62.5 % to S4)

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
Scenario <- "I3"

# Designs of S3 and S4
S3design <- "sparse" # 24 commons + 36 uniques per site (6 x 10)
nSitesS3 <- 6
nS3 <- 240
nCommonS3 <- 24
nRowS3 <- 6
nColS3 <- 10
S4design <- "pREP" # 4 sites x (20 + 5 replicated)
nSitesS4 <- 4

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
