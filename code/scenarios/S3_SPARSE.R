# Script to run the scenario S3_SPARSE: S3 sparse testing in 6 sites (540 entries, 54
# commons + 81 uniques per site, 810 plots); S4 p-rep in 1 site (25 plots).
# S3 + S4 = 240,000 EUR/year (16 % to S4); S2 = 1250 genotyped (50,000 EUR)

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
Scenario <- "S3_SPARSE"

# Designs of S3 and S4
S3design <- "sparse" # 54 commons + 81 uniques per site (9 x 15)
nSitesS3 <- 6
nS3 <- 540
nCommonS3 <- 54
nRowS3 <- 9
nColS3 <- 15
S4design <- "pREP" # 1 site x (20 + 5 replicated)
nSitesS4 <- 1

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
