# Script to run the scenario GENO_SPARSE: S3_SPARSE where the extra entries are new candidates
# (same 16 % from S2 to S3 as S3_MAX): S2 2500 genotyped (100,000 EUR); S3 sparse testing
# in 6 sites (400 entries, 40 commons + 60 uniques per site, 600 plots, 150,000 EUR); S4 p-rep
# in 1 site (25 plots, 37,500 EUR). Total S2 + S3 + S4 = 287,500 EUR/year (= S3_MAX)

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
Scenario <- "GENO_SPARSE"

# Designs of S3 and S4
S3design <- "sparse" # 40 commons + 60 uniques per site (10 x 10)
nSitesS3 <- 6
nS3 <- 400
nCommonS3 <- 40
nRowS3 <- 10
nColS3 <- 10
S4design <- "pREP" # 1 site x (20 + 5 replicated)
nSitesS4 <- 1

# S2 genotyped: 2500 instead of 1250 (also used by AdvanceYearGS.R). Select them again
# from the S1 of the last burn-in year
nS2GS <- 2500
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
