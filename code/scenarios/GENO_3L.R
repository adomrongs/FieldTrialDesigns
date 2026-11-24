# Script to run the scenario GENO_3L: S3_MAX with 3 sites instead of 4, the saved money to
# genotype S2. S3: the 200 entries once in each of 3 sites (600 plots, 150,000 EUR); S4 p-rep
# in 1 site (25 plots, 37,500 EUR); S2: 102,500 EUR / 40 = 2562 genotyped.
# Total S2 + S3 + S4 = 289,980 EUR/year

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
Scenario <- "GENO_3L"

# Designs of S3 and S4
S3design <- "RCBD" # 200 entries once in each of 3 sites (20 x 10 per site)
nSitesS3 <- 3
nRepS3 <- 1
nColS3 <- 20
nRowS3 <- 10
S4design <- "pREP" # 1 site x (20 + 5 replicated)
nSitesS4 <- 1

# S2 genotyped: 2562 instead of 1250 (also used by AdvanceYearGS.R). Select them again
# from the S1 of the last burn-in year
nS2GS <- 2562
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
