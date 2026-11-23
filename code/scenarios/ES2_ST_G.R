# Script to run the scenario ES2_ST_G: the 200 entries of ES2 in S3 sparse testing in 4
# sites (20 commons + 45 uniques per site, 260 plots); S4 p-rep in 1 site (25 plots).
# S3 + S4 = 65,000 + 37,500 = 102,500 EUR/year; the other 137,500 EUR of the S3 + S4
# budget (240,000) go to genotyping S2: 1250 + 3437 = 4687 S2 (187,480 EUR).
# Total S2 + S3 + S4 = 289,980 EUR/year (290,000 in the other scenarios)

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
Scenario <- "ES2_ST_G"

# Designs of S3 and S4
S3design <- "sparse" # 20 commons + 45 uniques per site (5 x 13)
nSitesS3 <- 4
nS3 <- 200
nCommonS3 <- 20
nRowS3 <- 5
nColS3 <- 13
S4design <- "pREP" # 1 site x (20 + 5 replicated)
nSitesS4 <- 1

# S2 genotyped: the 1250 of the other GS scenarios + 137,500 EUR / 40 EUR = 3437
# (also used by AdvanceYearGS.R). Select them again from the S1 of the last burn-in year
nS2GS <- 4687
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
