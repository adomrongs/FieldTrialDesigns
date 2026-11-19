# Script to run the scenario ES1: extreme 1 (minimum to S3): S3 p-rep in 1 site (200 entries + 10 replicated, 210 plots); S4 p-rep in 5 sites (125 plots).
# Budget of S3 + S4 = 240,000 EUR/year (78 % to S4)

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
Scenario <- "ES1"

# Designs of S3 and S4
S3design <- "pREP" # 1 site, 200 entries + 10 replicated (14 x 15)
nRepS3pREP <- 10
nRowS3 <- 14
nColS3 <- 15
S4design <- "pREP" # 5 sites x (20 + 5 replicated)
nSitesS4 <- 5

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
