# Script to run the phenotypic selection burn-in (one replicate) and save its
# workspace, the starting point of the future scenarios

# ------- Prepare Data -------
# Load packages
library(AlphaSimR)
library(FieldSimR)
library(dplyr)
source("code/processes/GlobalParameters.R")
source("code/processes/HelperFunctions.R")

# Define variables
args <- commandArgs(trailingOnly = TRUE)
Rep <- if (length(args) > 0) as.integer(args[1]) else 1
cat("\nWorking on replicate:", Rep, "\n")
Scenario <- "BurnIn"
S3design <- "RCBD"

# Build TPE (same seed in every replicate), then a seed per replicate
set.seed(2026)
source("code/processes/BuildGe.R")
set.seed(Rep)

# Simulate founders
source("code/processes/CreateFounders.R")

# Define output
source("code/processes/StoreVariables.R")
output$rep <- Rep

# Sample MET environments from the TPE (one sample per simulation year)
covs_tpe <- input_asr$cov.mat
samples_met <- sample_met(
  nenvs = p,
  nsamples = nStages + nCycles + nFuture,
  sample.size = maxnEnv,
  cov.mat = covs_tpe
)$sample

# Training data for GS (plot records of S3-S5 + genotyped pops)
train_data <- list(pheno = NULL, pop = NULL)

# ------- Fill Pipeline -------
source("code/processes/FillPipeline.R")

# ------- Run Burn-in -------
for (year in 1:nCycles) {
  cat("Working on year: ", year)
  simYear <- nStages + year

  source("code/processes/UpdateParents.R")

  source("code/processes/AdvanceYear.R")
}

# ------- Save Results -------
# Workspace of this replicate (loaded by the future scenarios)
dir.create(
  file.path("outputs", "burnin"),
  showWarnings = FALSE,
  recursive = TRUE
)
save.image(
  file.path("outputs", "burnin", sprintf("%s_rep%02d.RData", Scenario, Rep))
)

# Output of this replicate 
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
