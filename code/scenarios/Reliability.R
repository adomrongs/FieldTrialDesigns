# TEMPORARY script to track the reliability of the entry means of S2-S5 (to tune
# the error variances). Same as PS_RCBD (the only future scenario that
# phenotypes S2); delete it once varES2-varES5 are tuned

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
Scenario <- "Reliability"
S3design <- "RCBD"

# Define output: one row per year and stage
output <- data.frame(
  year = (nCycles + 1):(nCycles + nFuture),
  rep = Rep,
  scenario = Scenario
)
reliability <- NULL

# ------- Run Future scenario -------
for (year in (nCycles + 1):(nCycles + nFuture)) {
  cat("Working on year: ", year)
  simYear <- nStages + year

  source("code/processes/UpdateParents.R")

  source("code/processes/AdvanceYear.R")

  # Reliability of the entry means of each stage phenotyped this year:
  # genetic value of each entry = mean over the environments of its trial
  for (stage in c("S2", "S3", "S4", "S5")) {
    pop <- get(stage)
    envs <- met_envs(samples_met, simYear, get(paste0("nEnv", stage)))
    gv_met <- rowMeans(
      pop@gv[, 1] + pop@gv[, -1, drop = FALSE] %*% t(covs_tpe[envs, , drop = FALSE])
    )
    reliability <- bind_rows(
      reliability,
      data.frame(
        rep = Rep,
        year = year,
        stage = stage,
        n_ind = pop@nInd,
        n_reps = get(paste0("nRep", stage)),
        n_envs = length(envs),
        var_e = get(paste0("varE", stage)),
        var_g = var(gv_met),
        var_pheno = var(pop@pheno[, 1]),
        reliability = cor(pop@pheno[, 1], gv_met)^2
      )
    )
  }
}

# Save the results of this replicate
dir.create(
  file.path("outputs", Scenario),
  showWarnings = FALSE,
  recursive = TRUE
)
write.csv(
  reliability,
  file.path("outputs", Scenario, sprintf("%s_rep%02d.csv", Scenario, Rep)),
  row.names = FALSE
)
