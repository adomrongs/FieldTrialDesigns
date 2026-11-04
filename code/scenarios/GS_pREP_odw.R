# Script to run the GS scenario assuming a p-rep design with 25 % of the entries
# replicated (odw with the GRM, a new design every year)

# ------- Prepare Data -------
# Load packages
library(AlphaSimR)
library(FieldSimR)
library(SpATS)
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
Scenario <- "GS_pREP_odw"
S3design <- "pREP_odw"

# Define output
future_years <- (nCycles + 1):(nCycles + nFuture)
output <- data.frame(
  year = future_years,
  rep = Rep,
  scenario = Scenario,
  mean_gv_tpe = NA_real_,
  mean_gv_met = NA_real_,
  var_gv_tpe = NA_real_,
  var_gv_met = NA_real_,
  acc_tpe = NA_real_,
  acc_met = NA_real_,
  met_tpe = NA_real_,
  h2 = NA_real_,
  acc_gs = NA_real_,
  acc_S3_met = NA_real_,
  acc_S3_tpe = NA_real_,
  resp_S4 = NA_real_
)

# ------- Run Future scenario -------
for (year in future_years) {
  cat("Working on year: ", year)
  simYear <- nStages + year

  source("code/processes/FitModel.R")

  source("code/processes/UpdateParents.R")

  source("code/processes/AdvanceYearGS.R")
}

# Save the results of this replicate 
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
