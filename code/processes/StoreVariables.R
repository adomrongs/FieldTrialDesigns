# Store Variables
# Only S1 is tracked (population used to compare scenarios).
# Filled each year with save_output() (see HelperFunctions.R)
output <- data.frame(
  year = 1:nCycles,
  rep = numeric(nCycles),
  scenario = Scenario,
  mean_gv_tpe = numeric(nCycles),
  mean_gv_met = numeric(nCycles),
  var_gv_tpe = numeric(nCycles),
  var_gv_met = numeric(nCycles),
  acc_tpe = numeric(nCycles),
  acc_met = numeric(nCycles),
  met_tpe = numeric(nCycles),
  h2 = numeric(nCycles),
  acc_S3_met = numeric(nCycles),
  acc_S3_tpe = numeric(nCycles),
  resp_S4 = numeric(nCycles)
)
