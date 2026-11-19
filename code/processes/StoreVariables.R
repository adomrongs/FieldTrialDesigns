# Store Variables
# One row per year in 'output_years' (burn-in or future years), filled in
# AdvanceYear.R / AdvanceYearGS.R. Values in the TPE (environments of the TPE only):
#   mean_gv_S1, mean_gv_S5: mean genetic value of the S1 / S5 of the year (gain)
#   acc_S2, acc_S3, acc_S4: cor(selection criterion, genetic value) among the
#     candidates of that stage (GEBV in GS; phenotype mean in PS)
output <- data.frame(
  year = output_years,
  rep = Rep,
  scenario = Scenario,
  mean_gv_S1 = NA_real_,
  mean_gv_S5 = NA_real_,
  acc_S2 = NA_real_,
  acc_S3 = NA_real_,
  acc_S4 = NA_real_
)
