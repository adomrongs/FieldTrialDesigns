# Advance breeding program by 1 year in the GS scenarios: S3, S4 and S5 are
# selected on the GEBVs of FitModel.R (run at the start of the year)

cat(" Advance year GS \n")
index <- output$year == year

# ------- Select S5 -------
output$acc_S4[index] <- safe_cor(S4@ebv[, 1], tpe_value(S4, cov_tpe_on))
S5 <- selectInd(pop = S4, nInd = nS5, use = "ebv")
res <- phenotype_met(
  pop = S5,
  envs = met_envs(samples_met, simYear, nEnvS5),
  cov.mat = covs_tpe,
  nreps = nRepS5,
  varR = varES5,
  ncols = nColS5,
  nrows = nRowS5,
  return.pheno = TRUE
)
S5 <- res$pop
train_data <- store_training(train_data, res, year, "S5")
output$mean_gv_S5[index] <- mean(tpe_value(S5, cov_tpe_on))

# ------- Select S4 -------
output$acc_S3[index] <- safe_cor(S3@ebv[, 1], tpe_value(S3, cov_tpe_on))
S4 <- selectInd(pop = S3, nInd = nS4, use = "ebv")
source("code/processes/PhenotypeS4.R")
train_data <- store_training(train_data, res, year, "S4")
S4 <- res$pop[res$pop@id %in% S4@id]

# ------- Select S3 -------
output$acc_S2[index] <- safe_cor(S2@ebv[, 1], tpe_value(S2, cov_tpe_on))
S3 <- selectInd(pop = S2, nInd = nS3, use = "ebv")
source("code/processes/PhenotypeS3.R")
train_data <- store_training(train_data, res, year, "S3")
S3 <- res$pop[res$pop@id %in% S3@id]

# ------- Generate S1s -------
S1 <- randCross(Parents, nCrosses, nProgeny)
S1 <- phenotype_met(
  pop = S1,
  envs = met_envs(samples_met, simYear, nEnvS1),
  cov.mat = covs_tpe,
  nreps = nRepS1,
  varR = varES1,
  ncols = nColS1,
  nrows = nRowS1
)
output$mean_gv_S1[index] <- mean(tpe_value(S1, cov_tpe_on))

# ------- Select S2s (genotyped, not phenotyped) -------
S2 <- selectInd(pop = S1, nInd = nS2GS)

# ------- Store TRP data -------
first_year <- year - nYearsTrain
train_data$pheno <- train_data$pheno |>
  filter(year > first_year)
train_data$pop <- train_data$pop[which(
  train_data$pop@id %in% train_data$pheno$id
)]
