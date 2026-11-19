# Advance breeding program by 1 year (phenotypic selection: burn-in and PS_RCBD)
# Order S5 -> S1: each stage selects from last year's population before it is replaced
# (S2 is selected from S1 and evaluated in the same year, so it goes after S1)
cat(" Advance year \n")
index <- output$year == year

# S5
output$acc_S4[index] <- safe_cor(S4@pheno[, 1], tpe_value(S4, cov_tpe_on))
S5 <- selectInd(pop = S4, nInd = nS5) # mean of the plots of each entry
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

# S4: phenotyped with the design of the scenario (S4design, PhenotypeS4.R)
output$acc_S3[index] <- safe_cor(S3@pheno[, 1], tpe_value(S3, cov_tpe_on))
S4 <- selectInd(pop = S3, nInd = nS4) # mean of the plots of each entry
source("code/processes/PhenotypeS4.R")
train_data <- store_training(train_data, res, year, "S4")
S4 <- res$pop[res$pop@id %in% S4@id]

# S3: phenotyped with the design of the scenario (S3design, PhenotypeS3.R)
output$acc_S2[index] <- safe_cor(S2@pheno[, 1], tpe_value(S2, cov_tpe_on))
S3 <- selectInd(pop = S2, nInd = nS3)
source("code/processes/PhenotypeS3.R")
train_data <- store_training(train_data, res, year, "S3")
S3 <- res$pop[res$pop@id %in% S3@id]

# S1 (not genotyped: its records are not stored for GS)
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

# S2: seedlings that pass the visual selection, evaluated the same year
# (not stored for GS)
S2 <- selectInd(pop = S1, nInd = nS2)
S2 <- phenotype_met(
  pop = S2,
  envs = met_envs(samples_met, simYear, nEnvS2),
  cov.mat = covs_tpe,
  nreps = nRepS2,
  varR = varES2,
  ncols = nColS2,
  nrows = nRowS2
)

# Sliding window: keep only the last nYearsTrain years of training data
first_year <- year - nYearsTrain # records of this year and older are dropped
train_data$pheno <- train_data$pheno |>
  filter(year > first_year)
train_data$pop <- train_data$pop[which(
  train_data$pop@id %in% train_data$pheno$id
)]
