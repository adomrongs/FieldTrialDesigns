# Advance breeding program by 1 year
# Order S5 -> S1: each stage selects from last year's population before it is replaced
# (S2 is selected from S1 and evaluated in the same year, so it goes after S1)
cat(" Advance year \n")

# S5
S5 <- selectInd(pop = S4, nInd = nS5)
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

# S4
S4 <- selectInd(pop = S3, nInd = nS4) # mean of the plots of each entry

# S3 metrics: accuracy of the values used to select S4 (vs the genetic value in
# the S3 trial environment of last year and in the TPE) and response in S4
env_S3 <- met_envs(samples_met, simYear - 1, nEnvS3)
gv_S3_met <- as.numeric(S3@gv[, 1] + S3@gv[, -1, drop = FALSE] %*% covs_tpe[env_S3, ])
gv_S3_tpe <- as.numeric(S3@gv[, 1] + S3@gv[, -1, drop = FALSE] %*% colMeans(covs_tpe))
output$acc_S3_met[output$year == year] <- safe_cor(S3@pheno[, 1], gv_S3_met)
output$acc_S3_tpe[output$year == year] <- safe_cor(S3@pheno[, 1], gv_S3_tpe)
output$resp_S4[output$year == year] <- mean(gv_S3_tpe[S3@id %in% S4@id]) - mean(gv_S3_tpe)

res <- phenotype_met(
  pop = S4,
  envs = met_envs(samples_met, simYear, nEnvS4),
  cov.mat = covs_tpe,
  nreps = nRepS4,
  varR = varES4,
  ncols = nColS4,
  nrows = nRowS4,
  return.pheno = TRUE
)
S4 <- res$pop
train_data <- store_training(train_data, res, year, "S4")

# S3: phenotyped with the design of the scenario (S3design, PhenotypeS3.R)
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
output <- save_output(df = output, pop = S1, year = year)

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
