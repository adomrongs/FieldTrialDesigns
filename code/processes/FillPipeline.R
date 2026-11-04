# Fill the breeding pipeline: after nStages years each stage holds a different cohort
# (S1 and S2 are evaluated in the same year)
for (year in 1:nStages) {
  cat("FillPipeline year:", year, "of", nStages, "\n")
  simYear <- year # index into samples_met

  if (year <= 4) {
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
  }

  if (year <= 3) {
    S3 <- selectInd(pop = S2, nInd = nS3)
    source("code/processes/PhenotypeS3.R") # RCBD (S3design of the burn-in)
    S3 <- res$pop[res$pop@id %in% S3@id]
  }

  if (year <= 2) {
    S4 <- selectInd(pop = S3, nInd = nS4)
    S4 <- phenotype_met(
      pop = S4,
      envs = met_envs(samples_met, simYear, nEnvS4),
      cov.mat = covs_tpe,
      nreps = nRepS4,
      varR = varES4,
      ncols = nColS4,
      nrows = nRowS4
    )
  }

  if (year <= 1) {
    S5 <- selectInd(pop = S4, nInd = nS5)
    S5 <- phenotype_met(
      pop = S5,
      envs = met_envs(samples_met, simYear, nEnvS5),
      cov.mat = covs_tpe,
      nreps = nRepS5,
      varR = varES5,
      ncols = nColS5,
      nrows = nRowS5
    )
  }
}
