# Ge Parameters
p <- 1000 # Number of Env in the TPE
k <- 5 # Number of Multiplicative Terms
pOff <- 0.1 # Proportion of off-TPE environments (Werner et al. 2025)
rhoOn <- 0.6 # Base genetic correlation within the TPE
rhoOff <- 0.4 # Base genetic correlation within the off-TPE environments
rhoBetween <- -0.5 # Base genetic correlation between TPE and off-TPE
rangeCor <- 0.3 # Range of the correlations around the base values

#
nReps <- 30
nCycles <- 20 # burn-in years (PS)
nFuture <- 30 # additional years of the scenarios (loaded from the burn-in RData)
nStages <- 4 # years to fill the pipeline (S1 and S2 are evaluated the same year)
nParents <- 50
nQtl <- 100 # QTL per chromosome
nSnp <- 500 # SNPs per chromosome on the chip (no overlap with the QTL)
maxnEnv <- 6 # environments sampled per year (6 locations: S5 is tested in all)
nYearsTrain <- 5 # Years of plot records kept for GS training (sliding window)
nCrosses <- 100
nProgeny <- 200
residualBy <- "trial" # one residual variance per "trial" or per "stage" (FitModel.R)

# Costs (EUR). Annual budget of S3 + S4 = 600 x 250 + 60 x 1500 = 240,000 in
# every scenario; S2 = 50,000 in every scenario (500 phenotyped or 1250 genotyped)
costPlotS2 <- 100
costPlotS3 <- 250
costPlotS4 <- 1500
costGeno <- 40

# Stages: S1 seedlings -> S2 visual-selected seedlings (same year as S1) -> S3 -> S4 -> S5
nS1 <- 20000
nS2 <- 500 # PS: S2 phenotyped (500 x 100 EUR)
nS2GS <- 1250 # GS: S2 genotyped (1250 x 40 EUR)
nS3 <- 200 # entries in S3 (changed by the sparse testing scenarios)
nS4 <- 20
nS5 <- 5

nRepS1 <- 1
nRepS2 <- 1
nRepS3 <- 3 # replicates per site of the S3 RCBD
nRepS4 <- 3 # replicates per site of the S4 RCBD
nRepS5 <- 4

nEnvS1 <- 1
nEnvS2 <- 1
nSitesS3 <- 1 # sites of S3 (changed by the scenarios)
nSitesS4 <- 1 # sites of S4 (changed by the scenarios)
nEnvS5 <- 6

# Plot error variances, calibrated to a target reliability r of the entry means
# with Reliability.R (2 rounds, 05/11/2026): new varE = old varE * (error of the
# means needed for r) / (error of the means observed), error = var_g * (1 - r) / r
varES1 <- 18 # target h2 S1 = 0.1 (single plot)
varES2 <- 6.5 # reliability S2 ~ 0.15 (single plot; observed 0.14)
varES3 <- 1.7 # reliability S3 ~ 0.6 (mean of 3 reps; 0.55 with varE = 2.1)
varES4 <- 0.5 # reliability S4 ~ 0.8 (mean of 3 reps; 0.74 with varE = 0.7)
varES5 <- 0.5 # reliability S5 ~ 0.9 (mean of 4 reps x 4 locs; 0.84 with varE = 0.8)

# Field of one site (cols x rows); S3 and S4 are changed by the scenarios
nColS1 <- 200
nColS2 <- 25
nColS3 <- 60
nColS4 <- 12
nColS5 <- 4

nRowS1 <- 100
nRowS2 <- 20
nRowS3 <- 10
nRowS4 <- 5
nRowS5 <- 5

# S4 p-rep in each site: the 20 entries + 5 replicated (25 %) in 5 x 5 plots.
# Replicated entries: "rotating" (5 different per site), "fixed_random" or
# "fixed_GRM" (the same 5 in every site)
nRepS4pREP <- 5
nRowS4pREP <- 5
nColS4pREP <- 5
S4repAlloc <- "fixed_GRM" # rotating / fixed_random / fixed_GRM gave very similar results (19/11/2026)

# odw design model: genotypes random with the GRM, variance ratio
# gamma = sg2 / varE per plot, from the target reliability of the mean of nRep plots:
# h2 = gamma / (gamma + 1 / r) -> gamma = h2 / (r * (1 - h2))
h2S3 <- 0.6
gammaS3 <- h2S3 / (nRepS3 * (1 - h2S3)) # 0.6 / (3 * 0.4) = 0.5
h2S4 <- 0.8
gammaS4 <- h2S4 / (nRepS4 * (1 - h2S4)) # 0.8 / (3 * 0.2) = 1.33

# AR1 correlations of the plot errors of the designed trials (S3 and S4): the
# same values in the odw design model and in the FieldSimR errors (PhenotypeDesign.R)
corColS3 <- 0.5
corRowS3 <- 0.7
