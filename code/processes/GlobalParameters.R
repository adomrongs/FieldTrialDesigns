# Ge Parameters
p <- 1000 # Number of Env in the TPE
k <- 5 # Number of Multiplicative Terms
rho <- 0.4 # Base correlation
gamma <- -0.3 # Skewness of the distribution

#
nReps <- 30
nCycles <- 20 # burn-in years (PS)
nFuture <- 30 # additional years of the design/GS scenarios (loaded from the burn-in RData)
nStages <- 4 # years to fill the pipeline (S1 and S2 are evaluated the same year)
nParents <- 50
nQtl <- 100 # QTL per chromosome
nSnp <- 500 # SNPs per chromosome on the chip (no overlap with the QTL)
maxnEnv <- 5
nYearsTrain <- 5 # Years of plot records kept for GS training (sliding window)
nCrosses <- 100
nProgeny <- 200

# Stages: S1 seedlings -> S2 visual-selected seedlings (same year as S1) -> S3 -> S4 -> S5
nRepS1 <- 1
nRepS2 <- 1
nRepS3 <- 3
nRepS4 <- 3
nRepS5 <- 4

nEnvS1 <- 1
nEnvS2 <- 1
nEnvS3 <- 1
nEnvS4 <- 1
nEnvS5 <- 4

# Plot error variances. S1 and S2 calibrated from the realised h2 at the end of the
# burn-in (h2 = sg2 / (sg2 + varE) -> varE = sg2 * (1 - h2) / h2):
#   S1: h2 = 0.17 with varE = 10 -> sg2 ~ 2.0 -> target h2 ~ 0.1 -> varE ~ 18
#   S2: h2 = 0.16 with varE = 9  -> sg2 ~ 1.7 -> target h2 ~ 0.2 -> varE ~ 6.5
varES1 <- 18 # target h2 S1 = 0.1 (single plot)
varES2 <- 6.5 # target h2 S2 = 0.2 (single plot)
varES3 <- 8 # TODO: calibrate to target h2 S3 = 0.6 (mean of 3 reps)
varES4 <- 4 # TODO: calibrate to target h2 S4 = 0.6 (mean of 3 reps)
varES5 <- 2 # TODO: calibrate to target h2 S5 = 0.8 (mean of 4 reps x 4 locs)

nColS1 <- 200
nColS2 <- 50
nColS3 <- 60
nColS4 <- 12
nColS5 <- 4

nRowS1 <- 100
nRowS2 <- 20
nRowS3 <- 10
nRowS4 <- 5
nRowS5 <- 5

nS1 <- 20000
nS2 <- 1000
nS3 <- 200
nS4 <- 20
nS5 <- 5

# S3 p-rep design (odw): 25 % of the nS3 entries with 2 plots, the rest with 1
pRepS3 <- 0.25
nRepS3pREP <- pRepS3 * nS3 # 50 replicated entries -> 50 x 2 + 150 x 1 = 250 plots
nRowS3pREP <- 10 # two resolvable halves of 5 x 25 (rows 1-5 and 6-10), each
nColS3pREP <- 25 # replicated entry once per half

# S3 augmented design: nChecksS3 checks (best S4 at the end of the burn-in, fixed)
# in each of nBlocksAugS3 blocks -> nS3 + nChecksS3 * nBlocksAugS3 = 240 plots
nChecksS3 <- 4
nBlocksAugS3 <- 10
nRowS3Aug <- 12 # 10 blocks of 6 x 4 (2 bands of rows x 5 bands of columns)
nColS3Aug <- 20

# odw design model of S3 (p-rep and augmented, a new design every year):
# genotypes random with the genomic relationship matrix, variance ratio
# gamma = sg2 / varE per plot, from the target h2 of S3 on the mean of nRepS3 plots:
# h2 = gamma / (gamma + 1 / r) -> gamma = h2 / (r * (1 - h2)) = 0.6 / (3 * 0.4) = 0.5
h2S3 <- 0.6
gammaS3 <- h2S3 / (nRepS3 * (1 - h2S3))

# AR1 correlations of the plot errors of S3: the same values in the odw design
# model and in the FieldSimR errors of the designed trials (PhenotypeDesign.R)
corColS3 <- 0.5
corRowS3 <- 0.7
