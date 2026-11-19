# Base correlation matrix: TPE environments first, then the off-TPE ones
# (e.g. El Nino years), as in Werner et al. (2025)
nOff <- round(pOff * p)
nOn <- p - nOff
envs_on <- seq_len(nOn) # environments of the TPE (the target of selection)
base_mat <- matrix(rhoBetween, p, p)
base_mat[envs_on, envs_on] <- rhoOn
base_mat[-envs_on, -envs_on] <- rhoOff
diag(base_mat) <- 1

# Create diagonal variance matrix
De <- skew_diag_mat(
  n = p,
  shape = 1.5,
  scale = 1
)

# Create correlation matrix (structured noise around the base correlations)
Ce <- struc_cor_mat(
  base.mat = base_mat,
  range = rangeCor,
  rank = k,
  pos.def = TRUE
)
Ce <- (Ce + t(Ce)) / 2

# Generate Ge
Ge <- sqrt(De) %*% Ce %*% sqrt(De)
# cov2cor(Ge)[upper.tri(cov2cor(Ge), diag = F)] |>  hist(breaks = 200)
