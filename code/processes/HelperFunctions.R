# =============================================================================
# HelperFunctions.R
#
# All the project's functions live in this script.
# Conventions:
#   - Each scenario script (code/scenarios/*.R) sources this file at the start,
#     right after GlobalParameters.R. The process scripts (code/processes/*.R)
#     run in that same session, so they don't need to source it again.
#   - Each function is introduced by a header "# -------- name --------"
#     followed by at most two comment lines explaining what it does.
#   - Don't define functions in other scripts: if you need a new one, add it here.
#   - Plain data.frames + dplyr (group_by, summarise, left_join, bind_rows...);
#     no data.table, with() or tapply(). Pipe: always base R |> (not %>%).
# =============================================================================

# --------------------------------- met_envs ----------------------------------
# Returns the first 'nenvs' environments (rows of cov.mat) sampled for the
# simulation year 'sim_year' (one sample per year: fill pipeline + cycles).
met_envs <- function(samples_met, sim_year, nenvs) {
  samples_met[[sim_year]][seq_len(nenvs)]
}


# --------------------------------- safe_cor ----------------------------------
# Correlation that returns NA (with no warning) when there are fewer than
# 3 individuals, missing values or one of the vectors has zero variance.
safe_cor <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) < 3 || anyNA(x) || anyNA(y)) {
    return(NA_real_)
  }
  if (sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }
  cor(x, y)
}


# --------------------------------- tpe_value ---------------------------------
# Genetic value of each individual in the TPE: main effect + slopes x mean
# covariates of the TPE environments (cov_mean = cov_tpe_on, without off-TPE).
tpe_value <- function(pop, cov_mean) {
  as.numeric(pop@gv[, 1] + pop@gv[, -1, drop = FALSE] %*% cov_mean)
}


# -------------------------------- vanraden_grm --------------------------------
# Tetraploid VanRaden GRM from a 0-4 dosage matrix (rownames = ids), markers
# coded by the minor allele; scale_diag = TRUE scales it to mean diagonal 1.
vanraden_grm <- function(M, scale_diag = FALSE) {
  M <- M[, colMeans(M^2) - colMeans(M)^2 > 1e-8, drop = FALSE]
  freq <- colMeans(M) / 4
  flip <- freq > 0.5
  M[, flip] <- 4 - M[, flip]
  freq[flip] <- 1 - freq[flip]
  W <- sweep(M, 2, 4 * freq)
  G <- tcrossprod(W) / sum(4 * freq * (1 - freq))
  if (scale_diag) {
    G <- G / mean(diag(G))
  }
  G + diag(1e-3, nrow(G)) # make it invertible
}


# ------------------------------- phenotype_met -------------------------------
# Simulates a MET trial (GxE + AR1 plot errors, randomised by FieldSimR) in the
# 'envs' environments and writes the genotype means to pop@pheno.
phenotype_met <- function(
  pop,
  envs, # rows of cov.mat to evaluate, e.g. met_envs(samples_met, simYear, nEnvS5)
  cov.mat, # covs_tpe (all rows of the TPE)
  nreps,
  varR,
  ncols,
  nrows,
  spatial.model = "AR1",
  return.pheno = FALSE
) {
  # Check the environments before simulating anything
  if (length(envs) == 0 || anyNA(envs) || any(envs > nrow(cov.mat))) {
    stop(
      "phenotype_met: invalid 'envs' (",
      paste(envs, collapse = ", "),
      "). Check simYear and length(samples_met)."
    )
  }

  # Genetic values only in the sampled environments (cov.mat subset).
  # multi_asr_output requires nenvs > 1: with a single environment, add an
  # auxiliary one and drop it afterwards.
  envs_calc <- if (length(envs) == 1) {
    c(envs, setdiff(seq_len(nrow(cov.mat)), envs)[1])
  } else {
    envs
  }
  gv_met <- multi_asr_output(
    pop = pop,
    ntraits = 1,
    nenvs = length(envs_calc),
    nreps = nreps,
    cov.mat = cov.mat[envs_calc, , drop = FALSE]
  )
  gv_met <- gv_met |>
    filter(as.integer(as.character(env)) <= length(envs)) |>
    droplevels()

  # Plot errors
  error_met <- field_trial_error(
    nenvs = length(envs),
    nblocks = nreps,
    varR = varR,
    spatial.model = spatial.model,
    ncols = ncols,
    nrows = nrows
  )

  pheno_met <- make_phenotypes(
    gv.df = gv_met,
    error.df = error_met,
    randomise = TRUE
  )

  # Genotype means, in the order of pop@id
  pheno_means <- pheno_met |>
    mutate(id = sprintf("%.0f", as.numeric(as.character(id)))) |> # no "1e+05"
    group_by(id) |>
    summarise(mean_pheno = mean(y.Trait1), .groups = "drop")
  n_na <- sum(!(pop@id %in% pheno_means$id))
  if (n_na > 0) {
    stop(
      "phenotype_met: ", n_na, " of ", pop@nInd, " individuals have no phenotype.",
      " Example pop@id: ", paste(head(pop@id, 3), collapse = ", "),
      " | example id in the trial: ", paste(head(pheno_means$id, 3), collapse = ", ")
    )
  }
  pop@pheno[, 1] <- pheno_means$mean_pheno[match(pop@id, pheno_means$id)]

  if (return.pheno) {
    return(list(pop = pop, pheno = pheno_met))
  }
  pop
}


# ------------------------------ store_training -------------------------------
# Adds a trial's plot records to train_data$pheno (loc = "L1", "L2"...; env = one
# trial = "stage_year_loc") and its new individuals to train_data$pop (no repeats).
store_training <- function(train_data, res, year, stage) {
  pheno <- res$pheno |>
    mutate(
      loc = paste0("L", env), # site k of a year is the k-th environment sampled that year
      env = paste(stage, year, loc, sep = "_"),
      id = sprintf("%.0f", as.numeric(as.character(id))) # same format as pop@id (no "1e+05")
    )
  pheno$year <- year
  pheno$stage <- stage
  train_data$pheno <- bind_rows(train_data$pheno, pheno)

  if (is.null(train_data$pop)) {
    train_data$pop <- res$pop
  } else {
    new_ind <- which(!(res$pop@id %in% train_data$pop@id))
    if (length(new_ind) > 0) {
      train_data$pop <- mergePops(list(train_data$pop, res$pop[new_ind]))
    }
  }
  train_data
}
