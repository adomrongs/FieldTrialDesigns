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


# ------------------------------- phenotype_met -------------------------------
# Simulates a MET trial (GxE + AR1 plot errors) in the 'envs' environments,
# writes the genotype means to pop@pheno and the trial metrics to pop@miscPop.
phenotype_met <- function(
  pop,
  envs, # rows of cov.mat to evaluate, e.g. met_envs(samples_met, simYear, nEnvS2)
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

  # Genetic values only in the sampled environments (cov.mat subset):
  # avoids computing nInd x p x nreps values that would then be thrown away.
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

  # Genotype means
  pheno_means <- pheno_met |>
    mutate(id = sprintf("%.0f", as.numeric(as.character(id)))) |> # no "1e+05"
    group_by(id) |>
    summarise(mean_pheno_met = mean(y.Trait1), .groups = "drop")

  gv_means <- gv_met |>
    mutate(id = sprintf("%.0f", as.numeric(as.character(id)))) |> # no "1e+05"
    group_by(id) |>
    summarise(mean_gv_met = mean(gv.Trait1), .groups = "drop")

  # Mean genetic value across the whole TPE (same formula multi_asr_output uses):
  # main effect + slopes x mean covariates of the TPE
  tpe_means <- data.frame(
    id = pop@id,
    mean_gv_tpe = as.numeric(
      pop@gv[, 1] + pop@gv[, -1, drop = FALSE] %*% colMeans(cov.mat)
    )
  )

  # Everything in the order of pop@id
  means <- tpe_means |>
    left_join(gv_means, by = "id") |>
    left_join(pheno_means, by = "id")

  # Every individual must have a genetic value and a phenotype
  n_na_gv <- sum(is.na(means$mean_gv_met))
  n_na_pheno <- sum(is.na(means$mean_pheno_met))
  if (n_na_gv > 0 || n_na_pheno > 0) {
    stop(
      "phenotype_met: out of ",
      pop@nInd,
      " individuals, ",
      n_na_gv,
      " have no gv and ",
      n_na_pheno,
      " have no phenotype.",
      " Example pop@id: ",
      paste(head(pop@id, 3), collapse = ", "),
      " | example id in gv_met: ",
      paste(head(gv_met$id, 3), collapse = ", "),
      " | envs in gv_met: ",
      paste(unique(gv_met$env), collapse = ", ")
    )
  }

  pop@pheno[, 1] <- means$mean_pheno_met

  # Population-level metrics (selectInd wipes them; save them first)
  pop@miscPop$mean_gv_tpe <- mean(means$mean_gv_tpe)
  pop@miscPop$mean_gv_met <- mean(means$mean_gv_met)
  pop@miscPop$var_gv_tpe <- var(means$mean_gv_tpe)
  pop@miscPop$var_gv_met <- var(means$mean_gv_met)
  pop@miscPop$acc_tpe <- safe_cor(means$mean_gv_tpe, means$mean_pheno_met) # accuracy w.r.t. the TPE
  pop@miscPop$acc_met <- safe_cor(means$mean_gv_met, means$mean_pheno_met) # accuracy in the trial
  pop@miscPop$met_tpe <- safe_cor(means$mean_gv_met, means$mean_gv_tpe) # MET-TPE alignment
  pop@miscPop$h2 <- pop@miscPop$acc_met^2 # entry-mean heritability

  if (return.pheno) {
    return(list(pop = pop, pheno = pheno_met))
  }
  pop
}


# -------------------------------- save_output --------------------------------
# Copies the metrics that phenotype_met left in pop@miscPop into the row
# for 'year' of the 'output' table. Call it right after phenotype_met (S1 only).
save_output <- function(df, pop, year) {
  vars <- c(
    "mean_gv_tpe",
    "mean_gv_met",
    "var_gv_tpe",
    "var_gv_met",
    "acc_tpe",
    "acc_met",
    "met_tpe",
    "h2"
  )
  index <- which(df$year == year)
  df[index, vars] <- unlist(pop@miscPop[vars])
  df
}


# ------------------------------ store_training -------------------------------
# Adds a trial's plot records to train_data$pheno (loc = "L1", "L2"...; env = one
# trial = "stage_year_loc") and its new individuals to train_data$pop (no repeats).
store_training <- function(train_data, res, year, stage) {
  pheno <- res$pheno |>
    mutate(
      loc = paste0("L", env), # same locations every year (the GxE sample changes)
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


# -------------------------------- avg_grm_self --------------------------------
# TrainSel criterion (Fernandez-Gonzalez et al. 2023): minus the mean relationship
# of the selected entries, diagonal included (TrainSel maximises)
avg_grm_self <- function(soln, Data) {
  -mean(Data[["G"]][soln, soln])
}


# -------------------------------- cdmean_prep ---------------------------------
# TrainSel criterion: CDmean of the p-rep design given by soln (order of the plots);
# -(repeats) if a replicated entry appears twice in the same half
cdmean_prep <- function(soln, Data) {
  plot_names <- Data$items[soln]
  n_dup <- sum(duplicated(plot_names[Data$half == 1])) +
    sum(duplicated(plot_names[Data$half == 2]))
  if (n_dup > 0) {
    return(-n_dup)
  }

  Z <- model.matrix(~ factor(plot_names, levels = Data$genos) - 1)
  C12 <- crossprod(Data$RinvX, Z)
  C22 <- crossprod(Z, Data$Rinv %*% Z) + Data$Ginv
  C <- rbind(cbind(Data$tXRinvX, C12), cbind(t(C12), C22))
  PEV <- diag(chol2inv(chol(C)))[-seq_len(ncol(Data$tXRinvX))]
  mean((diag(Data$G) - PEV) / diag(Data$G))
}
