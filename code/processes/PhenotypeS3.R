# Script to phenotype the S3 trial with the design of the scenario.

if (S3design == "RCBD") {
  res <- phenotype_met(
    pop = S3,
    envs = met_envs(samples_met, simYear, nEnvS3),
    cov.mat = covs_tpe,
    nreps = nRepS3,
    varR = varES3,
    ncols = nColS3,
    nrows = nRowS3,
    return.pheno = TRUE
  )
}

if (S3design %in% c("pREP_odw", "Augmented")) {
  # Add checks if needed
  design_pop <- if (S3design == "Augmented") {
    mergePops(list(S3, checks_S3))
  } else {
    S3
  }
  # Build GRM for S3 pop
  M_s3 <- pullSnpGeno(design_pop)
  rownames(M_s3) <- design_pop@id
  var_marker_s3 <- colMeans(M_s3^2) - colMeans(M_s3)^2
  M_s3 <- M_s3[, var_marker_s3 > 1e-8, drop = FALSE]
  freq_s3 <- colMeans(M_s3) / 4
  flip_s3 <- freq_s3 > 0.5
  M_s3[, flip_s3] <- 4 - M_s3[, flip_s3]
  freq_s3[flip_s3] <- 1 - freq_s3[flip_s3]
  W_s3 <- sweep(M_s3, 2, 4 * freq_s3)
  G_s3 <- tcrossprod(W_s3) / sum(4 * freq_s3 * (1 - freq_s3))
  G_s3 <- G_s3 / mean(diag(G_s3)) + diag(1e-3, nrow(G_s3))
  Ginv_s3 <- solve(G_s3)
  attr(Ginv_s3, "INVERSE") <- TRUE

  # Fix odw GRM input
  ginv_stage2 <- paste("Ginv", S3design, "stage2", Rep, year, sep = "_")
  ginv_stage3 <- paste("Ginv", S3design, "stage3", Rep, year, sep = "_")
  assign(ginv_stage2, Ginv_s3)
  assign(ginv_stage3, Ginv_s3)
  genetic_stage2 <- as.formula(paste0("~ vm(name, ", ginv_stage2, ")"))
  genetic_stage3 <- as.formula(paste0("~ vm(name, ", ginv_stage3, ")"))
}

if (S3design == "pREP_odw") {
  # Create an initial data.frame
  init_stage2 <- data.frame(
    name = factor(S3@id, levels = rownames(G_s3))
  ) |>
    mutate(repF = factor(if_else(name %in% sample(name, nRepS3pREP), 2, 1))) |>
    arrange(repF)

  # Fix variance components
  sv_stage2 <- odw(
    fixed = ~1,
    random = genetic_stage2,
    residual = ~ dsum(~ units | repF),
    permute = genetic_stage2,
    data = init_stage2,
    start.values = TRUE
  )$vparameters.table
  sv_stage2$Value <- c(gammaS3, 1 / 1, 1 / 2)

  # Allocate replication to genotypes
  log_odw <- capture.output(
    stage2_s3 <- odw(
      fixed = ~1,
      random = genetic_stage2,
      residual = ~ dsum(~ units | repF),
      permute = genetic_stage2,
      G.param = sv_stage2,
      R.param = sv_stage2,
      search = "tabu+rw",
      maxit = 15,
      data = init_stage2
    )
  )
  ids_rep <- as.character(stage2_s3$design$name[stage2_s3$design$repF == 2])
  ids_unrep <- sample(setdiff(S3@id, ids_rep))

  # Prepare initial data.frame
  genotypes_s3 <- bind_rows(
    data.frame(name = rep(ids_rep, 2), block = rep(1:2, each = nRepS3pREP)),
    data.frame(name = ids_unrep, block = rep(1:2, each = length(ids_unrep) / 2))
  ) |>
    group_by(block) |>
    slice_sample(prop = 1) |>
    ungroup() |>
    arrange(block)
  init_stage3 <- expand.grid(
    row = seq_len(nRowS3pREP),
    col = seq_len(nColS3pREP)
  ) |>
    mutate(block = if_else(row <= nRowS3pREP / 2, 1, 2)) |>
    arrange(block) |>
    mutate(
      name = genotypes_s3$name,
      swp = if_else(name %in% ids_rep, paste0("R", block), "U") # step 1 groups
    )
  design_pop <- S3
  design_ncols <- nColS3pREP
  design_nrows <- nRowS3pREP
}

if (S3design == "Augmented") {
  genotypes_s3 <- bind_rows(
    data.frame(
      name = rep(checks_S3@id, nBlocksAugS3),
      block = rep(seq_len(nBlocksAugS3), each = nChecksS3)
    ),
    data.frame(
      name = sample(S3@id),
      block = rep(seq_len(nBlocksAugS3), each = nS3 / nBlocksAugS3)
    )
  ) |>
    group_by(block) |>
    slice_sample(prop = 1) |>
    ungroup() |>
    arrange(block)

  init_stage3 <- expand.grid(
    row = seq_len(nRowS3Aug),
    col = seq_len(nColS3Aug)
  ) |>
    mutate(
      block = (ceiling(row / (nRowS3Aug / 2)) - 1) *
        5 +
        ceiling(col / (nColS3Aug / 5))
    ) |>
    arrange(block) |>
    mutate(
      name = genotypes_s3$name,
      swp = if_else(name %in% checks_S3@id, paste0("C", block), "E") # step 1 groups
    )
  design_ncols <- nColS3Aug
  design_nrows <- nRowS3Aug
}

if (S3design %in% c("pREP_odw", "Augmented")) {
  # Prepare data.frame for modelling
  init_stage3 <- init_stage3 |>
    mutate(
      across(c(row, col, block, swp), as.factor),
      name = factor(name, levels = rownames(G_s3))
    ) |>
    arrange(col, row)

  # Fix variance components
  sv_stage3 <- odw(
    fixed = ~block,
    random = genetic_stage3,
    residual = ~ ar1(col):ar1(row),
    permute = genetic_stage3,
    swap = ~block,
    data = init_stage3,
    start.values = TRUE
  )$vparameters.table
  sv_stage3$Value <- c(gammaS3, 1, corColS3, corRowS3)

  # Allocate genotypes to plot (First step)
  log_odw <- capture.output(
    step1_s3 <- odw(
      fixed = ~block,
      random = genetic_stage3,
      residual = ~ ar1(col):ar1(row),
      permute = genetic_stage3,
      swap = ~swp,
      G.param = sv_stage3,
      R.param = sv_stage3,
      search = "tabu+rw",
      maxit = 10,
      data = init_stage3
    )
  )

  # Allocate genotypes to plot (Second step)
  log_odw <- capture.output(
    step2_s3 <- odw(
      fixed = ~block,
      random = genetic_stage3,
      residual = ~ ar1(col):ar1(row),
      permute = genetic_stage3,
      swap = ~block,
      G.param = sv_stage3,
      R.param = sv_stage3,
      search = "tabu+rw",
      maxit = 100,
      data = step1_s3$design
    )
  )

  # Extract final layout
  design_layout <- step2_s3$design |>
    transmute(
      row = as.numeric(as.character(row)),
      col = as.numeric(as.character(col)),
      block = as.numeric(as.character(block)),
      id = as.character(name)
    )
  design_envs <- met_envs(samples_met, simYear, nEnvS3)
  design_varR <- varES3
  source("code/processes/PhenotypeDesign.R")

  rm(list = c(ginv_stage2, ginv_stage3)) # free the copies of this year
}
