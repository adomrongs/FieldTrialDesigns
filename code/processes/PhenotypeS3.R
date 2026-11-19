# Script to phenotype the S3 trial with the design of the scenario (S3design):
#   "RCBD"   all entries in each of the nSitesS3 sites, nRepS3 blocks per site
#            (randomised by FieldSimR)
#   "sparse" sparse testing in nSitesS3 sites: nCommonS3 entries (random) in every
#            site and the rest once in a single site, no replication within a site;
#            sites and positions optimised with odw (AR1 x AR1 in each site)
#   "pREP"   one site, nRepS3pREP random entries with a second plot; positions
#            optimised with odw (DesignPrepSites.R)
# Output: res = list(pop, pheno), ready for store_training()

if (S3design == "RCBD") {
  res <- phenotype_met(
    pop = S3,
    envs = met_envs(samples_met, simYear, nSitesS3),
    cov.mat = covs_tpe,
    nreps = nRepS3,
    varR = varES3,
    ncols = nColS3,
    nrows = nRowS3,
    return.pheno = TRUE
  )
}

if (S3design == "sparse") {
  # GRM of the S3 entries (design model, mean diagonal 1)
  M_s3 <- pullSnpGeno(S3)
  rownames(M_s3) <- S3@id
  G_s3 <- vanraden_grm(M_s3, scale_diag = TRUE)
  Ginv_s3 <- solve(G_s3)
  attr(Ginv_s3, "INVERSE") <- TRUE

  # odw keeps information linked to the object name: a new name per call
  ginv_a <- paste("Ginv_S3_sparse_a", Rep, year, sep = "_")
  ginv_b <- paste("Ginv_S3_sparse_b", Rep, year, sep = "_")
  assign(ginv_a, Ginv_s3)
  assign(ginv_b, Ginv_s3)
  genetic_a <- as.formula(paste0("~ vm(name, ", ginv_a, ")"))
  genetic_b <- as.formula(paste0("~ vm(name, ", ginv_b, ")"))

  # Random initial design: commons in every site, uniques split evenly
  ids_common <- sample(S3@id, nCommonS3)
  ids_unique <- sample(setdiff(S3@id, ids_common))
  n_unique_site <- length(ids_unique) / nSitesS3
  if (nCommonS3 + n_unique_site != nRowS3 * nColS3) {
    stop("PhenotypeS3.R: plots per site do not match nRowS3 x nColS3")
  }
  init_s3 <- bind_rows(lapply(seq_len(nSitesS3), function(s) {
    uniq <- ids_unique[((s - 1) * n_unique_site + 1):(s * n_unique_site)]
    expand.grid(row = seq_len(nRowS3), col = seq_len(nColS3)) |>
      mutate(
        site = s,
        name = sample(c(ids_common, uniq)),
        swp = if_else(name %in% ids_common, paste0("C", s), "U") # step 2a groups
      )
  })) |>
    mutate(
      across(c(site, row, col, swp), factor),
      name = factor(name, levels = S3@id)
    ) |>
    arrange(site, col, row)

  # Fix variance components (matched by name)
  sv_a <- odw(
    fixed = ~site,
    random = genetic_a,
    residual = ~ dsum(~ ar1(col):ar1(row) | site),
    permute = genetic_a,
    swap = ~swp,
    data = init_s3,
    start.values = TRUE
  )$vparameters.table
  sv_a$Value[grepl("^vm", sv_a$Component)] <- gammaS3
  sv_a$Value[grepl("!R$", sv_a$Component)] <- 1
  sv_a$Value[grepl("col!cor$", sv_a$Component)] <- corColS3
  sv_a$Value[grepl("row!cor$", sv_a$Component)] <- corRowS3
  sv_b <- sv_a
  sv_b$Component <- sub(ginv_a, ginv_b, sv_b$Component, fixed = TRUE)

  # Step 2a: uniques move between sites, commons within their site
  log_odw <- capture.output(
    step_a_s3 <- odw(
      fixed = ~site,
      random = genetic_a,
      residual = ~ dsum(~ ar1(col):ar1(row) | site),
      permute = genetic_a,
      swap = ~swp,
      G.param = sv_a,
      R.param = sv_a,
      search = "tabu+rw",
      maxit = 10,
      data = init_s3
    )
  )

  # Step 2b: positions within each site
  log_odw <- capture.output(
    step_b_s3 <- odw(
      fixed = ~site,
      random = genetic_b,
      residual = ~ dsum(~ ar1(col):ar1(row) | site),
      permute = genetic_b,
      swap = ~site,
      G.param = sv_b,
      R.param = sv_b,
      search = "tabu+rw",
      maxit = 10,
      data = step_a_s3$design |> arrange(site, col, row)
    )
  )

  # Extract final layout and phenotype it
  design_layout <- step_b_s3$design |>
    transmute(
      site = as.numeric(as.character(site)),
      row = as.numeric(as.character(row)),
      col = as.numeric(as.character(col)),
      id = as.character(name)
    )
  rm(list = c(ginv_a, ginv_b)) # free the copies of this year
}

if (S3design == "pREP") {
  prep_pop <- S3
  prep_reps <- list(sample(S3@id, nRepS3pREP))
  prep_nrow <- nRowS3
  prep_ncol <- nColS3
  prep_gamma <- gammaS3
  prep_stage <- "S3"
  source("code/processes/DesignPrepSites.R")
}

if (S3design %in% c("sparse", "pREP")) {
  design_pop <- S3
  design_envs <- met_envs(samples_met, simYear, nSitesS3)
  design_varR <- varES3
  design_ncols <- nColS3
  design_nrows <- nRowS3
  source("code/processes/PhenotypeDesign.R")
}
