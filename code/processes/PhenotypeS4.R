# Script to phenotype the S4 trial with the design of the scenario (S4design):
#   "RCBD" all entries in each of the nSitesS4 sites, nRepS4 blocks per site
#          (randomised by FieldSimR)
#   "pREP" the nS4 entries in each of the nSitesS4 sites + nRepS4pREP replicated
#          entries per site (S4repAlloc: "rotating" = different ones in each site,
#          "fixed_random" / "fixed_GRM" = the same ones in every site); positions
#          optimised with odw in each site (DesignPrepSites.R)
# Output: res = list(pop, pheno), ready for store_training()

if (S4design == "RCBD") {
  res <- phenotype_met(
    pop = S4,
    envs = met_envs(samples_met, simYear, nSitesS4),
    cov.mat = covs_tpe,
    nreps = nRepS4,
    varR = varES4,
    ncols = nColS4,
    nrows = nRowS4,
    return.pheno = TRUE
  )
}

if (S4design == "pREP") {
  ids_s4 <- sample(S4@id)

  # Replicated entries of each site
  if (S4repAlloc == "rotating") {
    # consecutive blocks of the shuffled entries: each entry replicated in a
    # different site (with 5 sites, 5 entries twice, in different sites)
    slots <- (seq_len(nSitesS4 * nRepS4pREP) - 1) %% length(ids_s4) + 1
    prep_reps <- split(ids_s4[slots], rep(seq_len(nSitesS4), each = nRepS4pREP))
  }
  if (S4repAlloc == "fixed_random") {
    prep_reps <- rep(list(ids_s4[seq_len(nRepS4pREP)]), nSitesS4)
  }
  if (S4repAlloc == "fixed_GRM") {
    # odw chooses the replicated entries with the GRM: unreplicated entries have
    # nSitesS4 plots, replicated entries 2 x nSitesS4 (residual of their means)
    M_rep <- pullSnpGeno(S4)
    rownames(M_rep) <- S4@id
    Ginv_rep <- solve(vanraden_grm(M_rep, scale_diag = TRUE))
    attr(Ginv_rep, "INVERSE") <- TRUE
    ginv_rep <- paste("Ginv_S4_rep", Rep, year, sep = "_")
    assign(ginv_rep, Ginv_rep)
    genetic_rep <- as.formula(paste0("~ vm(name, ", ginv_rep, ")"))
    init_rep <- data.frame(name = factor(ids_s4, levels = S4@id)) |>
      mutate(
        repF = factor(if_else(name %in% ids_s4[seq_len(nRepS4pREP)], "2", "1"))
      ) |>
      arrange(repF)
    sv_rep <- odw(
      fixed = ~1,
      random = genetic_rep,
      residual = ~ dsum(~ units | repF),
      permute = genetic_rep,
      data = init_rep,
      start.values = TRUE
    )$vparameters.table
    sv_rep$Value[grepl("^vm", sv_rep$Component)] <- gammaS4
    sv_rep$Value[grepl("repF_1!R", sv_rep$Component)] <- 1 / nSitesS4
    sv_rep$Value[grepl("repF_2!R", sv_rep$Component)] <- 1 / (2 * nSitesS4)
    log_odw <- capture.output(
      rep_s4 <- odw(
        fixed = ~1,
        random = genetic_rep,
        residual = ~ dsum(~ units | repF),
        permute = genetic_rep,
        G.param = sv_rep,
        R.param = sv_rep,
        search = "tabu+rw",
        maxit = 15,
        data = init_rep
      )
    )
    ids_rep <- as.character(rep_s4$design$name[rep_s4$design$repF == "2"])
    prep_reps <- rep(list(ids_rep), nSitesS4)
    rm(list = ginv_rep)
  }

  # Positions in each site (odw) and phenotypes
  prep_pop <- S4
  prep_nrow <- nRowS4pREP
  prep_ncol <- nColS4pREP
  prep_gamma <- gammaS4
  prep_stage <- "S4"
  source("code/processes/DesignPrepSites.R")

  design_pop <- S4
  design_envs <- met_envs(samples_met, simYear, nSitesS4)
  design_varR <- varES4
  design_ncols <- nColS4pREP
  design_nrows <- nRowS4pREP
  source("code/processes/PhenotypeDesign.R")
}
