# p-rep design in one or more sites: every entry once in each site plus the
# replicated entries of that site; odw optimises the positions within each
# site (AR1 x AR1 per site, no correlation between sites).
# Input:  prep_pop    pop with the entries
#         prep_reps   list (one element per site) with the ids replicated in it
#         prep_nrow, prep_ncol  field of one site (rows x cols = entries + replicated)
#         prep_gamma  sigma2_g / sigma2_e per plot of the design model
#         prep_stage  "S3" or "S4" (only for the name of the Ginv object)
# Output: design_layout data.frame site, row, col, id (one row per plot)

prep_ids <- prep_pop@id
prep_sites <- length(prep_reps)
if (length(prep_ids) + length(prep_reps[[1]]) != prep_nrow * prep_ncol) {
  stop("DesignPrepSites.R: plots per site do not match prep_nrow x prep_ncol")
}

# GRM of the entries (design model, mean diagonal 1)
M_prep <- pullSnpGeno(prep_pop)
rownames(M_prep) <- prep_ids
Ginv_prep <- solve(vanraden_grm(M_prep, scale_diag = TRUE))
attr(Ginv_prep, "INVERSE") <- TRUE

# odw keeps information linked to the object name: a new name per call
ginv_prep <- paste("Ginv", prep_stage, "prep", Rep, year, sep = "_")
assign(ginv_prep, Ginv_prep)
genetic_prep <- as.formula(paste0("~ vm(name, ", ginv_prep, ")"))

# Random initial design
init_prep <- bind_rows(lapply(seq_len(prep_sites), function(s) {
  expand.grid(row = seq_len(prep_nrow), col = seq_len(prep_ncol)) |>
    mutate(site = s, name = sample(c(prep_ids, prep_reps[[s]])))
})) |>
  mutate(
    across(c(site, row, col), factor),
    name = factor(name, levels = prep_ids)
  ) |>
  arrange(site, col, row)

if (prep_sites == 1) {
  sv_prep <- odw(
    fixed = ~1,
    random = genetic_prep,
    residual = ~ ar1(col):ar1(row),
    permute = genetic_prep,
    data = init_prep,
    start.values = TRUE
  )$vparameters.table
} else {
  sv_prep <- odw(
    fixed = ~site,
    random = genetic_prep,
    residual = ~ dsum(~ ar1(col):ar1(row) | site),
    permute = genetic_prep,
    swap = ~site,
    data = init_prep,
    start.values = TRUE
  )$vparameters.table
}
sv_prep$Value[grepl("^vm", sv_prep$Component)] <- prep_gamma
sv_prep$Value[grepl("!R$", sv_prep$Component)] <- 1
sv_prep$Value[grepl("col!cor$", sv_prep$Component)] <- corColS3
sv_prep$Value[grepl("row!cor$", sv_prep$Component)] <- corRowS3

# Positions within each site
if (prep_sites == 1) {
  log_odw <- capture.output(
    prep_design <- odw(
      fixed = ~1,
      random = genetic_prep,
      residual = ~ ar1(col):ar1(row),
      permute = genetic_prep,
      G.param = sv_prep,
      R.param = sv_prep,
      search = "tabu+rw",
      maxit = 100,
      data = init_prep
    )
  )
} else {
  log_odw <- capture.output(
    prep_design <- odw(
      fixed = ~site,
      random = genetic_prep,
      residual = ~ dsum(~ ar1(col):ar1(row) | site),
      permute = genetic_prep,
      swap = ~site,
      G.param = sv_prep,
      R.param = sv_prep,
      search = "tabu+rw",
      maxit = 100,
      data = init_prep
    )
  )
}

design_layout <- prep_design$design |>
  transmute(
    site = as.numeric(as.character(site)),
    row = as.numeric(as.character(row)),
    col = as.numeric(as.character(col)),
    id = as.character(name)
  )
rm(list = ginv_prep) # free the copy of this year
