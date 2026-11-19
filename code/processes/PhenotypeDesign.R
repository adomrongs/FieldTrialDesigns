# Phenotype a trial in one or more sites whose layout is already given (odw
# design), instead of letting FieldSimR randomise it as phenotype_met() does.
# Input:  design_pop    pop with every id of the layout
#         design_layout data.frame site, row, col, id (one row per plot)
#         design_envs   environment of each site (rows of covs_tpe), in site order
#         design_varR, design_ncols, design_nrows (field of one site), covs_tpe
# Output: res = list(pop, pheno), as phenotype_met(return.pheno = TRUE), so it can
#         go straight to store_training() (env = site); pop@pheno = mean of the
#         plots of each id over all sites

n_sites_design <- length(design_envs)

# Genetic values in the environment of each site (multi_asr_output needs
# nenvs > 1: auxiliary environment with a single site, dropped afterwards)
envs_calc <- if (n_sites_design == 1) {
  c(design_envs, setdiff(seq_len(nrow(covs_tpe)), design_envs)[1])
} else {
  design_envs
}
gv_design <- multi_asr_output(
  pop = design_pop,
  ntraits = 1,
  nenvs = length(envs_calc),
  nreps = 1,
  cov.mat = covs_tpe[envs_calc, , drop = FALSE]
) |>
  mutate(
    site = as.integer(as.character(env)),
    id = sprintf("%.0f", as.numeric(as.character(id)))
  ) |>
  filter(site <= n_sites_design) |>
  select(site, id, gv.Trait1)

# Spatially correlated plot errors of each site (one block: the design is in the
# layout, not in FieldSimR), independent between sites, with the AR1
# correlations of the odw design model
err_design <- field_trial_error(
  nenvs = n_sites_design,
  nblocks = 1,
  varR = design_varR,
  spatial.model = "AR1",
  col.cor = corColS3,
  row.cor = corRowS3,
  ncols = design_ncols,
  nrows = design_nrows
) |>
  mutate(
    site = as.integer(as.character(env)),
    col = as.numeric(as.character(col)),
    row = as.numeric(as.character(row))
  ) |>
  select(site, col, row, e.Trait1)

# Plot phenotypes: y = gv + plot error. Same column types as the FieldSimR output
# (factors) so that bind_rows() in store_training() works with the other stages
pheno_design <- design_layout |>
  left_join(gv_design, by = c("site", "id")) |>
  left_join(err_design, by = c("site", "col", "row")) |>
  mutate(
    y.Trait1 = gv.Trait1 + e.Trait1,
    env = factor(site),
    block = factor(1),
    col = factor(col),
    row = factor(row)
  )

if (anyNA(pheno_design$y.Trait1)) {
  stop(
    "PhenotypeDesign.R: ",
    sum(is.na(pheno_design$y.Trait1)),
    " plots without phenotype (check ids and grid dimensions)."
  )
}

# Mean of the plots of each individual, in the order of design_pop@id
means_design <- pheno_design |>
  group_by(id) |>
  summarise(mean_pheno = mean(y.Trait1), .groups = "drop")
design_pop@pheno[, 1] <- means_design$mean_pheno[match(
  design_pop@id,
  means_design$id
)]

res <- list(pop = design_pop, pheno = pheno_design)
