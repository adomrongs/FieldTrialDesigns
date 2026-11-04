# Phenotype a single-environment trial whose layout is already given (an odw
# design), instead of letting FieldSimR randomise it as phenotype_met() does.
# Input:  design_pop    pop with every id of the layout
#         design_layout data.frame row, col, block, id (one row per plot)
#         design_envs   environment of the trial (one row of covs_tpe)
#         design_varR, design_ncols, design_nrows, covs_tpe
# Output: res = list(pop, pheno), as phenotype_met(return.pheno = TRUE), so it can
#         go straight to store_training(); pop@pheno = mean of the plots of each id

# Genetic values in the trial environment (multi_asr_output needs nenvs > 1:
# auxiliary environment, dropped afterwards, as in phenotype_met)
envs_calc <- c(design_envs, setdiff(seq_len(nrow(covs_tpe)), design_envs)[1])
gv_design <- multi_asr_output(
  pop = design_pop,
  ntraits = 1,
  nenvs = length(envs_calc),
  nreps = 1,
  cov.mat = covs_tpe[envs_calc, , drop = FALSE]
) |>
  filter(as.integer(as.character(env)) == 1) |>
  mutate(id = sprintf("%.0f", as.numeric(as.character(id)))) |>
  select(id, gv.Trait1)

# Spatially correlated plot errors on the whole grid (one block: the design is in
# the layout, not in FieldSimR), with the AR1 correlations of the odw design model
err_design <- field_trial_error(
  nenvs = 1,
  nblocks = 1,
  varR = design_varR,
  spatial.model = "AR1",
  col.cor = corColS3,
  row.cor = corRowS3,
  ncols = design_ncols,
  nrows = design_nrows
) |>
  mutate(
    col = as.numeric(as.character(col)),
    row = as.numeric(as.character(row))
  ) |>
  select(col, row, e.Trait1)

# Plot phenotypes: y = gv + plot error. Same column types as the FieldSimR output
# (factors) so that bind_rows() in store_training() works with the other stages
pheno_design <- design_layout |>
  left_join(gv_design, by = "id") |>
  left_join(err_design, by = c("col", "row")) |>
  mutate(
    y.Trait1 = gv.Trait1 + e.Trait1,
    env = factor(1),
    block = factor(block),
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
