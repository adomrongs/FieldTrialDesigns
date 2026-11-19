# Script to obtain the genetic merit (GEBV) of the S2, S3 and S4 individuals.
# One-stage model, as in Werner et al. (2025), with spatial terms per trial:
#   y = trial + id (GRM) + row(trial) + col(trial) + e,   e ~ N(0, diag(sigma2_e per trial))
# trial = stage_year_loc (fixed). The residual of each trial also holds its GxE.

cat(" Fit genomic model \n")

# ------- Build GRM and its inverse -------
geno_pop <- mergePops(list(train_data$pop, S2[!(S2@id %in% train_data$pop@id)]))
M <- pullSnpGeno(geno_pop)
rownames(M) <- geno_pop@id
G <- vanraden_grm(M)

id_train <- unique(train_data$pheno$id)
id_cand <- S2@id
G_train <- G[id_train, id_train]
Gi <- solve(G_train)
attr(Gi, "inverse") <- TRUE

# ------- Prepare data -------
d_fit <- train_data$pheno |>
  mutate(
    trial = factor(env),
    stage = factor(stage),
    id = factor(id, levels = id_train),
    trial_row = factor(paste(env, row, sep = "_")),
    trial_col = factor(paste(env, col, sep = "_")),
    resid_group = if (residualBy == "trial") trial else stage
  )

# ------- One-stage model -------
fit <- mmes(
  y.Trait1 ~ 1,
  random = ~ vsm(ism(id), Gu = Gi) +
    vsm(ism(trial)) +
    vsm(ism(trial_row)) +
    vsm(ism(trial_col)),
  rcov = ~ vsm(dsm(resid_group), ism(units)),
  data = d_fit,
  verbose = FALSE
)
u_train <- fit$uList[[1]][id_train, 1]

# ------- GEBVs -------
# S2 (not phenotyped): predicted from the training individuals through G
gebv <- G[id_cand, id_train] %*% solve(G_train, u_train)
S2@ebv <- matrix(as.numeric(gebv), ncol = 1)

# S3 and S4 of last year (in the training data): used to select S4 and S5
S3@ebv <- matrix(as.numeric(u_train[S3@id]), ncol = 1)
S4@ebv <- matrix(as.numeric(u_train[S4@id]), ncol = 1)
