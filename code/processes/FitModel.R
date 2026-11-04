# Script to obtain the genetic merit of S2 individuals

cat(" Fit genomic model \n")

# ------- Stage1 -------
# Run model per environment
blues <- NULL
vcov_blues <- list()
for (e in unique(train_data$pheno$env)) {
  d <- train_data$pheno |>
    filter(env == e) |>
    mutate(
      id = factor(id),
      block = factor(block),
      col = as.numeric(as.character(col)),
      row = as.numeric(as.character(row))
    )

  fit1 <- SpATS(
    response = "y.Trait1",
    genotype = "id",
    genotype.as.random = FALSE,
    spatial = ~ SAP(col, row),
    random = ~block,
    data = d,
    control = list(monitoring = 0)
  )
  pred <- predict(
    fit1,
    which = "id",
    predFixed = "marginal",
    return.vcov.matrix = TRUE
  )

  blues <- bind_rows(
    blues,
    data.frame(
      env = e,
      loc = unique(d$loc),
      id = as.character(pred$id),
      BLUE = pred$predicted.values,
      se = pred$standard.errors
    )
  )
  vcov_blues[[e]] <- as.matrix(attr(pred, "vcov"))
}

# Build GRM and its inverse
geno_pop <- mergePops(list(train_data$pop, S2[!(S2@id %in% train_data$pop@id)]))
M <- pullSnpGeno(geno_pop)
rownames(M) <- geno_pop@id
var_marker <- colMeans(M^2) - colMeans(M)^2
M <- M[, var_marker > 1e-8, drop = FALSE]
freq <- colMeans(M) / 4
flip <- freq > 0.5 # code every marker by its minor allele
M[, flip] <- 4 - M[, flip]
freq[flip] <- 1 - freq[flip]
W <- sweep(M, 2, 4 * freq)
G <- tcrossprod(W) / sum(4 * freq * (1 - freq))
G <- G + diag(1e-3, nrow(G)) # make it invertible

id_train <- unique(blues$id)
id_cand <- S2@id
G_train <- G[id_train, id_train]
Gi <- solve(G_train) # G already has 1e-3 added to the diagonal
attr(Gi, "inverse") <- TRUE

# Build inputs for stage2
blues <- blues |>
  mutate(
    env.id = paste(env, id, sep = "_"),
    env = factor(env),
    loc = factor(loc),
    id = factor(id, levels = id_train),
    env.id = factor(env.id, levels = env.id)
  )
Oi <- Matrix::bdiag(lapply(vcov_blues[as.character(unique(blues$env))], solve))
dimnames(Oi) <- list(levels(blues$env.id), levels(blues$env.id))
attr(Oi, "inverse") <- TRUE

# ------- Stage2 -------
fit2 <- mmes(
  BLUE ~ env,
  random = ~ vsm(ism(id), Gu = Gi) +
    vsm(ism(env.id), Gu = Oi, sigma2 = 1, fixedSigma2 = TRUE),
  rcov = ~units,
  data = blues
)
u_train <- fit2$uList[[1]][id_train, 1]

# Extract GEBVs for the test set and compute its accuracy
gebv <- G[id_cand, id_train] %*% solve(G_train, u_train)
S2@ebv <- matrix(as.numeric(gebv), ncol = 1)
gv_tpe_S2 <- as.numeric(
  S2@gv[, 1] + S2@gv[, -1, drop = FALSE] %*% colMeans(covs_tpe)
)
acc_gs <- safe_cor(S2@ebv[, 1], gv_tpe_S2)

# GEBVs of the S3 entries (trial of last year, used to select S4)
S3@ebv <- matrix(as.numeric(u_train[S3@id]), ncol = 1)
