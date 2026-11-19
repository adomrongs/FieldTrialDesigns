# Timing test of odw for the S3 sparse testing designs and the S4 p-rep per site.
# Standalone (outside the repo): simulates a tetraploid GRM with family structure,
# then times (1) the choice of common entries and (2) the allocation of the
# unique entries to sites, for the largest and smallest S3 scenarios.
# Run on the cluster: module load r-bio; Rscript odw_sparse_timing.R

# ------- Prepare Data -------
# Load packages
library(odw)
library(dplyr)

Sys.setenv(OMP_NUM_THREADS = 1)
set.seed(2026)

# Design parameters (same as the simulation)
gamma <- 0.5 # sigma2_g / sigma2_e per plot
maxit_test <- 10 # iterations per odw call (time scales roughly linearly)

# S3 scenarios to time: sites, plots, entries, common entries
scenarios <- data.frame(
  scenario = c("ES2-ST", "I1", "I3", "I4"),
  n_sites = c(6, 6, 6, 3),
  n_plots = c(810, 660, 360, 210),
  n_entries = c(540, 440, 240, 176),
  n_common = c(54, 44, 24, 17)
)

# ------- Simulate a tetraploid GRM with family structure -------
# 50 parents, 100 crosses, progeny = two random tetraploid gametes (2 of 4 alleles)
n_snp <- 6000
n_par <- 50
parents <- matrix(rbinom(n_par * n_snp, 4, runif(n_snp, 0.05, 0.95)),
  nrow = n_par, byrow = TRUE
)
crosses <- t(replicate(100, sample(n_par, 2)))
max_entries <- max(scenarios$n_entries)
cross_id <- rep(seq_len(100), length.out = max_entries)
M <- t(sapply(cross_id, function(k) {
  p1 <- parents[crosses[k, 1], ]
  p2 <- parents[crosses[k, 2], ]
  rhyper(n_snp, p1, 4 - p1, 2) + rhyper(n_snp, p2, 4 - p2, 2)
}))
rownames(M) <- paste0("G", seq_len(max_entries))
M <- M[, apply(M, 2, var) > 0]
freq <- colMeans(M) / 4
W <- sweep(M, 2, 4 * freq)
G_all <- tcrossprod(W) / sum(4 * freq * (1 - freq))
G_all <- G_all / mean(diag(G_all)) + diag(1e-3, nrow(G_all))

# ------- Run timing test -------
results <- NULL
for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  cat("\n=====", sc$scenario, "=====\n")
  ids <- rownames(G_all)[seq_len(sc$n_entries)]
  Ginv <- solve(G_all[ids, ids])
  attr(Ginv, "INVERSE") <- TRUE

  # odw keeps information linked to the object name: a new name per call
  ginv_step1 <- paste0("Ginv_step1_", i)
  ginv_step2 <- paste0("Ginv_step2_", i)
  assign(ginv_step1, Ginv)
  assign(ginv_step2, Ginv)
  genetic_step1 <- as.formula(paste0("~ vm(name, ", ginv_step1, ")"))
  genetic_step2 <- as.formula(paste0("~ vm(name, ", ginv_step2, ")"))

  # Step 1: which entries are common (in all sites) -----------------------
  init_rep <- data.frame(name = factor(ids, levels = ids)) |>
    mutate(repF = factor(if_else(name %in% sample(name, sc$n_common), "2", "1"))) |>
    arrange(repF)

  sv1 <- odw(
    fixed = ~1,
    random = genetic_step1,
    residual = ~ dsum(~units | repF),
    permute = genetic_step1,
    data = init_rep,
    start.values = TRUE
  )$vparameters.table
  print(sv1)
  # genetic ratio, unique entries (1 plot), common entries (mean of n_sites plots)
  sv1$Value <- c(gamma, 1, 1 / sc$n_sites)

  t1 <- system.time(
    log1 <- capture.output(
      step1 <- odw(
        fixed = ~1,
        random = genetic_step1,
        residual = ~ dsum(~units | repF),
        permute = genetic_step1,
        G.param = sv1,
        R.param = sv1,
        search = "tabu+rw",
        maxit = maxit_test,
        data = init_rep
      )
    )
  )
  ids_common <- as.character(step1$design$name[step1$design$repF == "2"])
  ids_unique <- sample(setdiff(ids, ids_common))

  # Step 2: allocate the unique entries to sites ---------------------------
  n_unique_site <- length(ids_unique) / sc$n_sites
  init_sites <- bind_rows(lapply(seq_len(sc$n_sites), function(s) {
    uniq <- ids_unique[((s - 1) * n_unique_site + 1):(s * n_unique_site)]
    data.frame(
      site = s,
      name = c(ids_common, uniq),
      swp = c(rep(paste0("C", s), length(ids_common)), rep("U", length(uniq)))
    )
  })) |>
    mutate(
      site = factor(site),
      swp = factor(swp),
      name = factor(name, levels = ids)
    ) |>
    arrange(site)

  sv2 <- odw(
    fixed = ~site,
    random = genetic_step2,
    residual = ~ dsum(~units | site),
    permute = genetic_step2,
    swap = ~swp,
    data = init_sites,
    start.values = TRUE
  )$vparameters.table
  print(sv2)
  sv2$Value <- c(gamma, rep(1, sc$n_sites))

  t2 <- system.time(
    log2 <- capture.output(
      step2 <- odw(
        fixed = ~site,
        random = genetic_step2,
        residual = ~ dsum(~units | site),
        permute = genetic_step2,
        swap = ~swp,
        G.param = sv2,
        R.param = sv2,
        search = "tabu+rw",
        maxit = maxit_test,
        data = init_sites
      )
    )
  )

  # Checks: plots per site, commons once per site, uniques once in total
  chk <- step2$design |>
    count(name) |>
    mutate(type = if_else(name %in% ids_common, "common", "unique"))
  dup_site <- step2$design |>
    count(site, name) |>
    filter(n > 1) |>
    nrow()
  cat("Plots per site:", table(step2$design$site), "\n")
  cat("Commons with", sc$n_sites, "plots:", all(chk$n[chk$type == "common"] == sc$n_sites), "\n")
  cat("Uniques with 1 plot:", all(chk$n[chk$type == "unique"] == 1), "\n")
  cat("Entries repeated within a site:", dup_site, "\n")
  cat("Last lines of the step 2 log:\n")
  cat(tail(log2, 5), sep = "\n")

  results <- bind_rows(
    results,
    data.frame(
      scenario = sc$scenario,
      n_plots = sc$n_plots,
      n_entries = sc$n_entries,
      maxit = maxit_test,
      step1_sec = t1[["elapsed"]],
      step2_sec = t2[["elapsed"]]
    )
  )
  rm(list = c(ginv_step1, ginv_step2))
}

# ------- S4: p-rep in one site (20 entries + 5 replicated, 5 x 5) -------
ids4 <- rownames(G_all)[1:20]
Ginv4 <- solve(G_all[ids4, ids4])
attr(Ginv4, "INVERSE") <- TRUE
genetic_s4 <- ~ vm(name, Ginv4)
init_s4 <- expand.grid(row = factor(1:5), col = factor(1:5)) |>
  mutate(name = factor(sample(c(ids4, sample(ids4, 5))), levels = ids4)) |>
  arrange(col, row)
sv4 <- odw(
  fixed = ~1,
  random = genetic_s4,
  residual = ~ ar1(col):ar1(row),
  permute = genetic_s4,
  data = init_s4,
  start.values = TRUE
)$vparameters.table
print(sv4)
sv4$Value <- c(gamma, 1, 0.5, 0.7)
t4 <- system.time(
  log4 <- capture.output(
    step4 <- odw(
      fixed = ~1,
      random = genetic_s4,
      residual = ~ ar1(col):ar1(row),
      permute = genetic_s4,
      G.param = sv4,
      R.param = sv4,
      search = "tabu+rw",
      maxit = 100,
      data = init_s4
    )
  )
)
results <- bind_rows(
  results,
  data.frame(
    scenario = "S4 one site", n_plots = 25, n_entries = 20, maxit = 100,
    step1_sec = NA, step2_sec = t4[["elapsed"]]
  )
)

# ------- Results -------
print(results)
write.csv(results, "odw_sparse_timing.csv", row.names = FALSE)
