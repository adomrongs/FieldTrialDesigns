# Test: does a spatial component (AR1 x AR1 in each site) make odw worthwhile for
# the S3 sparse testing designs? Compares the mean PEV (avalue) of a random design
# with the odw design, and times each step. Standalone (outside the repo).
# Run on the cluster: module load r-bio; Rscript odw_sparse_spatial.R

# ------- Prepare Data -------
# Load packages
library(odw)
library(dplyr)

Sys.setenv(OMP_NUM_THREADS = 1)
set.seed(2026)

# Design parameters (same as the simulation)
gamma <- 0.5 # sigma2_g / sigma2_e per plot
cor_col <- 0.5 # AR1 correlation between columns (corColS3)
cor_row <- 0.7 # AR1 correlation between rows (corRowS3)
maxit_test <- 10

# S3 scenarios: sites, entries, commons and field of each site (rows x cols)
scenarios <- data.frame(
  scenario = c("I4", "I3", "ES2-ST"),
  n_sites = c(3, 6, 6),
  n_entries = c(176, 240, 540),
  n_common = c(17, 24, 54),
  n_row = c(7, 6, 9), # 7 x 10 = 70, 6 x 10 = 60, 9 x 15 = 135 plots per site
  n_col = c(10, 10, 15)
)

# ------- Simulate a tetraploid GRM with family structure -------
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

# ------- Run test -------
results <- NULL
for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  cat("\n=====", sc$scenario, "=====\n")
  ids <- rownames(G_all)[seq_len(sc$n_entries)]
  Ginv <- solve(G_all[ids, ids])
  attr(Ginv, "INVERSE") <- TRUE

  # A new object name per odw call (odw links information to the name)
  ginv_a <- paste0("Ginv_a_", i)
  ginv_b <- paste0("Ginv_b_", i)
  assign(ginv_a, Ginv)
  assign(ginv_b, Ginv)
  genetic_a <- as.formula(paste0("~ vm(name, ", ginv_a, ")"))
  genetic_b <- as.formula(paste0("~ vm(name, ", ginv_b, ")"))

  # Random initial design: commons in every site, uniques split evenly
  ids_common <- sample(ids, sc$n_common)
  ids_unique <- sample(setdiff(ids, ids_common))
  n_unique_site <- length(ids_unique) / sc$n_sites
  init_sites <- bind_rows(lapply(seq_len(sc$n_sites), function(s) {
    uniq <- ids_unique[((s - 1) * n_unique_site + 1):(s * n_unique_site)]
    expand.grid(row = seq_len(sc$n_row), col = seq_len(sc$n_col)) |>
      mutate(
        site = s,
        name = sample(c(ids_common, uniq)),
        swp = if_else(name %in% ids_common, paste0("C", s), "U")
      )
  })) |>
    mutate(
      across(c(site, row, col, swp), factor),
      name = factor(name, levels = ids)
    ) |>
    arrange(site, col, row)

  # Start values: matched by name (genetic ratio, residual variance, AR1 correlations)
  sv <- odw(
    fixed = ~site,
    random = genetic_a,
    residual = ~ dsum(~ ar1(col):ar1(row) | site),
    permute = genetic_a,
    swap = ~swp,
    data = init_sites,
    start.values = TRUE
  )$vparameters.table
  print(sv)
  sv$Value[grepl("^vm", sv$Component)] <- gamma
  sv$Value[grepl("!R$", sv$Component)] <- 1
  sv$Value[grepl("col!cor$", sv$Component)] <- cor_col
  sv$Value[grepl("row!cor$", sv$Component)] <- cor_row
  print(sv)

  # Step 2a: allocate the uniques to sites (spatial model) ----------------
  t_a <- system.time(
    log_a <- capture.output(
      step_a <- odw(
        fixed = ~site,
        random = genetic_a,
        residual = ~ dsum(~ ar1(col):ar1(row) | site),
        permute = genetic_a,
        swap = ~swp,
        G.param = sv,
        R.param = sv,
        search = "tabu+rw",
        maxit = maxit_test,
        data = init_sites
      )
    )
  )

  # Step 2b: positions within each site ------------------------------------
  t_b <- system.time(
    log_b <- capture.output(
      step_b <- odw(
        fixed = ~site,
        random = genetic_b,
        residual = ~ dsum(~ ar1(col):ar1(row) | site),
        permute = genetic_b,
        swap = ~site,
        G.param = sv,
        R.param = sv,
        search = "tabu+rw",
        maxit = maxit_test,
        data = step_a$design |> arrange(site, col, row)
      )
    )
  )

  # Checks and mean PEV (avalue): random start, after 2a, after 2b
  dup_site <- step_b$design |>
    count(site, name) |>
    filter(n > 1) |>
    nrow()
  cat("Entries repeated within a site:", dup_site, "\n")
  av_a <- log_a[grepl("avalue", log_a)]
  av_b <- log_b[grepl("avalue", log_b)]
  cat("Step 2a:\n")
  cat(head(av_a, 2), tail(av_a, 1), sep = "\n")
  cat("Step 2b:\n")
  cat(head(av_b, 2), tail(av_b, 1), sep = "\n")
  av_start <- as.numeric(sub(".*= ([0-9.]+).*", "\\1", av_a[1]))
  av_end <- as.numeric(sub(".*: ([0-9.]+).*", "\\1", tail(av_b, 1)))

  results <- bind_rows(
    results,
    data.frame(
      scenario = sc$scenario,
      n_plots = nrow(init_sites),
      step_a_sec = t_a[["elapsed"]],
      step_b_sec = t_b[["elapsed"]],
      avalue_random = av_start,
      avalue_odw = av_end,
      improvement_pct = 100 * (av_start - av_end) / av_start
    )
  )
  rm(list = c(ginv_a, ginv_b))
}

# ------- Results -------
print(results)
write.csv(results, "odw_sparse_spatial.csv", row.names = FALSE)
