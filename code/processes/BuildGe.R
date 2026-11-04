
# Create diagonal variance matrix
De <- skew_diag_mat(
  n = p,
  shape = 1.5,
  scale = 1
)

# Create correlation matrix
Ce <- struc_cor_mat(
  n = p,
  base.cor = rho,
  skew = gamma,
  rank = k
)

# Generate Ge
Ge <- sqrt(De) %*% Ce %*% sqrt(De)
