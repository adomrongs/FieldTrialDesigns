# Create founders
cat(" Simulating founder parents \n")

# Generate initial haplotypes
founderPop <- runMacs2(
  nInd = nParents,
  nChr = 12,
  segSites = nQtl + nSnp,
  inbred = FALSE,
  ploidy = 4L
)
SP <- SimParam$new(founderPop)
# QTL and SNPs on different sites: the SNPs only capture the QTL through LD
SP$restrSegSites(minQtlPerChr = nQtl, minSnpPerChr = nSnp, overlap = FALSE)

# Create an input object with FieldSimR for AlphaSimR
# to obtain desired genotype slopes
input_asr <- multi_asr_input(
  nenvs = p,
  mean = 0,
  var = diag(De),
  corA = Ce,
  nterms = k
)

# Add additive trait
SP$addTraitA(
  nQtlPerChr = nQtl,
  mean = input_asr$mean,
  var = input_asr$var,
  corA = input_asr$corA
)

# SNP chip used for every genomic relationship matrix (GS and odw)
SP$addSnpChip(nSnpPerChr = nSnp)

# Create founder parents
Parents <- newPop(founderPop)

#-------------
rm(founderPop)
