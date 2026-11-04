# Script removes the 20 oldest parents each year and replaces them
# with the 20 genotypes of the S4 stage to form a new crossing block
cat(" Update parents \n")

# Drop 20 parents
Parents <- Parents[21:nParents]

# Update with the 20 new parents from the S4 stage
Parents <- c(Parents, S4)
