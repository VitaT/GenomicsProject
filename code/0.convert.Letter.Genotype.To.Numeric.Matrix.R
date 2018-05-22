### converting Genotype with letters to numeric matrix ####
#################################################################################
library(data.table)
library(dplyr)

## paths
input <- "../input/"
output <- "../output/"
pic_path <- "../graphs/"

## import genotype data
dAF <- fread(paste0(input, "genotypes_AF"), header = FALSE)
dEA <- fread(paste0(input, "genotypes_EA"), header = FALSE)
dSA <- fread(paste0(input, "genotypes_SA"), header = FALSE)
dWE <- fread(paste0(input, "genotypes_WE"), header = FALSE)
dAM <- fread(paste0(input, "genotypes_AM"), header = FALSE)
dO <- fread(paste0(input, "genotypes_O"), header = FALSE)
dCAS <- fread(paste0(input, "genotypes_CAS"), header = FALSE)

## Snp data
dSnp <- fread(paste0(input, "snps_filtered")) 
colnames(dSnp) <- c("ID", "chrom", "pos", "base", "alt")


#################################################################################
## function to convert the genotypes from letters to numbers
## derived allele is coded 0, ancestrial allele is coded 1
convertHaplotype01 <- function(dSnp, dGenotypes, save_name, ...) {
  dComb <- cbind(dSnp[, 3:5], dGenotypes)
  dConvMatrix <- apply(dComb, 1, function(x) {
    a <- x[3]
    b <- x[2]
    n_A <- gsub(a, 0, x[-c(1:3)]) 
    HapCovRow <- gsub(b, 1, n_A) %>% as.numeric()
  })
  saveRDS(dConvMatrix, file = paste0(output, save_name, "_MT_Genotypes.RDS"))
  return(dConvMatrix)
}

mEA <- convertHaplotype01(dSnp, dEA, "EA")
mAF <- convertHaplotype01(dSnp, dAF, "AF")
mWE <- convertHaplotype01(dSnp, dWE, "WE")
mSA <- convertHaplotype01(dSnp, dSA, "SA")
mO <- convertHaplotype01(dSnp, dO, "O")
mCAS <- convertHaplotype01(dSnp, dCAS, "CAS")
mAM <- convertHaplotype01(dSnp, dAM, "AM")
## takes a bit of time....
## SNP end up in the columns and individuals in the rows. But that is ok


