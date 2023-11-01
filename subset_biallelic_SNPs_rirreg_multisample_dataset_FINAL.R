
library(vcfR)

vcf <- read.vcfR("rirreg.w3.fb.p10.21nuc_filtered_cds_polymorphicONLY.vcf4", verbose = FALSE)
localchromR <- create.chromR(vcf)
AD_rirreg1 <-extract.gt(vcf, element="AD", as.numeric = F, IDtoRowNames = F, extract = TRUE)
C_name<- getCHROM(vcf)
C_POS<- getPOS(vcf)
File <- cbind(C_name, C_POS, AD_rirreg1)

File <- as.data.frame(File, stringsAsFactors = FALSE)

# Add a new column called "counts" at the end of the File data frame
File$counts <- NA

# Loop through the columns and populate the counts column based on the number of alleles
num_snps <- nrow(File)
num_samples <- ncol(File) - 2

for (j in 1:num_snps) {
  for (i in 3:(num_samples + 2)) {
    alleles <- File[j, i]
    if (!is.na(alleles)) {
      num_unique_alleles <- length(na.omit(strsplit(alleles, ",")[[1]]))
      File$counts[j] <- num_unique_alleles
      break
    }
  }
}

# Subset the dataframe to include only biallelic, triallelic, quadruplelic, and quintuplelic SNP sites
biallelic_snps <- File[File$counts == 2, ]
nrow(biallelic_snps)
triallelic_snps <- File[File$counts == 3, ]
nrow(triallelic_snps)
quadruplelic_snps <- File[File$counts == 4, ]
nrow(quadruplelic_snps)
quintuplelic_snps <- File[File$counts == 5, ]
nrow(quintuplelic_snps)

subset_file<-biallelic_snps[, 1:2]
nrow(subset_file)

write.table(subset_file, file = "rirreg.w3.fb.p10_filtered_cds_polymorphicONLY.vcf4_biallelic_sites_subset_file.txt",  sep = "\t", quote=FALSE, col.names = F, row.names = F)
