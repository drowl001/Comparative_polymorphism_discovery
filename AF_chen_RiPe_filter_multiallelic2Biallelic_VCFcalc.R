
library(vcfR)
library(ggplot2)

vcf <- read.vcfR("rirreg.chen.fb.RiPE_filtered_polymorphic.vcf", verbose = FALSE)
localchromR <- create.chromR(vcf)
AD_rirreg1 <-extract.gt(vcf, element="AD", as.numeric = F, IDtoRowNames = F, extract = TRUE)
C_name<- getCHROM(vcf)
C_POS<- getPOS(vcf)
FA <- extract.info(vcf, element="AF", as.numeric = T, mask = FALSE)

#Check for multiple alleles and subset only biallelic sites
# Function to check the number of alleles in AD field
get_allele_count <- function(x) lengths(strsplit(x, ","))

# Check the number of alleles at each site
allele_counts <- sapply(AD_rirreg1, get_allele_count)

# Count the occurrences of sites with 2, 3, 4, and 5 alleles
count_alleles <- table(allele_counts)

# Print the counts
print(count_alleles)

# Check for biallelic sites (sites with exactly 2 values in the AD field)
biallelic_sites <- allele_counts == 2

# Filter  biallelic sites
AD_rirreg2 <- AD_rirreg1[biallelic_sites]
C_name <- getCHROM(vcf)[biallelic_sites]
C_POS <- getPOS(vcf)[biallelic_sites]

# Create a new data frame with only biallelic sites
File <- data.frame(C_name = C_name, C_POS = C_POS, AD_rirreg2 = AD_rirreg2)

# Print the updated data frame containing only biallelic sites
print(File)

# Split the AD_rirreg2 column into two separate columns
AD_values <- strsplit(as.character(File$AD_rirreg2), ",")

# Convert the split values to numeric
AD_value1 <- as.numeric(sapply(AD_values, function(x) x[1]))
AD_value2 <- as.numeric(sapply(AD_values, function(x) x[2]))

# Calculate the sum of AD_value1 and AD_value2
NR <- AD_value1 + AD_value2

# Print NR to check the values
print(NR)

count_less_than_5 <- sum(NR < 5)
# Print the result
print(count_less_than_5)

NR[NR < 5] <- 0 #5 reads
(out1<- which(NR == 0))
AD_value1<-replace(AD_value1, out1, 0)
AD_value2<-replace(AD_value2, out1, 0)
freq <- AD_value2/(AD_value1 + AD_value2)
Af<-round(freq, digits = 3)
Af[is.na(Af)] <-NA
Final_File <- data.frame(C_name = File$C_name, C_POS = File$C_POS, Af = Af)


# Exclude rows with missing data from the Final_File dataframe
Final_File <- Final_File[complete.cases(Final_File), ]

# Get the number of rows in the final dataframe
nrow(Final_File)

subset_file <- Final_File[,1:2]

write.table(subset_file, file = "rirreg.chen.fb.RiPE_filtered_polymorphic.vcf_final_subset_file.txt", sep = "\t", quote=FALSE, col.names = FALSE, row.names = F)#this is the file I use to subset the original vcf file to make the new filtered vcf
