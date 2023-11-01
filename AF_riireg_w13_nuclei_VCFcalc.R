
library(vcfR)
library(ggplot2)

vcf <- read.vcfR("rirreg.w3.fb.p10_filtered_cds_polymorphicONLY_biallelic.only.vcf", verbose = FALSE)
localchromR <- create.chromR(vcf)
AD_rirreg1 <-extract.gt(vcf, element="AD", as.numeric = F, IDtoRowNames = F, extract = TRUE)
C_name<- getCHROM(vcf)
C_POS<- getPOS(vcf)

ad1 <- masplit(AD_rirreg1, delim = ",",  sort = 0, record = 1)
ad2 <- masplit(AD_rirreg1, delim = ",",  sort = 0,record = 2)

#to mask SNPs by fraction
NR <- ad1 + ad2

NR[NR < 5] <- 0 # 5 reads
(out1<- which(NR == 0))
ad1<-replace(ad1, out1, 0)
ad2<-replace(ad2, out1, 0)
Afreq <- ad2/(ad1+ad2)
Af<-round(Afreq, digits = 3)
Af[is.na(Af)] <-"."
Af[Af > 0.1 & Af < 0.9] <- "."
Af[Af >=  0.9] <- 1
Af[Af <= 0.1 & Af != "."] <- 0
File <- cbind(C_name, C_POS, Af)

Af <- File[rowSums(Af == ".") <= 4, ] #33 pc max missing data

#calculate Allele frequencies
Alt <- rowSums(Af == 1)
Ref<- rowSums(Af == 0)
AN <- Alt + Ref
AAF <- Alt/AN # to calculate Alternate allele frequency by dividing the number of ones by the number of rows/sample (24 in this case)
#RAF<-Ref/AN
Alt_allele_Frequency<-round(AAF, digits = 3)
print(Alt_allele_Frequency)


alt_freq_df <- data.frame(freq = c(Alt_allele_Frequency))

subset_file <- Final_file[,1:2]
nrow(subset_file)

write.table(subset_file, file = "rirreg.w3.fb.p10_filtered_cds_polymorphicONLY_biallelic.only.vcf_subset_file.txt", sep = "\t", quote=FALSE, col.names = FALSE, row.names = F)#this is the file I use to subset the original vcf file to make the new filtered vcf
