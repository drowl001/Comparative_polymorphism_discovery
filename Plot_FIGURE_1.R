
library(vcfR)
library(ggplot2)
library(cowplot) 

# Load and process data for the first plot
vcf1 <- read.vcfR("rirreg.chen.fb.RiPE_filtered_cds_polymorphicONLY_biallelic.only_finalfiltrd.vcf", verbose = FALSE)
AD_rirreg1_1 <- extract.gt(vcf1, element = "AD", as.numeric = FALSE, IDtoRowNames = FALSE, extract = TRUE)
ad1_1 <- as.numeric(masplit(AD_rirreg1_1, delim = ",", sort = 0, record = 1))
ad2_1 <- as.numeric(masplit(AD_rirreg1_1, delim = ",", sort = 0, record = 2))
Afreq_1 <- ad2_1 / (ad1_1 + ad2_1)
Afreq_1 <- Afreq_1[!is.na(Afreq_1)]
Afreq_1 <- round(Afreq_1, digits = 3)
all_freqs_df_1 <- data.frame(freq = c(Afreq_1))

# Load and process data for the second plot
vcf2 <- read.vcfR("rirreg.w3.fb.p10_filtered_cds_polymorphicONLY_biallelic.only_finalfiltrd.vcf", verbose = FALSE)
AD_rirreg1_2 <- extract.gt(vcf2, element = "AD", as.numeric = FALSE, IDtoRowNames = FALSE, extract = TRUE)
parsed_AD <- matrix(NA, nrow = nrow(AD_rirreg1_2), ncol = ncol(AD_rirreg1_2))
AAF_per_sample <- matrix(NA, nrow = nrow(AD_rirreg1_2), ncol = ncol(AD_rirreg1_2))
for (i in 1:nrow(AD_rirreg1_2)) {
  for (j in 1:ncol(AD_rirreg1_2)) {
    if (!is.na(AD_rirreg1_2[i, j])) {
      ad_values <- strsplit(AD_rirreg1_2[i, j], ",")[[1]]
      ref_count <- as.numeric(ad_values[1])
      alt_count <- as.numeric(ad_values[2])
      dp <- ref_count + alt_count
      parsed_AD[i, j] <- alt_count
      AAF_per_sample[i, j] <- round(alt_count / dp, 3)
    }
  }
}
average_AAF_2 <- rowMeans(AAF_per_sample, na.rm = TRUE)
average_AAF_2 <- round(average_AAF_2, 3)
all_freqs_df_2 <- data.frame(freq = c(average_AAF_2))

# Set the DPI for the combined  output
pdf(file = "combined_plots.pdf", width = 16, height = 6, units = "in", res = 350)

# Create the first plot
plot1 <- ggplot(all_freqs_df_1, aes(x = freq)) + 
  geom_histogram(binwidth = 0.05, fill = "#4B0076", color = "white") +
  labs(x = "Alternate allele frequency", y = "Number of SNPs") +
  scale_y_continuous(limits = c(0, 125)) +
  theme(axis.title.x = element_text(colour = "black", size = 12, face = "bold"),  
        axis.title.y = element_text(colour = "black", size = 12, face = "bold"), 
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.2, colour = "black"),
        text = element_text(size = 12, face = "bold"))  

# Create the second plot
plot2 <- ggplot(all_freqs_df_2, aes(x = freq)) + 
  geom_histogram(binwidth = 0.05, fill = "darkgoldenrod", color = "white") +
  labs(x = "Alternate allele frequency", y = "Number of SNPs") +
  scale_y_continuous(limits = c(0, 125)) +
  theme(axis.title.x = element_text(colour = "black", size = 12, face = "bold"),  
        axis.title.y = element_text(colour = "black", size = 12, face = "bold"), 
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.2, colour = "black"),
        text = element_text(size = 12, face = "bold"))  

# Arrange and label the plots
combined_plot <- plot_grid(plot1 + labs(title = "A", title.bold = TRUE), plot2 + labs(title = "B", title.bold = TRUE), labels = NULL, ncol = 2)  

# Print the combined plot
print(combined_plot)

# Save the combined plot to a PNG file
ggsave("combined_plot_FIGURE1.pdf", combined_plot)

# Close the device
dev.off()
