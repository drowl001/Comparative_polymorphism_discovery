
library(ggplot2)
library(cowplot)

# List of dataset filenames
dataset_filenames <- c(
  "261_sites_bcftulsmpileup_FINAL.txt",
  "261_random_RiPE.bam_sites1_bcftulsmpileup.txt",
  "261_random_RiPE.bam_sites2_bcftulsmpileup.txt",
  "261_random_RiPE.bam_sites3_bcftulsmpileup.txt",
  "261_random_RiPE.bam_sites4_bcftulsmpileup.txt",
  "261_random_RiPE.bam_sites5_bcftulsmpileup.txt",
  "261_random_RiPE.bam_sites6_bcftulsmpileup.txt",
  "261_random_RiPE.bam_sites7_bcftulsmpileup.txt",
  "261_random_RiPE.bam_sites8_bcftulsmpileup.txt",
  "261_random_RiPE.bam_sites9_bcftulsmpileup.txt",
  "261_random_RiPE.bam_sites10_bcftulsmpileup.txt"
)

# List of dataset labels
dataset_labels <- c(
  "Non-detected 13n dataset",
  "Random dataset 1",
  "Random dataset 2",
  "Random dataset 3",
  "Random dataset 4",
  "Random dataset 5",
  "Random dataset 6",
  "Random dataset 7",
  "Random dataset 8",
  "Random dataset 9",
  "Random dataset 10"
)

# Create a data frame to store the counts
counts_data <- data.frame(Dataset = character(0), 
                          "Bi-allelic sites" = numeric(0), 
                          "Multi-allelic sites" = numeric(0), 
                          "Invariant sites" = numeric(0))

# Create a data frame to store the data for the plot
plot_data <- data.frame()

# Defining the colors
colors <- c("#4B0076" ,"grey","#ECE7F2")

for (i in 1:length(dataset_filenames)) {
  # Read the text files
  data <- read.table(dataset_filenames[i], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  # Function to remove the last digit from column 6
  remove_last_digit <- function(s) {
    parts <- strsplit(s, ",")
    sapply(parts, function(x) {
      if (length(x) > 2) {
        x[length(x)] <- sub("\\d$", "", x[length(x)])
      }
      paste(x, collapse = ",")
    })
  }
  
  # Apply the function to column 6
  data$V6 <- remove_last_digit(data$V6)
  
  # Print to verify the modified data
  print(data)
  
  # Create a new column 'allele_count' to store the number of alleles
  data$allele_count <- sapply(data$V6, function(allele_string) {
    alleles <- unlist(strsplit(allele_string, ","))
    return(length(alleles))
  })
  
  # Create subsets based on allele counts
  quintuplelic_snps <- data[data$allele_count == 5, ]
  nrow(quintuplelic_snps)
  
  quadruplelic_snps <- data[data$allele_count == 4, ]
  nrow(quadruplelic_snps)
  
  triallelic_snps <- data[data$allele_count == 3, ]
  nrow(triallelic_snps)
  
  biallelic_snps <- data[data$allele_count == 2, ]
  nrow(biallelic_snps)
  
  
  # Function to calculate AAF
  calculate_AAF <- function(row) {
    alleles <- unlist(strsplit(row, ","))
    if (length(alleles) >= 2) {
      digit1 <- as.numeric(alleles[1])
      digit2 <- as.numeric(alleles[2])
      if (is.na(digit1) || is.na(digit2)) {
        return(NA)
      }
      aaf <- digit2 / (digit1 + digit2)
      return(aaf)
    } else {
      return(NA)
    }
  }
  
  # Calculate AAF for biallelic_snps
  biallelic_snps$AAF <- sapply(biallelic_snps$V6, calculate_AAF)
  
  # Add a column for the status
  biallelic_snps$Status <- ifelse(biallelic_snps$AAF == 0, "Invariant sites", "Bi-allelic sites")
  
  #print out the numbers for the different categories in the dataset
  Biallelic <- biallelic_snps[biallelic_snps$AAF  > 0.000, ]
  nrow(Biallelic)
  
  Invariant <- biallelic_snps[biallelic_snps$AAF == 0.000, ]
  nrow(Invariant)
  
  
  # Create a data frame to store counts of each type
  site_counts <- data.frame(
    Site_Type = factor(c("Bi-allelic sites", "Multi-allelic sites", "Invariant sites"), 
                       levels = c("Invariant sites", "Bi-allelic sites", "Multi-allelic sites")),
    Count = c(
      nrow(biallelic_snps[biallelic_snps$AAF != 0, ]), 
      nrow(triallelic_snps) + nrow(quadruplelic_snps) + nrow(quintuplelic_snps),
      nrow(biallelic_snps[biallelic_snps$AAF == 0, ])
    )
  )
  
  #Get multi-allelic sites
  MultiAllelic <- nrow(data[data$allele_count > 2, ])
  
  # Calculate counts for different categories within each dataset
  Biallelic <- nrow(biallelic_snps[biallelic_snps$AAF > 0.000, ])
  Invariant <- nrow(biallelic_snps[biallelic_snps$AAF == 0.000, ])
  
  # Create a data frame for counts
  counts_entry <- data.frame(
    Dataset = dataset_labels[i],
    "Bi-allelic sites" = Biallelic,
    "Multi-allelic sites" = MultiAllelic,
    "Invariant sites" = Invariant
  )
  
  # Append to the counts_data data frame
  counts_data <- rbind(counts_data, counts_entry)
  
  # Create a data frame for the plot
  plot_entry <- data.frame(
    Dataset = dataset_labels[i],
    Category = c("Bi-allelic sites", "Multi-allelic sites", "Invariant sites"),
    Count = c(Biallelic, MultiAllelic, Invariant)
  )
  
  # Append to the plot_data data frame
  plot_data <- rbind(plot_data, plot_entry)
}

# Make order of the x-axis labels
plot_data$Dataset <- factor(plot_data$Dataset, levels = dataset_labels)

# Print the counts
print(counts_data)

hist <- ggplot(plot_data, aes(x = Dataset, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "Counts") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, face = "bold", size = 12),  
        axis.text.y = element_text(face = "bold", size = 12),  
        legend.position = "right",
        axis.title = element_text(face = "bold", size = 14, color = "black")) +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.2, colour = "black")) +
  guides(fill = guide_legend(title = NULL)) +
  theme(legend.text = element_text(face = "bold")) 


hist




# Create an empty data frame to store the counts of different site categories
all_counts_data <- data.frame(
  File = character(),
  Category = character(),
  Count = numeric()
)



# Create a list of file names
file_names <- c(
  "470_sites_unique_inCHEN_bcftulsmpileup_FINAL.txt",
  "random_sites_from13_bams_473_bcftulsmpileup.txt1",
  "random_sites_from13_bams_473_bcftulsmpileup.txt2",
  "random_sites_from13_bams_473_bcftulsmpileup.txt3",
  "random_sites_from13_bams_473_bcftulsmpileup.txt4",
  "random_sites_from13_bams_473_bcftulsmpileup.txt5",
  "random_sites_from13_bams_473_bcftulsmpileup.txt6",
  "random_sites_from13_bams_473_bcftulsmpileup.txt7",
  "random_sites_from13_bams_473_bcftulsmpileup.txt8",
  "random_sites_from13_bams_473_bcftulsmpileup.txt9",
  "random_sites_from13_bams_473_bcftulsmpileup.txt10"
)


# Create corresponding labels
file_labels <- c(
  "Non-detected whole organism dataset",
  "Random dataset 1",
  "Random dataset 2",
  "Random dataset 3",
  "Random dataset 4",
  "Random dataset 5",
  "Random dataset 6",
  "Random dataset 7",
  "Random dataset 8",
  "Random dataset 9",
  "Random dataset 10"
)

# Process each file
for (i in 1:length(file_names)) {
  file_name <- file_names[i]
  file_label <- file_labels[i]
  # Read the input file
  df <- read.delim(file_name, sep = "\t", header = FALSE)
  
  # Loop through columns from 6 to the last
  for (i in 6:ncol(df)) {
    # Check if the column contains comma-separated values
    if (all(grepl(",", df[[i]]))) {
      # Split the column values by commas
      split_values <- strsplit(df[[i]], ",")
      
      # Filter and modify values with more than 2 elements
      split_values <- lapply(split_values, function(x) {
        if (length(x) > 2) {
          x <- x[-length(x)]  # Remove the last element
        }
        return(x)
      })
      
      # Combine the modified values back into a single string
      df[[i]] <- sapply(split_values, paste, collapse = ",")
    }
  }
  
  df
  data<-df
  
  
  
  #####
  # Create a new column 'counts' to store the number of unique alleles
  data$counts <- 0
  
  # Loop through the columns and populate the counts column based on the number of alleles
  num_snps <- nrow(data)
  num_samples <- ncol(data) - 5  # Adjust the starting column index
  
  for (j in 1:num_snps) {
    for (i in 6:(num_samples + 5)) {  # Adjust the starting column index
      alleles <- data[j, i]
      if (!is.na(alleles)) {
        num_unique_alleles <- length(unlist(strsplit(alleles, ",")))
        data$counts[j] <- num_unique_alleles
        break
      }
    }
  }
  
  # Create subsets for different categories
  biallelic_snps <- data[data$counts == 2, ]
  nrow(biallelic_snps)
  triallelic_snps <- data[data$counts == 3, ]
  nrow(triallelic_snps)
  quadruplelic_snps <- data[data$counts == 4, ]
  nrow(quadruplelic_snps)
  quintuplelic_snps <- data[data$counts == 5, ]
  nrow(quintuplelic_snps)
  
  # Create a subset for the remaining rows
  remaining <- data[data$counts != 2 & data$counts != 3 & data$counts != 4 & data$counts != 5, ]
  
  # Print or manipulate the resulting subsets as needed
  print("Biallelic SNPs:")
  print(biallelic_snps)
  nrow(biallelic_snps)
  
  print("Triallelic SNPs:")
  print(triallelic_snps)
  nrow(triallelic_snps)
  
  print("Quadruplelic SNPs:")
  print(quadruplelic_snps)
  nrow(quadruplelic_snps)
  
  print("Quintuplelic SNPs:")
  print(quintuplelic_snps)
  nrow(quintuplelic_snps)
  
  # Remove the count column (column 6 to the end) and create a new data frame
  new_DF <- biallelic_snps[, -ncol(biallelic_snps)]

  
  # Create new columns for REF and ALT counts in new_DF
  for (col in 6:18) {
    col_name <- paste0("V", col)
    
    # Split the comma-separated values
    split_values <- strsplit(new_DF[, col_name], ",")
    
    # Extract REF and ALT values
    ref_values <- sapply(split_values, function(x) as.numeric(x[1]))
    alt_values <- sapply(split_values, function(x) as.numeric(x[2]))
    
    # Create new column names for REF and ALT
    ref_col_name <- paste0(col_name, "_REF")
    alt_col_name <- paste0(col_name, "_ALT")
    
    # Add REF and ALT columns to new_DF
    new_DF[, ref_col_name] <- ref_values
    new_DF[, alt_col_name] <- alt_values
  }
  
  # Remove the original columns with comma-separated values
  new_DF <- new_DF[, -c(3:18)]
  nrow(new_DF)
  
  # Define the prefix for the REF and ALT columns in new_DF
  column_prefix <- c("V6_", "V7_", "V8_", "V9_", "V10_", "V11_", "V12_", "V13_", "V14_", "V15_", "V16_", "V17_", "V18_")
  
  # Create ref_columns and alt_columns based on the column_prefix
  ref_columns <- paste0(column_prefix, "REF")
  alt_columns <- paste0(column_prefix, "ALT")
  
  
  # Extract the REF and ALT columns from new_DF based on ref_columns and alt_columns
  ref_data <- new_DF[, ref_columns]
  alt_data <- new_DF[, alt_columns]
  
  
  # Calculate the sum of REF and ALT counts for each column
  total_counts <- ref_data + alt_data
  
  # Replace values where the total count is less than 5 with 0
  total_counts[total_counts < 5] <- 0
  
  # Calculate the Alternate Allele Frequency (AAF) with three decimal places
  AAF <- alt_data / total_counts
  AAF <- round(AAF, digits = 3)
  
  # Replace NA values with '.'
  AAF[is.na(AAF)] <- '.'
  
  # Assign '.' to positions with AAF > 0.1 and AAF < 0.9
  AAF[(AAF > 0.1 & AAF < 0.9) | total_counts < 5] <- '.'
  
  # Assign '1' to positions with AAF >= 0.9
  AAF[AAF >= 0.9] <- '1'
  
  # Assign '0' to positions with AAF <= 0.1 (excluding '.')
  AAF[(AAF <= 0.1 & AAF != '.')] <- '0'
  
  
  File <- cbind(C_name = new_DF$V1, C_POS = new_DF$V2, AAF)
  str(File)
  File
  
  # Select the columns of interest (from the third column to the second last column)
  selected_columns <- File[, -(1:2)]
  
  # Create a new column ALT_nuclei
  File$ALT_nuclei <- rowSums(selected_columns == "1", na.rm = TRUE)
  File
  
  
  DT <- File[, -c(1, 2, ncol(File))]
  
  
  # Count the number of columns with '1' or '0' only (excluding '.')
  counts_ALT <- rowSums(DT == '1'| DT == '.')
  nrow(counts_ALT)
  
  # Count the number of columns with '1' or '0' only (excluding '.')
  counts <- rowSums(DT == '1' | DT == '0', na.rm = TRUE)
  
  # Add the 'counts' column to the 'counts_new' dataframe
  File$counts <- counts
  
  new_DF$counts <- File$counts
  
  new_DF$ALT_nuclei <- File$ALT_nuclei

  
  # Create a new column "status" based on ALT columns
  new_DF$status <- ifelse(rowSums(new_DF[, seq(4, ncol(new_DF), by = 2)] != 0) > 0, "Polymorphic sites", "Invariant sites")
  
  # View the updated data frame
  new_DF
  
  
  Pass_Mask_sites <- new_DF[new_DF$counts >= 9, ]
  nrow(Pass_Mask_sites)
  
  
  Fail_Mask_sites <-new_DF[new_DF$counts < 9, ]
  nrow(Fail_Mask_sites)
  
  REF_sites<-new_DF[new_DF$status == "Invariant sites", ]
  nrow(REF_sites)
  
  Polymorphic_sites<-new_DF[new_DF$status == "Polymorphic sites", ]
  nrow(Polymorphic_sites)
  
  Polymorphic_Pass<-Pass_Mask_sites[Pass_Mask_sites$status == "Polymorphic sites", ]
  nrow(Polymorphic_Pass)
  
  Alt_nuc<- new_DF[new_DF$ALT_nuclei >= 1, ]
  nrow(Alt_nuc)
  
  No_nuc<-Polymorphic_sites[Polymorphic_sites$counts == 0, ]
  nrow(No_nuc)
  
  With_nuc<-Polymorphic_sites[Polymorphic_sites$counts > 0, ]
  nrow(With_nuc)
  
  With_nuc_pass<-Polymorphic_sites[Polymorphic_sites$counts >= 9, ]
  nrow(With_nuc_pass)

  # Create a new data frame combining Triallelic and Quadruplelic SNPs as multiallelic
  multiallelic_snps <- rbind(triallelic_snps, quadruplelic_snps)
  
  
  # Create a data frame for the counts
  counts_data <- data.frame(
    Category = factor(c("Polymorphic sites", "Invariant sites", "Multi-allelic sites"), 
                      levels = c("Polymorphic sites", "Invariant sites", "Multi-allelic sites")),
    Count = c(nrow(Polymorphic_sites), nrow(REF_sites), nrow(multiallelic_snps))
  )
  
  # Add the file label to the counts_data
  counts_data$File <- factor(file_label, levels = file_labels)
  
  # Append counts_data to the all_counts_data
  all_counts_data <- rbind(all_counts_data, counts_data)
}

# Create a vertical stacked bar plot with custom x-axis labels at a 60-degree angle
stacked_bar_plot <- ggplot(data = all_counts_data, aes(x = File, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Counts") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, face = "bold", size = 12),  
        axis.text.y = element_text(face = "bold", size = 12),  
        axis.title = element_text(colour = "black", size = 12, face = "bold"),  
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.2, colour = "black")) +
  theme(legend.position = "right", legend.title = element_blank()) +
  theme(legend.text = element_text(face = "bold")) +
  scale_fill_manual(values = c("Polymorphic sites" = "darkgoldenrod", "Invariant sites" = "gray", "Multi-allelic sites" = "lightgoldenrod")) +
  guides(fill = guide_legend(title = "Category")) +
  labs(fill = "Category")


# Print the stacked bar plot
print(stacked_bar_plot)

# Add labels to the plots
hist <- hist + labs(title = "A", title.size = 14) + theme(plot.title = element_text(face = "bold"))
stacked_bar_plot <- stacked_bar_plot + labs(title = "B", title.size = 14) + theme(plot.title = element_text(face = "bold"))

# Adjust the margin for the hist plot to make sure x-axis label is not cut off
hist <- hist + theme(plot.margin = margin(b = 40))

# Create a plot without labels
combined_plot <- plot_grid(
  hist,
  stacked_bar_plot
)

print(combined_plot) 

# Save the combined plot to a PNG file
ggsave("combined_plot_FIGURE3.png", combined_plot, width = 16, height = 8)
