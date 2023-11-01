library(ggplot2)

Single_nuclei_dataset <- c(0.961, 0.684, 0.68, 0.68, 0.675, 0.657, 0.605, 0.604, 0.603, 0.602, 0.547, 0.528, 0.518, 0.518, 0.514, 0.509, 0.501, 0.477, 0.473, 0.461, 0.436, 0.434, 0.153, 0.149, 0.104)
Whole_organism_dataset <- c(0.502, 0.514, 0.51, 0.497, 0.495, 0.487, 0.479, 0.471, 0.463, 0.458, 0.395, 0.391, 0.38, 0.377, 0.344, 0.34, 0.338, 0.305, 0.285, 0.282, 0.241, 0.235, 0.228, 0.226, 0.136)

# Create a data frame 
df <- data.frame(
  Single_nuclei_dataset = Single_nuclei_dataset,
  Whole_organism_dataset = Whole_organism_dataset
)



# Create a scatter plot
plot<- ggplot(df, aes(x = Single_nuclei_dataset, y = Whole_organism_dataset)) +
  geom_point()+
  labs(x = "13n single nuclei dataset", y = "Whole organism dataset") 

print(plot)



ggsave("FIGURE_S1_correlation.png", plot,  width = 6, height = 6, units = "in", dpi = 300, device='png')


# Perform a Pearson correlation test
correlation_result <- cor.test(Single_nuclei_dataset, Whole_organism_dataset, method = "pearson")


correlation_result
