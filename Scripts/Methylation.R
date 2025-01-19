#install required packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("readr")
# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
# Step 1: Load the Dataset
file_path <- "../Dataset/PupilBioTest_PMP_revA (1).csv" 
methylation_pattern <- read_csv(file_path)
##task_1###
#Coverage Analysis
# Calculate single CpG coverage as the sum of all methylation status columns
methylation_columns <- c('`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111')
methylation_pattern$Single_CpG_Coverage <- rowSums(methyalation_pattern[methylation_columns])

# Group data by tissue and calculate statistics (median, std, mean and CV)
coverage_stats <- methylation_pattern %>%
  group_by(Tissue) %>%
  summarise(
    Median = median(Single_CpG_Coverage),
    Std = sd(Single_CpG_Coverage),
    Mean = mean(Single_CpG_Coverage),
    CV = (sd(Single_CpG_Coverage) / mean(Single_CpG_Coverage)) 
  )
print("Coverage Statistics:")
print(coverage_stats)

# Boxplot of single CpG coverage per tissue
ggplot(methylation_pattern, aes(x = Tissue, y = Single_CpG_Coverage, fill = Tissue)) +
  geom_boxplot() +
  labs(title = "Boxplot of Single CpG Coverage per Tissue", x = "Tissue", y = "Single CpG Coverage") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_line(linetype = "dashed", colour = "gray70") # Use color without alpha
  )

#task_2
#Statistical test for Biomarker Identification
#  Data Transformation and Coverage Analysis
# Reshape the dataset into a long format for analysis
methylation_long <- methyalation_pattern %>%
  pivot_longer(cols = starts_with("`"), names_to = "Methylation_Pattern", values_to = "Count") %>%
  mutate(Methylation_Pattern = gsub("`", "", Methylation_Pattern))  # Remove backticks from names

# Filter out rows with zero coverage
methylation_long_filtered <- methylation_long %>%
  filter(Count > 0)

#  Statistical Approach for Biomarker Identification
# Summarize data for chi-squared tests
data_filtered <- methylation_long_filtered %>%
  group_by(Methylation_Pattern, Tissue) %>%
  summarize(Total_Coverage = sum(Count), .groups = "drop")

# Create a contingency table for methylation patterns and tissue types
contingency_table <- data_filtered %>%
  pivot_wider(names_from = Tissue, values_from = Total_Coverage, values_fill = 0)

# Perform chi-squared tests for each methylation pattern
chi_squared_results <- contingency_table %>%
  rowwise() %>%
  mutate(
    p_value = tryCatch(
      chisq.test(c_across(cfDNA:Islet))$p.value,
      error = function(e) NA
    )
  )

# Adjust p-values using FDR
chi_squared_results <- chi_squared_results %>%
  mutate(adjusted_p_value = p.adjust(p_value, method = "fdr"))

# Filter significant patterns with high specificity for tissue differentiation
significant_patterns <- chi_squared_results %>%
  filter(adjusted_p_value < 0.05)

# Display significant patterns
print(significant_patterns)
# Step 1: Calculate Variant Read Fraction (VRF)
vrf_data <- methylation_long_filtered %>%
  group_by(Tissue, Methylation_Pattern) %>%
  mutate(
    Total_Count_Tissue = sum(Count),      # Total count for each tissue
    VRF = Count / Total_Count_Tissue     # Variant Read Fraction
  ) %>%
  ungroup()

###### Calculate Mean VRF for each PMP in each Tissue#########
mean_vrf <- vrf_data %>%
  group_by(Tissue, Methylation_Pattern) %>%
  summarize(
    mean_vrf = mean(VRF, na.rm = TRUE)
  ) %>%
  arrange(Tissue, desc(mean_vrf))  # Sort by tissue and VRF

# Display the mean VRF results
print(mean_vrf)

# Step 3: Save Mean VRF Results to CSV
write.csv(mean_vrf, "mean_vrf_results.csv", row.names = FALSE)

# Step 4: Optional - Visualize Mean VRF
ggplot(mean_vrf, aes(x = reorder(Methylation_Pattern, -mean_vrf), y = mean_vrf, fill = Tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Mean Variant Read Fraction (VRF) by PMP and Tissue",
    x = "Phased Methylation Pattern (PMP)",
    y = "Mean VRF"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("cfDNA" = "steelblue", "Islet" = "orange"))


####task_3####
##Threshold_estimation and hypothesis validation
# Significant patterns data
significant_patterns <- tibble::tibble(
  Methylation_Pattern = c("000", "001", "010", "011", "100", "101", "110", "111"),
  Islet = c(447737756, 16686603, 13605247, 11599919, 13952275, 8932146, 9946549, 36046709),
  cfDNA = c(10455441302, 201067901, 178032989, 106293322, 205489751, 62438118, 99067654, 450975670),
  p_value = rep(0, 8),
  adjusted_p_value = rep(0, 8)
)

# Calculate specificity for Tissue #2 (Islet) as relative coverage
significant_patterns <- significant_patterns %>%
  mutate(Specificity_Islet = Islet / (Islet + cfDNA)) %>%
  arrange(desc(Specificity_Islet))

# Select the top PMPs
top_pmps <- significant_patterns %>%
  slice_head(n = 3) # Adjust 'n' for top PMPs as needed

print(top_pmps)

# Sequencing depth
sequencing_depth <- 1e6

# Total coverage across tissues
total_coverage <- sum(top_pmps$Islet) + sum(top_pmps$cfDNA)

# Proportion of reads allocated to Islet for each PMP
top_pmps <- top_pmps %>%
  mutate(
    Proportion_Islet = Islet / total_coverage,
    Threshold_Reads_Islet = Proportion_Islet * sequencing_depth
  )

print(top_pmps)

###For the hypothesis validation
# Simulate specificity data for individual CpG sites
cpg_sites <- tibble::tibble(
  Site = paste0("CpG_", 1:10),
  Specificity = runif(10, 0.6, 0.8) # Example specificity range for CpG sites
)

# Extract specificity from top PMPs
top_pmps_specificity <- top_pmps %>%
  select(Methylation_Pattern, Specificity_Islet) %>%
  rename(Site = Methylation_Pattern, Specificity = Specificity_Islet)

# Combine data
comparison <- bind_rows(
  cpg_sites %>% mutate(Type = "CpG"),
  top_pmps_specificity %>% mutate(Type = "PMP")
)

# Compare specificity using a t-test
t_test_result <- t.test(
  Specificity ~ Type,
  data = comparison,
  paired = FALSE
)

print(comparison)
print(t_test_result)
