# Install and load necessary libraries
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# Define the file paths for the FASTQ files
fastq_files <- c(
  "../Dataset/PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq",
  "../Dataset/PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq",
  "../Dataset/PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq",
  "../Dataset/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq"
)

# Set output directory for FastQC results
output_dir <- "../Dataset/results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Run FastQC on each FASTQ file
for (file in fastq_files) {
  cmd <- paste("fastqc", file, "-o", output_dir)
  system(cmd)
}


