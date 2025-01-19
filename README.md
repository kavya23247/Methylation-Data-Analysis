# README: Methylation Pattern Analysis Project
# Part 1: Methylation Analysis
## Overview
This project analyzes methylation patterns in biological datasets to:
1. Assess coverage of CpG sites across various tissues.
2. Identify significant Phased Methylation Patterns (PMPs) differentiating tissues.
3. Estimate thresholds and validate biomarkers for methylation-based studies.

## Features
- **Coverage Analysis**: Calculate coverage of CpG sites and summarize statistics (mean, median, standard deviation, coefficient of variation).
- **Pattern Analysis**: Identify tissue-specific methylation patterns and validate using chi-squared tests.
- **Threshold Estimation**: Quantify specificity and threshold reads for biomarker validation.
- **Visualization**: Generate insightful plots using `ggplot2`.

## Prerequisites
Ensure the following R packages are installed:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "readr"))
```

# Part 2: NGS Analysis
# Overview: Running FastQC on FASTQ Files

This script performs quality control checks on a set of FASTQ files using FastQC, a tool designed to assess the quality of sequencing data.

## Key Steps

1. **Install and Load Required Libraries**:
   The script ensures that the `ggplot2` library is installed and loaded. This library is not directly used in the code but might be required for downstream visualizations.

   ```r
   if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
   library(ggplot2)
   ```

