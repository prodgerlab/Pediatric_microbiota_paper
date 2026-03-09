## Pathological Phimosis is Associated with Foreskin Immune Cell Infiltration but Not Microbiota Composition
## Microbial diversity analysis - alpha diversiy
## Prepared by: Rachel Penney
## Reviewed by: Jorge Rojas-Vargas


# Load libraries
library(microbiome)
library(phyloseq)
library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)
library(vegan)


# Load data
counts<-read.table("cutadapt_counts.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
rownames(counts)[rownames(counts)=="NE-077_pre"] <- "NE_077_pre"
tax<-read.table("SILVA_v138.2_tax.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
adults <- read_excel("Shannon_div_peds_vs_adults.xlsx")

# Metadata
mtdt <- read.csv("metadata.csv", header = T, row.names = 1)



## -----------------------------
##
## ALPHA DIVERSITY
##
## -----------------------------

#Now combine samples from all the cohorts- but only used paired samples
dm <- counts[c("EL_003_pre", "EL_007_un", "EL_025_un1", "EL_036_pre", "EL_050_pre", "EL_051_pre", "EL_052_pre", "EL_059_pre", "EL_070_pre", "EL_074_pre", "EL_075_pre", "EL_078_pre", "NE_002_pre2", "NE_009_un", "NE_010_un1", "NE_018_pre", "NE_034_pre", "NE_035_pre", "NE_041_pre", "NE_056_pre", "NE_058_pre", "NE_063_pre", "NE_067_pre", "NE_069_pre", "NE_079_pre", "PA_001_pre", "PA_008_un", "PA_011_un1", "PA_013_pre", "PA_014_pre", "PA_015_pre", "PA_016_un1", "PA_029_pre2", "PA_030_pre2", "PA_033_pre", "PA_042_pre", "PA_044_pre", "PA_045_pre", "PA_046_pre", "PA_048_pre", "PA_049_pre", "PA_053_pre", "PA_054_pre", "PA_055_pre", "PA_057_pre", "PA_060_pre", "PA_061_pre2", "PA_064_pre", "PA_065_pre", "PA_076_pre", "PA_080_pre", "EL_003_pw", "EL_007_un1", "EL_025_un", "EL_036_pw", "EL_050_pw", "EL_051_pw", "EL_052_pw", "EL_059_pw", "EL_070_pw", "EL_074_pw", "EL_075_pw", "EL_078_pw", "NE_002_pw", "NE_009_un1", "NE_010_un", "NE_018_pw", "NE_034_pw", "EL_035_pw", "NE_041_pw", "NE_056_pw", "NE_058_pw", "NE_063_pw", "NE_067_pw", "NE_069_pw", "NE_079_pw", "PA_001_pw", "PA_008_un1", "PA_011_un", "PA_013_pw", "PA_014_pw", "PA_015_pw", "PA_016_un", "PA_029_pw", "PA_030_pw", "PA_033_pw", "PA_042_pw", "PA_044_pw", "PA_045_pw", "PA_046_pw", "PA_048_pw2", "PA_049_pw", "PA_053_pw", "PA_054_pw", "PA_055_pw", "PA_057_pw", "PA_060_pw", "PA_061_pw", "PA_064_pw", "PA_065_pw", "PA_076_pw", "PA_080_pw"),]

meta <- mtdt[c("EL_003_pre", "EL_007_un", "EL_025_un1", "EL_036_pre", "EL_050_pre", "EL_051_pre", "EL_052_pre", "EL_059_pre", "EL_070_pre", "EL_074_pre", "EL_075_pre", "EL_078_pre", "NE_002_pre2", "NE_009_un", "NE_010_un1", "NE_018_pre", "NE_034_pre", "NE_035_pre", "NE_041_pre", "NE_056_pre", "NE_058_pre", "NE_063_pre", "NE_067_pre", "NE_069_pre", "NE_079_pre", "PA_001_pre", "PA_008_un", "PA_011_un1", "PA_013_pre", "PA_014_pre", "PA_015_pre", "PA_016_un1", "PA_029_pre2", "PA_030_pre2", "PA_033_pre", "PA_042_pre", "PA_044_pre", "PA_045_pre", "PA_046_pre", "PA_048_pre", "PA_049_pre", "PA_053_pre", "PA_054_pre", "PA_055_pre", "PA_057_pre", "PA_060_pre", "PA_061_pre2", "PA_064_pre", "PA_065_pre", "PA_076_pre", "PA_080_pre", "EL_003_pw", "EL_007_un1", "EL_025_un", "EL_036_pw", "EL_050_pw", "EL_051_pw", "EL_052_pw", "EL_059_pw", "EL_070_pw", "EL_074_pw", "EL_075_pw", "EL_078_pw", "NE_002_pw", "NE_009_un1", "NE_010_un", "NE_018_pw", "NE_034_pw", "EL_035_pw", "NE_041_pw", "NE_056_pw", "NE_058_pw", "NE_063_pw", "NE_067_pw", "NE_069_pw", "NE_079_pw", "PA_001_pw", "PA_008_un1", "PA_011_un", "PA_013_pw", "PA_014_pw", "PA_015_pw", "PA_016_un", "PA_029_pw", "PA_030_pw", "PA_033_pw", "PA_042_pw", "PA_044_pw", "PA_045_pw", "PA_046_pw", "PA_048_pw2", "PA_049_pw", "PA_053_pw", "PA_054_pw", "PA_055_pw", "PA_057_pw", "PA_060_pw", "PA_061_pw", "PA_064_pw", "PA_065_pw", "PA_076_pw", "PA_080_pw"),]

# Calculate alpha diversity for each individual sample
shannon <- diversity(t(dm), index = "shannon")   # per sample

# Build dataframe
merged_df <- data.frame(
  Shannon = shannon,
  Time = meta$Time,
  Cohort = meta$Cohort
)

colnames(merged_df) <- c("Shannon", "Time", "Cohort")

merged_df$Time <- factor(merged_df$Time, levels = c("pre", "post"))

# Make a new color palette for 6 groups
mycolours <- c("red", "blue")

# Plot the boxplot of Shannon
Shannon <- boxplot(
  Shannon ~ Time, 
  data = merged_df, 
  xlab = "Time", 
  ylab = "Shannon's Index", 
  main = NULL, 
  col = mycolours, 
  las = 2, # Rotate x-axis labels for better visibility
  outline = FALSE # Hide outliers to avoid duplication with data points
)

# Add data points using stripchart
stripchart(
  Shannon ~ Time, 
  data = merged_df, 
  method = "jitter",  # Jitter the points to avoid overlap
  pch = 16,  # Use solid dots
  col = "black", # Color of the points
  vertical = TRUE,  # Align points with boxplot
  add = TRUE # Overlay on the existing boxplot
)



# Performed Shapiro wilk and p=0.09 so use parametric t-test
t_test_result <- t.test(Shannon ~ Time, data = merged_df, var.equal = TRUE)
print(t_test_result)

#p=0.0001363, Pre group mean: 3.367, Post group mean: 2.939

#Medians for the paper
tapply(merged_df$Shannon, merged_df$Time, median, na.rm = TRUE)[c("pre", "post")]
# pre     post 
#3.373256 3.150362

#Calculate the median values of pre vs post for the abstract of the paper

# Calculate median for each cohort
div_medians <- merged_df %>%
  group_by(Time) %>%
  summarise(Median_Shannon = median(Shannon, na.rm = TRUE))

print(div_medians)

#  Group Median_Shannon
#<chr>          <dbl>
#1 post            3.15
#2 pre             3.37


#Now to make the boxplot

# Ensure 'Cohort' is a factor with levels ordered as 'P' and then 'A'
adults$Cohort <- factor(adults$Cohort, levels = c("P", "A"))

mycolours <- c("orange", "limegreen")

# Plot the boxplot grouped by the reordered Group
Shannon <- boxplot(
  Shannon ~ Cohort, 
  data = adults, 
  xlab = "Cohort", 
  ylab = "Shannon's Index", 
  main = NULL, 
  col = mycolours, 
  las = 2, # Rotate x-axis labels for better visibility
  outline = FALSE # Hide outliers to avoid duplication with data points
)

# Add data points using stripchart
stripchart(
  Shannon ~ Cohort, 
  data = adults, 
  method = "jitter",  # Jitter the points to avoid overlap
  pch = 16,  # Use solid dots
  col = "black", # Color of the points
  vertical = TRUE,  # Align points with boxplot
  add = TRUE # Overlay on the existing boxplot
)


# Perform the Wilcoxon rank-sum test (unpaired) for the current cohort
test_results <- t.test(Shannon ~ Cohort, data = adults, var.equal = TRUE)

# Print the results
test_results
#p-value = <2.2e-16 Group P mean: 3.33, Group A meanL 2.61

#Median Shannon diversity for peds vs adults
# Base R approach
tapply(adults$Shannon, adults$Cohort, median, na.rm = TRUE)[c("P", "A")]
#       P        A 
#3.363856 2.594174
