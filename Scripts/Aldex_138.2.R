## Pathological Phimosis is Associated with Foreskin Immune Cell Infiltration but Not Microbiota Composition
## Differential analysis
## Prepared by: Rachel Penney
## Reviewed by: Jorge Rojas-Vargas

# Load libraries
library(tidyr)
library(dplyr)
library(ALDEx2) 
library(ggrepel)
library(ggplot2)

# Load data
counts <- read.table("agg_counts_138.2.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)
rel <- read.table("rel_abund_138.2.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

subset <- counts[, c("EL_003_pre", "EL_007_un", "EL_025_un1", "EL_036_pre", "EL_050_pre", "EL_051_pre", "EL_052_pre", "EL_059_pre", "EL_070_pre", "EL_074_pre", "EL_075_pre", "EL_078_pre", "NE_002_pre2", "NE_009_un", "NE_010_un1", "NE_018_pre", "NE_034_pre", "NE_035_pre", "NE_041_pre", "NE_056_pre", "NE_058_pre", "NE_063_pre", "NE_067_pre", "NE_069_pre", "NE_079_pre", "PA_001_pre", "PA_008_un", "PA_011_un1", "PA_013_pre", "PA_014_pre", "PA_015_pre", "PA_016_un1", "PA_029_pre2", "PA_030_pre2", "PA_033_pre", "PA_042_pre", "PA_044_pre", "PA_045_pre", "PA_046_pre", "PA_048_pre", "PA_049_pre", "PA_053_pre", "PA_054_pre", "PA_055_pre", "PA_057_pre", "PA_060_pre", "PA_061_pre2", "PA_064_pre", "PA_065_pre", "PA_076_pre", "PA_080_pre", "EL_003_pw", "EL_007_un1", "EL_025_un", "EL_036_pw", "EL_050_pw", "EL_051_pw", "EL_052_pw", "EL_059_pw", "EL_070_pw", "EL_074_pw", "EL_075_pw", "EL_078_pw", "NE_002_pw", "NE_009_un1", "NE_010_un", "NE_018_pw", "NE_034_pw", "EL_035_pw", "NE_041_pw", "NE_056_pw", "NE_058_pw", "NE_063_pw", "NE_067_pw", "NE_069_pw", "NE_079_pw", "PA_001_pw", "PA_008_un1", "PA_011_un", "PA_013_pw", "PA_014_pw", "PA_015_pw", "PA_016_un", "PA_029_pw", "PA_030_pw", "PA_033_pw", "PA_042_pw", "PA_044_pw", "PA_045_pw", "PA_046_pw", "PA_048_pw2", "PA_049_pw", "PA_053_pw", "PA_054_pw", "PA_055_pw", "PA_057_pw", "PA_060_pw", "PA_061_pw", "PA_064_pw", "PA_065_pw", "PA_076_pw", "PA_080_pw")]

counts_t <-t(subset)


# Generate a relative abundance table and remove SVs accounting for < 1% of reads in every sample
# If you care about very rare things you could lower this to 0.1% (0.001)
props <- apply(counts_t, 1, function(x) {x/sum(x)})
filt_by_props <- counts_t[, apply(props, 1, max) >= 0.01]

tc <- t(filt_by_props)

colnames(tc)[colnames(tc) == "EL_007_un"] <- "EL_007_pre"
colnames(tc)[colnames(tc) == "EL_007_un1"] <- "EL_007_pw"
colnames(tc)[colnames(tc) == "EL_025_un1"] <- "EL_025_pre"
colnames(tc)[colnames(tc) == "EL_025_un"] <- "EL_025_pw"
colnames(tc)[colnames(tc) == "NE_009_un"] <- "NE_009_pre"
colnames(tc)[colnames(tc) == "NE_009_un1"] <- "NE_009_pw"
colnames(tc)[colnames(tc) == "NE_010_un1"] <- "NE_010_pre"
colnames(tc)[colnames(tc) == "NE_010_un"] <- "NE_010_pw"
colnames(tc)[colnames(tc) == "PA_008_un"] <- "PA_008_pre"
colnames(tc)[colnames(tc) == "PA_008_un1"] <- "PA_008_pw"
colnames(tc)[colnames(tc) == "PA_011_un1"] <- "PA_011_pre"
colnames(tc)[colnames(tc) == "PA_011_un"] <- "PA_011_pw"
colnames(tc)[colnames(tc) == "PA_016_un1"] <- "PA_016_pre"
colnames(tc)[colnames(tc) == "PA_016_un"] <- "PA_016_pw"

# specify your groups
pre<-grep("pre", colnames(tc), value = TRUE)
pw<-grep("pw", colnames(tc), value = TRUE)


aldex.in<-tc[,c(pre, pw)]

conds<-c(rep("pre", length(pre)), rep("pw", length(pw)))

# get the clr values
# this is the main ALDEx function for all downstream analyses
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)

# perform t-test (both Welches and Wilcoxon, plus a Benjamini-Hochberg multiple test correction)
x.tt <- aldex.ttest(x, paired.test=FALSE)

# estimate effect size and the within and between condition values
# include indiv. samples or not
x.effect <- aldex.effect(x)

# merge the data
x.all <- data.frame(x.tt, x.effect)


# Making volcano plot

# colours based on significance!
# add a column of NAs
x.all$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
x.all$diffexpressed[x.all$effect > 0.49] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
x.all$diffexpressed[x.all$effect < -0.48] <- "DOWN"

#labels based on significance!
x.all$delabel <- NA
x.all$delabel[x.all$diffexpressed != "NO"] <- rownames(x.all)[x.all$diffexpressed != "NO"]


p1 <- ggplot(data = x.all, aes(x = effect, y = -log10(wi.eBH), col = diffexpressed)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("red", "black", "blue"), 
                     labels = c("Pre Circumcision", "Not Significant", "Post Circumcision")) +
  ggtitle("ALDEx2 Pre vs Post") +
  labs(y = "-log10(p-value)", x = "Effect Size", color = "Legend") +
  geom_text_repel(aes(label = delabel), size = 5, fontface = "italic", show.legend = FALSE) +
  theme(
    legend.position = c(0.9, 0.14),
    legend.title = element_text(size = 10, face = "bold"),
    legend.background = element_rect(
      fill = "white",          # Background color of the legend
      size = 0.5,              # Border line thickness
      linetype = "solid",      # Solid border line
      colour = "black"         # Border color
    ),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),          # Remove all grid lines
    axis.line = element_line(color = "black"), # Keep the axes lines
    axis.ticks = element_line(color = "black") # Add axis ticks
  ) +
  scale_x_continuous(limits = c(-1, 1)) # Set x-axis range


p1 


#Now collecting data for my table
pre <- rel[, c("EL_003_pre", "EL_007_un", "EL_025_un1", "EL_036_pre", "EL_050_pre", "EL_051_pre", "EL_052_pre", "EL_059_pre", "EL_070_pre", "EL_074_pre", "EL_075_pre", "EL_078_pre", "NE_002_pre2", "NE_009_un", "NE_010_un1", "NE_018_pre", "NE_034_pre", "NE_035_pre", "NE_041_pre", "NE_056_pre", "NE_058_pre", "NE_063_pre", "NE_067_pre", "NE_069_pre", "NE_079_pre", "PA_001_pre", "PA_008_un", "PA_011_un1", "PA_013_pre", "PA_014_pre", "PA_015_pre", "PA_016_un1", "PA_029_pre2", "PA_030_pre2", "PA_033_pre", "PA_042_pre", "PA_044_pre", "PA_045_pre", "PA_046_pre", "PA_048_pre", "PA_049_pre", "PA_053_pre", "PA_054_pre", "PA_055_pre", "PA_057_pre", "PA_060_pre", "PA_061_pre2", "PA_064_pre", "PA_065_pre", "PA_076_pre", "PA_080_pre")]

# Count the number of samples with more than 100 reads in the "Campylobacter" row
samples_with_high_reads <- sum(pre["Campylobacter", ] > 1)
print(samples_with_high_reads)
#48

samples_with_high_reads <- sum(pre["Peptoniphilus", ] > 1)
print(samples_with_high_reads)
#51

samples_with_high_reads <- sum(pre["Varibaculum", ] > 1)
print(samples_with_high_reads)
#51

samples_with_high_reads <- sum(pre["Negativicoccus", ] > 1)
print(samples_with_high_reads)
#31

samples_with_high_reads <- sum(pre["Porphyromonas", ] > 1)
print(samples_with_high_reads)
#48

samples_with_high_reads <- sum(pre["Hoylesella", ] > 1)
print(samples_with_high_reads)
#51

samples_with_high_reads <- sum(pre["Ezakiella", ] > 1)
print(samples_with_high_reads)
#47

samples_with_high_reads <- sum(pre["Fastidiosipila", ] > 1)
print(samples_with_high_reads)
#6

samples_with_high_reads <- sum(pre["Corynebacterium", ] > 1)
print(samples_with_high_reads)
#43

samples_with_high_reads <- sum(pre["Staphylococcus", ] > 1)
print(samples_with_high_reads)
#45

samples_with_high_reads <- sum(pre["Arcanobacterium", ] > 1)
print(samples_with_high_reads)
#8

samples_with_high_reads <- sum(pre["Actinotignum", ] > 1)
print(samples_with_high_reads)
#34


post <- rel[, c("EL_003_pw", "EL_007_un1", "EL_025_un", "EL_036_pw", "EL_050_pw", "EL_051_pw", "EL_052_pw", "EL_059_pw", "EL_070_pw", "EL_074_pw", "EL_075_pw", "EL_078_pw", "NE_002_pw", "NE_009_un1", "NE_010_un", "NE_018_pw", "NE_034_pw", "EL_035_pw", "NE_041_pw", "NE_056_pw", "NE_058_pw", "NE_063_pw", "NE_067_pw", "NE_069_pw", "NE_079_pw", "PA_001_pw", "PA_008_un1", "PA_011_un", "PA_013_pw", "PA_014_pw", "PA_015_pw", "PA_016_un", "PA_029_pw", "PA_030_pw", "PA_033_pw", "PA_042_pw", "PA_044_pw", "PA_045_pw", "PA_046_pw", "PA_048_pw2", "PA_049_pw", "PA_053_pw", "PA_054_pw", "PA_055_pw", "PA_057_pw", "PA_060_pw", "PA_061_pw", "PA_064_pw", "PA_065_pw", "PA_076_pw", "PA_080_pw")]

# Count the number of samples with more than 100 reads in the "Campylobacter" row
samples_with_high_reads1 <- sum(post["Campylobacter", ] > 1)
print(samples_with_high_reads1)
#28

samples_with_high_reads1 <- sum(post["Peptoniphilus", ] > 1)
print(samples_with_high_reads1)
#50

samples_with_high_reads1 <- sum(post["Varibaculum", ] > 1)
print(samples_with_high_reads1)
#47

samples_with_high_reads1 <- sum(post["Negativicoccus", ] > 1)
print(samples_with_high_reads1)
#11

samples_with_high_reads1 <- sum(post["Porphyromonas", ] > 1)
print(samples_with_high_reads1)
#45

samples_with_high_reads1 <- sum(post["Hoylesella", ] > 1)
print(samples_with_high_reads1)
#44

samples_with_high_reads1 <- sum(post["Ezakiella", ] > 1)
print(samples_with_high_reads1)
#41

samples_with_high_reads1 <- sum(post["Fastidiosipila", ] > 1)
print(samples_with_high_reads1)
#0

samples_with_high_reads1 <- sum(post["Corynebacterium", ] > 1)
print(samples_with_high_reads1)
#50

samples_with_high_reads1 <- sum(post["Staphylococcus", ] > 1)
print(samples_with_high_reads1)
#49

samples_with_high_reads1 <- sum(post["Arcanobacterium", ] > 1)
print(samples_with_high_reads1)
#0

samples_with_high_reads1 <- sum(post["Actinotignum", ] > 1)
print(samples_with_high_reads1)
#21



# Calculate the median value for each taxon
rel <- read.table("rel_abund_138.2.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

pre1 <- rel[, c("EL_003_pre", "EL_007_un", "EL_025_un1", "EL_036_pre", "EL_050_pre", "EL_051_pre", "EL_052_pre", "EL_059_pre", "EL_070_pre", "EL_074_pre", "EL_075_pre", "EL_078_pre", "NE_002_pre2", "NE_009_un", "NE_010_un1", "NE_018_pre", "NE_034_pre", "NE_035_pre", "NE_041_pre", "NE_056_pre", "NE_058_pre", "NE_063_pre", "NE_067_pre", "NE_069_pre", "NE_079_pre", "PA_001_pre", "PA_008_un", "PA_011_un1", "PA_013_pre", "PA_014_pre", "PA_015_pre", "PA_016_un1", "PA_029_pre2", "PA_030_pre2", "PA_033_pre", "PA_042_pre", "PA_044_pre", "PA_045_pre", "PA_046_pre", "PA_048_pre", "PA_049_pre", "PA_053_pre", "PA_054_pre", "PA_055_pre", "PA_057_pre", "PA_060_pre", "PA_061_pre2", "PA_064_pre", "PA_065_pre", "PA_076_pre", "PA_080_pre")]

taxa_medians <- apply(pre1, 1, median)
print(taxa_medians)

# Calculate the standard deviation for each sample (column)
sample_sd <- apply(pre1, 1, sd)

# Print the result
print(sample_sd)



post1 <- rel[, c("EL_003_pw", "EL_007_un1", "EL_025_un", "EL_036_pw", "EL_050_pw", "EL_051_pw", "EL_052_pw", "EL_059_pw", "EL_070_pw", "EL_074_pw", "EL_075_pw", "EL_078_pw", "NE_002_pw", "NE_009_un1", "NE_010_un", "NE_018_pw", "NE_034_pw", "EL_035_pw", "NE_041_pw", "NE_056_pw", "NE_058_pw", "NE_063_pw", "NE_067_pw", "NE_069_pw", "NE_079_pw", "PA_001_pw", "PA_008_un1", "PA_011_un", "PA_013_pw", "PA_014_pw", "PA_015_pw", "PA_016_un", "PA_029_pw", "PA_030_pw", "PA_033_pw", "PA_042_pw", "PA_044_pw", "PA_045_pw", "PA_046_pw", "PA_048_pw2", "PA_049_pw", "PA_053_pw", "PA_054_pw", "PA_055_pw", "PA_057_pw", "PA_060_pw", "PA_061_pw", "PA_064_pw", "PA_065_pw", "PA_076_pw", "PA_080_pw")]

taxa_medians <- apply(post1, 1, median)
print(taxa_medians)

# Calculate the standard deviation for each sample (column)
sample_sd1 <- apply(post1, 1, sd)

# Print the result
print(sample_sd1)




#Now collecting data for my adults vs peds table
#Looking for prevalence- samples with greater than 1% relative abundance

pre <- rel[, c("PA_001_pre", "FS_006_pre1", "EL_007_un", "PA_008_un", "NE_009_un", "PA_011_un1", "PA_013_pre", "PA_014_pre", "PA_015_pre", "NE_018_pre", "PA_019_pre", "NE_021_pre", "EL_022_pre", "EL_023_pre", "FS_024_pre", "EL_026_pre", "EL_027_pre2", "EL_028_pre2", "PA_029_pre2", "PA_030_pre2", "EL_031_pre", "EL_032_pre", "PA_033_pre", "NE_034_pre", "NE_035_pre", "EL_036_pre", "NE_037_pre", "NE_038_pre", "PA_039_pre", "PA_040_pre", "NE_041_pre", "PA_042_pre", "NE_043_pre", "PA_044_pre", "PA_045_pre", "EL_047_pre", "PA_048_pre", "PA_049_pre", "EL_050_pre", "EL_051_pre", "EL_052_pre", "PA_053_pre", "PA_054_pre", "PA_055_pre", "NE_056_pre", "PA_057_pre", "NE_058_pre", "EL_059_pre", "PA_060_pre", "PA_061_pre2", "PA_062_pre", "NE_063_pre", "PA_065_pre", "NE_066_pre", "EL_068_pre2", "NE_069_pre", "NE_071_pre", "NE_072_pre", "EL_073_pre", "EL_074_pre", "EL_075_pre", "PA_076_pre", "NE.077_pre", "NE_079_pre", "PA_080_pre")]

samples_with_high_reads <- sum(pre["Acinetobacter", ] > 1)
print(samples_with_high_reads)
#0

samples_with_high_reads <- sum(pre["Actinotignum", ] > 1)
print(samples_with_high_reads)
#45

samples_with_high_reads <- sum(pre["Aerococcus", ] > 1)
print(samples_with_high_reads)
#10

samples_with_high_reads <- sum(pre["Anaerococcus", ] > 1)
print(samples_with_high_reads)
#45

samples_with_high_reads <- sum(pre["Campylobacter", ] > 1)
print(samples_with_high_reads)
#71

samples_with_high_reads <- sum(pre["Corynebacterium", ] > 1)
print(samples_with_high_reads)
#63

samples_with_high_reads <- sum(pre["Dialister", ] > 1)
print(samples_with_high_reads)
#33

samples_with_high_reads <- sum(pre["Ezakiella", ] > 1)
print(samples_with_high_reads)
#70

samples_with_high_reads <- sum(pre["Fenollaria", ] > 1)
print(samples_with_high_reads)
#53

samples_with_high_reads <- sum(pre["Finegoldia", ] > 1)
print(samples_with_high_reads)
#38

samples_with_high_reads <- sum(pre["Fusobacterium", ] > 1)
print(samples_with_high_reads)
#10

samples_with_high_reads <- sum(pre["Hoylesella", ] > 1)
print(samples_with_high_reads)
#75

samples_with_high_reads <- sum(pre["Jonquetella", ] > 1)
print(samples_with_high_reads)
#3

samples_with_high_reads <- sum(pre["Lactobacillus", ] > 1)
print(samples_with_high_reads)
#0

samples_with_high_reads <- sum(pre["Mobiluncus", ] > 1)
print(samples_with_high_reads)
#43

samples_with_high_reads <- sum(pre["Murdochiella", ] > 1)
print(samples_with_high_reads)
#8

samples_with_high_reads <- sum(pre["Negativicoccus", ] > 1)
print(samples_with_high_reads)
#45

samples_with_high_reads <- sum(pre["Peptoniphilus", ] > 1)
print(samples_with_high_reads)
#75

samples_with_high_reads <- sum(pre["Porphyromonas", ] > 1)
print(samples_with_high_reads)
#71

samples_with_high_reads <- sum(pre["Prevotella", ] > 1)
print(samples_with_high_reads)
#54

samples_with_high_reads <- sum(pre["Schaalia", ] > 1)
print(samples_with_high_reads)
#12

samples_with_high_reads <- sum(pre["Staphylococcus", ] > 1)
print(samples_with_high_reads)
#67

samples_with_high_reads <- sum(pre["Streptococcus", ] > 1)
print(samples_with_high_reads)
#9

samples_with_high_reads <- sum(pre["Varibaculum", ] > 1)
print(samples_with_high_reads)
#75

samples_with_high_reads <- sum(pre["Winkia", ] > 1)
print(samples_with_high_reads)
#3

#Now for adults
load("abund_PC_genus.RData")
pre <-abund_PC_genus*100

samples_with_high_reads <- sum(pre["Acinetobacter", ] > 1)
print(samples_with_high_reads)
#4

samples_with_high_reads <- sum(pre["Actinotignum", ] > 1)
print(samples_with_high_reads)
#6

samples_with_high_reads <- sum(pre["Aerococcus", ] > 1)
print(samples_with_high_reads)
#3

samples_with_high_reads <- sum(pre["Anaerococcus", ] > 1)
print(samples_with_high_reads)
#23

samples_with_high_reads <- sum(pre["Campylobacter", ] > 1)
print(samples_with_high_reads)
#34

samples_with_high_reads <- sum(pre["Corynebacterium", ] > 1)
print(samples_with_high_reads)
#39

samples_with_high_reads <- sum(pre["Dialister", ] > 1)
print(samples_with_high_reads)
#31

samples_with_high_reads <- sum(pre["Ezakiella", ] > 1)
print(samples_with_high_reads)
#22

samples_with_high_reads <- sum(pre["Fenollaria", ] > 1)
print(samples_with_high_reads)
#23

samples_with_high_reads <- sum(pre["Finegoldia", ] > 1)
print(samples_with_high_reads)
#40

samples_with_high_reads <- sum(pre["Fusobacterium", ] > 1)
print(samples_with_high_reads)
#8

samples_with_high_reads <- sum(pre["Hoylesella", ] > 1)
print(samples_with_high_reads)
#48

samples_with_high_reads <- sum(pre["Jonquetella", ] > 1)
print(samples_with_high_reads)
#5

samples_with_high_reads <- sum(pre["Lactobacillus", ] > 1)
print(samples_with_high_reads)
#10

samples_with_high_reads <- sum(pre["Mobiluncus", ] > 1)
print(samples_with_high_reads)
#22

samples_with_high_reads <- sum(pre["Murdochiella", ] > 1)
print(samples_with_high_reads)
#12

samples_with_high_reads <- sum(pre["Negativicoccus", ] > 1)
print(samples_with_high_reads)
#35

samples_with_high_reads <- sum(pre["Peptoniphilus", ] > 1)
print(samples_with_high_reads)
#54

samples_with_high_reads <- sum(pre["Porphyromonas", ] > 1)
print(samples_with_high_reads)
#37

samples_with_high_reads <- sum(pre["Prevotella", ] > 1)
print(samples_with_high_reads)
#46

samples_with_high_reads <- sum(pre["Schaalia", ] > 1)
print(samples_with_high_reads)
#4

samples_with_high_reads <- sum(pre["Staphylococcus", ] > 1)
print(samples_with_high_reads)
#9

samples_with_high_reads <- sum(pre["Streptococcus", ] > 1)
print(samples_with_high_reads)
#4

samples_with_high_reads <- sum(pre["Varibaculum", ] > 1)
print(samples_with_high_reads)
#13

samples_with_high_reads <- sum(pre["Winkia", ] > 1)
print(samples_with_high_reads)
#17







taxa_medians <- apply(pre, 1, median)
print(taxa_medians)

# Calculate the standard deviation for each sample (column)
sample_sd <- apply(pre, 1, sd)

# Print the result
print(sample_sd)


load("abund_PC_genus.RData")
PC <- (abund_PC_genus)*100

taxa_medians <- apply(PC, 1, median)
print(taxa_medians)

# Calculate the standard deviation for each sample (column)
sample_sd <- apply(PC, 1, sd)

# Print the result
print(format(sample_sd, scientific = FALSE))




#ALDEx2 comparing ages 0-4 to 12-18
counts <- read.table("agg_counts_138.2.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

subset <- counts[, c("PA_001_pre", "NE_002_pre2", "EL_003_pre", "FS_006_pre1", "NE_010_un1", "PA_013_pre", "PA_016_un1", "PA_019_pre", "EL_022_pre", "FS_024_pre", "EL_027_pre2", "EL_028_pre2", "EL_031_pre", "EL_032_pre", "NE_037_pre", "NE_038_pre", "PA_039_pre", "NE_041_pre", "PA_042_pre", "PA_045_pre", "PA_046_pre", "EL_050_pre", "PA_053_pre", "EL_059_pre", "PA_061_pre2", "PA_064_pre", "NE_067_pre", "EL_070_pre", "EL_073_pre", "EL_074_pre", "EL_075_pre", "NE-077_pre", "EL_078_pre", "NE_079_pre")]


counts_t <-t(subset)
# Generate a relative abundance table and remove SVs accounting for < 1% of reads in every sample
# If you care about very rare things you could lower this to 0.1% (0.001)
props <- apply(counts_t, 1, function(x) {x/sum(x)})
filt_by_props <- counts_t[, apply(props, 1, max) >= 0.01]

tc <- t(filt_by_props)

colnames(tc)[colnames(tc) == "PA_001_pre"] <- "PA_001_y"
colnames(tc)[colnames(tc) == "NE_002_pre2"] <- "NE_002_o"
colnames(tc)[colnames(tc) == "EL_003_pre"] <- "EL_003_o"
colnames(tc)[colnames(tc) == "FS_006_pre1"] <- "NE_006_y"
colnames(tc)[colnames(tc) == "NE_010_un1"] <- "NE_010_o"
colnames(tc)[colnames(tc) == "PA_013_pre"] <- "PA_013_o"
colnames(tc)[colnames(tc) == "PA_016_un1"] <- "PA_016_o"
colnames(tc)[colnames(tc) == "PA_019_pre"] <- "PA_019_y"
colnames(tc)[colnames(tc) == "EL_022_pre"] <- "EL_022_y"
colnames(tc)[colnames(tc) == "FS_024_pre"] <- "NE_024_y"
colnames(tc)[colnames(tc) == "EL_027_pre2"] <- "EL_027_o"
colnames(tc)[colnames(tc) == "EL_028_pre2"] <- "EL_028_o"
colnames(tc)[colnames(tc) == "EL_031_pre"] <- "EL_031_y"
colnames(tc)[colnames(tc) == "EL_032_pre"] <- "EL_032_y"
colnames(tc)[colnames(tc) == "NE_037_pre"] <- "NE_037_o"
colnames(tc)[colnames(tc) == "NE_038_pre"] <- "NE_038_y"
colnames(tc)[colnames(tc) == "PA_039_pre"] <- "PA_039_o"
colnames(tc)[colnames(tc) == "NE_041_pre"] <- "NE_041_o"
colnames(tc)[colnames(tc) == "PA_042_pre"] <- "PA_042_o"
colnames(tc)[colnames(tc) == "PA_045_pre"] <- "PA_045_y"
colnames(tc)[colnames(tc) == "PA_046_pre"] <- "PA_046_o"
colnames(tc)[colnames(tc) == "EL_050_pre"] <- "EL_050_o"
colnames(tc)[colnames(tc) == "EL_053_pre"] <- "EL_053_o"
colnames(tc)[colnames(tc) == "EL_059_pre"] <- "EL_059_y"
colnames(tc)[colnames(tc) == "PA_061_pre2"] <- "PA_061_o"
colnames(tc)[colnames(tc) == "PA_064_pre"] <- "PA_064_o"
colnames(tc)[colnames(tc) == "NE_067_pre"] <- "NE_067_o"
colnames(tc)[colnames(tc) == "EL_070_pre"] <- "EL_070_o"
colnames(tc)[colnames(tc) == "EL_073_pre"] <- "EL_073_y"
colnames(tc)[colnames(tc) == "EL_074_pre"] <- "EL_074_y"
colnames(tc)[colnames(tc) == "EL_075_pre"] <- "EL_075_y"
colnames(tc)[colnames(tc) == "NE-077_pre"] <- "NE_077_o"
colnames(tc)[colnames(tc) == "EL_078_pre"] <- "EL_078_o"
colnames(tc)[colnames(tc) == "NE_079_pre"] <- "NE_079_y"

#specify your groups
young<-grep("y", colnames(tc), value = TRUE)
old<-grep("o", colnames(tc), value = TRUE)


aldex.in<-tc[,c(young, old)]

conds<-c(rep("young", length(young)), rep("old", length(old)))

#get the clr values
#this is the main ALDEx function for all downstream analyses
#mc.samples=128 is often sufficient
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)

#perform t-test (both Welches and Wilcoxon, plus a Benjamini-Hochberg multiple test correction)
x.tt <- aldex.ttest(x, paired.test=FALSE)

#estimate effect size and the within and between condition values
#include indiv. samples or not
x.effect <- aldex.effect(x)

#merge the data
x.all <- data.frame(x.tt, x.effect)

#write a .txt with the results
write.table(x.all, file="aldex_ages_0-4_12-18_138.2.txt", sep="\t", quote=F, col.names=NA)



#Making volcano plot

#colours based on significance!
# add a column of NAs
x.all$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
x.all$diffexpressed[x.all$effect > 0.5] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
x.all$diffexpressed[x.all$effect < -0.5] <- "DOWN"

#labels based on significance!
x.all$delabel <- NA
x.all$delabel[x.all$diffexpressed != "NO"] <- rownames(x.all)[x.all$diffexpressed != "NO"]

library(ggrepel)
library(ggplot2)

p1 <- ggplot(data=x.all, aes(x=effect, y=-log10(wi.eBH), col=diffexpressed)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values=c("red", "black", "blue"), 
                     labels=c("Ages 12-18", "Not Significant", "Ages 0-4")) +
  ggtitle("ALDEx2 Ages 12-18 vs 0-4") +
  labs(y= "-log10(p-value)", x = "Effect Size", color = "Legend") +
  geom_text_repel(aes(label = delabel), size = 3) +
  theme(legend.position = c(0.9, 0.14), legend.title = element_text(size=10, 
                                                                    face="bold"), legend.background = element_rect(fill="white",
                                                                                                                   size=0.5, linetype="solid", colour ="black"), plot.title = element_text(hjust = 0.5)) 


p1 

#Now age-related data for my table
counts_table <- read.table("rel_abund_138.2.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)
#Age 0-4
young <- counts_table[, c("PA_001_pre", "FS_006_pre1", "PA_019_pre", "EL_022_pre", "FS_024_pre", "EL_031_pre", "EL_032_pre", "NE_038_pre", "PA_045_pre", "EL_059_pre", "EL_073_pre", "EL_074_pre", "EL_075_pre", "NE_079_pre")]

# Calculate the median of each row
young$Median <- apply(young, 1, median, na.rm = TRUE)

young$SD <- apply(young, 1, sd, na.rm = TRUE)

#Ages 12-18
old <- counts_table[, c("NE_002_pre2", "EL_003_pre", "NE_010_un1", "PA_013_pre", "PA_016_un1", "EL_027_pre2", "EL_028_pre2", "NE_037_pre", "PA_039_pre", "NE_041_pre", "PA_042_pre", "PA_046_pre", "EL_050_pre", "PA_053_pre", "PA_061_pre2", "PA_064_pre", "NE_067_pre", "EL_070_pre", "NE.077_pre", "EL_078_pre")]

# Calculate the median of each row
old$Median <- apply(old, 1, median, na.rm = TRUE)

old$SD <- apply(old, 1, sd, na.rm = TRUE)


#Now for prevalence
young1 <- counts[, c("PA_001_pre", "FS_006_pre1", "PA_019_pre", "EL_022_pre", "FS_024_pre", "EL_031_pre", "EL_032_pre", "NE_038_pre", "PA_045_pre", "EL_059_pre", "EL_073_pre", "EL_074_pre", "EL_075_pre", "NE_079_pre")]

old1 <- counts[, c("NE_002_pre2", "EL_003_pre", "NE_010_un1", "PA_013_pre", "PA_016_un1", "EL_027_pre2", "EL_028_pre2", "NE_037_pre", "PA_039_pre", "NE_041_pre", "PA_042_pre", "PA_046_pre", "EL_050_pre", "PA_053_pre", "PA_061_pre2", "PA_064_pre", "NE_067_pre", "EL_070_pre", "NE-077_pre", "EL_078_pre")]



#Check correlations for ages
library(readxl)
library(ALDEx2)

meta1 <- read_excel("Meta_Age.xlsx")

rownames(meta1) <- meta1[[1]]  # Assign the first column as row names

meta <- meta1[c("EL_003_pre", "EL_007_un", "EL_022_pre", "EL_023_pre", "EL_025_un1", "EL_026_pre", "EL_027_pre2", "EL_028_pre2", "EL_031_pre", "EL_032_pre", "EL_036_pre", "EL_047_pre", "EL_050_pre", "EL_051_pre", "EL_052_pre", "EL_059_pre", "EL_068_pre2", "EL_070_pre", "EL_073_pre", "EL_074_pre", "EL_075_pre", "EL_078_pre", "FS_006_pre1", "FS_024_pre", "NE.077_pre", "NE_002_pre2", "NE_009_un", "NE_010_un1", "NE_018_pre", "NE_021_pre", "NE_034_pre", "NE_035_pre", "NE_037_pre", "NE_038_pre", "NE_041_pre", "NE_043_pre", "NE_056_pre", "NE_058_pre", "NE_063_pre", "NE_066_pre", "NE_067_pre", "NE_069_pre", "NE_071_pre", "NE_072_pre", "NE_079_pre", "PA_001_pre", "PA_008_un", "PA_011_un1", "PA_013_pre", "PA_014_pre", "PA_015_pre", "PA_016_un1", "PA_019_pre",  "PA_029_pre2", "PA_030_pre2", "PA_033_pre", "PA_039_pre", "PA_040_pre", "PA_042_pre", "PA_044_pre", "PA_045_pre", "PA_046_pre", "PA_048_pre", "PA_049_pre", "PA_053_pre", "PA_054_pre", "PA_055_pre", "PA_057_pre", "PA_060_pre", "PA_061_pre2", "PA_062_pre", "PA_064_pre", "PA_065_pre", "PA_076_pre", "PA_080_pre"),]

clr <- aldex.clr(pre, mc.samples = 128)
corr_age <- aldex.corr(clr, meta$Age)

#order by smallest BH-corrected P-value for Spearman
corr_order1<-corr_age[order(corr_age$spearman.eBH),]
#write to file
write.table(corr_order1, file= "AGE_correlation_taxa.txt", sep= "\t", quote=F, col.names=NA)
