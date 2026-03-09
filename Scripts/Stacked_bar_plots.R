## Pathological Phimosis is Associated with Foreskin Immune Cell Infiltration but Not Microbiota Composition
## Abundance correlations and comparisons
## Prepared by: Rachel Penney
## Reviewed by: Jorge Rojas-Vargas

# Load libraries
library(tidyr)
library(dplyr)
library(ggplot2)


# Read in raw counts and taxonomic information files
counts <- read.table("rel_abund_138.2.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

#All pre samples- one per patient
pre <- counts[, c("EL_003_pre", "EL_007_un", "EL_022_pre", "EL_023_pre", "EL_025_un1", "EL_026_pre", "EL_027_pre2", "EL_028_pre2", "EL_031_pre", "EL_032_pre", "EL_036_pre", "EL_047_pre", "EL_050_pre", "EL_051_pre", "EL_052_pre", "EL_059_pre", "EL_068_pre2", "EL_070_pre", "EL_073_pre", "EL_074_pre", "EL_075_pre", "EL_078_pre", "FS_006_pre1", "FS_024_pre", "NE.077_pre", "NE_002_pre2", "NE_009_un", "NE_010_un1", "NE_018_pre", "NE_021_pre", "NE_034_pre", "NE_035_pre", "NE_037_pre", "NE_038_pre", "NE_041_pre", "NE_043_pre", "NE_056_pre", "NE_058_pre", "NE_063_pre", "NE_066_pre", "NE_067_pre", "NE_069_pre", "NE_071_pre", "NE_072_pre", "NE_079_pre", "PA_001_pre", "PA_008_un", "PA_011_un1", "PA_013_pre", "PA_014_pre", "PA_015_pre", "PA_016_un1", "PA_019_pre",  "PA_029_pre2", "PA_030_pre2", "PA_033_pre", "PA_039_pre", "PA_040_pre", "PA_042_pre", "PA_044_pre", "PA_045_pre", "PA_046_pre", "PA_048_pre", "PA_049_pre", "PA_053_pre", "PA_054_pre", "PA_055_pre", "PA_057_pre", "PA_060_pre", "PA_061_pre2", "PA_062_pre", "PA_064_pre", "PA_065_pre", "PA_076_pre", "PA_080_pre")]

# Extract the Hoylesella row
prevotella_values <- unlist(pre["Hoylesella", ])

# Order the columns by the decreasing values of Hoylesella
ordered_columns <- order(prevotella_values, decreasing = FALSE)

# Reorder the counts data frame
counts_ordered <- pre[, ordered_columns]

#Create subset for samples we actually care about in order to determine how prevalent the SVs are in the samples (>100)

#Order taxa how you want
subset_pre_20 <- counts_ordered[ c("Hoylesella", "Peptoniphilus", "Ezakiella", "Varibaculum", "Porphyromonas", "Campylobacter", "Staphylococcus", "Corynebacterium", "Prevotella", "Fenollaria", "Anaerococcus", "Finegoldia", "Actinotignum", "Mobiluncus", "Negativicoccus", "Dialister", "Streptococcus", "Fusobacterium", "Aerococcus", "Schaalia", "Winkia", "Lactobacillus", "Murdochiella", "Jonquetella", "Acinetobacter"),]


#Add Other as the 21st row

column_sums <- colSums(subset_pre_20[, 1:ncol(subset_pre_20)])

# Calculate the values for the new row
new_row <- data.frame(Genera = "Other",
                      `PA_015_pre` = 100 - column_sums[1],
                      `PA_057_pre` = 100 - column_sums[2],
                      `EL_070_pre` = 100 - column_sums[3],
                      `PA_016_un1` = 100 - column_sums[4],
                      `EL_075_pre` = 100 - column_sums[5],
                      `NE_010_un1` = 100 - column_sums[6],
                      `PA_014_pre` = 100 - column_sums[7],
                      `PA_076_pre` = 100 - column_sums[8], 
                      `PA_045_pre` = 100 - column_sums[9], 
                      `PA_030_pre2` = 100 - column_sums[10], 
                      `PA_062_pre` = 100 - column_sums[11], 
                      `FS_024_pre` = 100 - column_sums[12],                                          `EL_073_pre` = 100 - column_sums[13], 
                      `EL_023_pre` = 100 - column_sums[14], 
                      `EL_036_pre` = 100 - column_sums[15],
                      `NE_002_pre2` = 100 - column_sums[16], 
                      `NE_079_pre` = 100 - column_sums[17], 
                      `NE_071_pre` = 100 - column_sums[18],
                      `PA_061_pre2` = 100 - column_sums[19],
                      `EL_032_pre` = 100 - column_sums[20],
                      `PA_064_pre` = 100 - column_sums[21],
                      `PA_039_pre` = 100 - column_sums[22],
                      `EL_003_pre` = 100 - column_sums[23],
                      `EL_028_pre2` = 100 - column_sums[24],
                      `EL_052_pre` = 100 - column_sums[25], 
                      `NE_069_pre` = 100 - column_sums[26], 
                      `PA_029_pre2` = 100 - column_sums[27], 
                      `EL_007_un` = 100 - column_sums[28], 
                      `PA_054_pre` = 100 - column_sums[29],
                      `PA_044_pre` = 100 - column_sums[30],
                      `PA_060_pre` = 100 - column_sums[31], 
                      `NE_037_pre` = 100 - column_sums[32],
                      `EL_026_pre` = 100 - column_sums[33],
                      `NE.077_pre` = 100 - column_sums[34],
                      `NE_072_pre` = 100 - column_sums[35],                                             `PA_019_pre` = 100 - column_sums[36],
                      `PA_049_pre` = 100 - column_sums[37],
                      `PA_046_pre` = 100 - column_sums[38],
                      `FS_006_pre1` = 100 - column_sums[39],
                      `PA_065_pre` = 100 - column_sums[40],
                      `EL_025_un1` = 100 - column_sums[41], 
                      `EL_027_pre2` = 100 - column_sums[42],                                            `NE_021_pre` = 100 - column_sums[43],
                      `PA_013_pre` = 100 - column_sums[44], 
                      `NE_041_pre` = 100 - column_sums[45],
                      `NE_038_pre` = 100 - column_sums[46], 
                      `NE_066_pre` = 100 - column_sums[47],
                      `NE_043_pre` = 100 - column_sums[48],
                      `PA_055_pre` = 100 - column_sums[49],
                      `NE_035_pre` = 100 - column_sums[50],
                      `EL_068_pre2` = 100 - column_sums[51],
                      `PA_080_pre` = 100 - column_sums[52],
                      `EL_078_pre` = 100 - column_sums[53],
                      `PA_053_pre` = 100 - column_sums[54],
                      `EL_022_pre` = 100 - column_sums[55],                                             `PA_040_pre` = 100 - column_sums[56],
                      `PA_008_un` = 100 - column_sums[57],
                      `PA_033_pre` = 100 - column_sums[58],
                      `NE_067_pre` = 100 - column_sums[59],
                      `EL_031_pre` = 100 - column_sums[60],
                      `PA_001_pre` = 100 - column_sums[61], 
                      `NE_009_un` = 100 - column_sums[62],
                      `EL_059_pre` = 100 - column_sums[63],
                      `PA_042_pre` = 100 - column_sums[64], 
                      `EL_050_pre` = 100 - column_sums[65],
                      `EL_047_pre` = 100 - column_sums[66],
                      `NE_063_pre` = 100 - column_sums[67],
                      `NE_034_pre` = 100 - column_sums[68],
                      `EL_051_pre` = 100 - column_sums[69],
                      `NE_018_pre` = 100 - column_sums[70],
                      `PA_048_pre` = 100 - column_sums[71], 
                      `EL_074_pre` = 100 - column_sums[72],
                      `NE_058_pre` = 100 - column_sums[73],
                      `NE_056_pre` = 100 - column_sums[74],
                      `PA_011_un1` = 100 - column_sums[75])

rownames(new_row) <- new_row[,1]
new_row <- new_row[, -1]

# Add the new row to the bottom of the data frame
subset_pre_21 <- rbind(subset_pre_20, new_row)

#Make the genera their own column again instead of the index
subset_pre_21$Genus <- row.names(subset_pre_21)
row.names(subset_pre_21) <- NULL
subset_pre_21 <- subset_pre_21[, c(76, 1:75)]


#Convert my data into long form using number of reads

# Gather data into long format
long_pre_21 <- gather(subset_pre_21, key = "Sample", value = "Value", -Genus)

# Reorder columns
long_pre_21 <- long_pre_21[, c(2, 1, 3)]

# Rename columns
colnames(long_pre_21) <- c("Sample", "Genus", "Value")

# Reset row names
rownames(long_pre_21) <- NULL


#Make a stacked bar plot

#create your taxa colours list 
colour<-c("#0000FF" ,"#CC33FF", "#336600", "#FF99CC", "#330099",  
          "#990000", "#FF66CC", "#FFFF66", "#33cc33", "#663300", "#FF9900", "#00CCFF", "#CCFFCC", "#CC3399", "#CC9900",  "#009999", "#FF0033",  "#99FFFF",          "#996600","#CC99CC", "#660099","#33CC99", "#CC6633", "#6633FF", "#FFCCCC", "#999999")

#assign colours to the taxa
names(colour) = levels(long_pre_21$Genus)

#create sample order data frame and keep only the order and UID column 
sample_order <- long_pre_21 %>%
  filter(Genus == "Hoylesella") %>%   # Filter rows where Genus is Hoylesella
  mutate(order = row_number()) %>%    # Create an order column
  dplyr::select(Sample, order)

#creating the barplot
# create a vector of genera in the order they appear in the data frame
genus_order <- unique(long_pre_21$Genus)

# create a factor with the levels in the desired order
Taxa <- factor(long_pre_21$Genus, levels = genus_order)

# create the plot for pre swabs
inner_join(long_pre_21, sample_order, by = "Sample") %>%
  mutate(Sample = factor(Sample), 
         Sample = reorder(Sample, order)) %>%
  ggplot(aes(x = Sample, y = Value, fill = Taxa)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colour) + 
  labs(x="Patients",
       y="Relative Abundance Percentage",
       title = "Relative Abundance of Taxa in the Pre Swabs") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())



#Pre paired samples
pre <- counts[, c("EL_003_pre", "EL_007_un", "EL_025_un1", "EL_036_pre", "EL_050_pre", "EL_051_pre", "EL_052_pre", "EL_059_pre", "EL_070_pre", "EL_074_pre", "EL_075_pre", "EL_078_pre", "NE_002_pre2", "NE_009_un", "NE_010_un1", "NE_018_pre", "NE_034_pre", "NE_035_pre", "NE_041_pre", "NE_056_pre", "NE_058_pre", "NE_063_pre", "NE_067_pre", "NE_069_pre", "NE_079_pre", "PA_001_pre", "PA_008_un", "PA_011_un1", "PA_013_pre", "PA_014_pre", "PA_015_pre", "PA_016_un1", "PA_029_pre2", "PA_030_pre2", "PA_033_pre", "PA_042_pre", "PA_044_pre", "PA_045_pre", "PA_046_pre", "PA_048_pre", "PA_049_pre", "PA_053_pre", "PA_054_pre", "PA_055_pre", "PA_057_pre", "PA_060_pre", "PA_061_pre2", "PA_064_pre", "PA_065_pre", "PA_076_pre", "PA_080_pre")]

# Extract the Prevotella row
prevotella_values <- unlist(pre["Hoylesella", ])

# Order the columns by the decreasing values of Prevotella
ordered_columns <- order(prevotella_values, decreasing = FALSE)

# Reorder the counts data frame
counts_ordered <- pre[, ordered_columns]

#Create subset for samples we actually care about in order to determine how prevalent the SVs are in the samples (>100)

#Order taxa how you want
subset_pre_20 <- counts_ordered[ c("Hoylesella", "Peptoniphilus", "Ezakiella", "Varibaculum", "Porphyromonas", "Campylobacter", "Staphylococcus", "Corynebacterium", "Prevotella", "Fenollaria", "Anaerococcus", "Finegoldia", "Actinotignum", "Mobiluncus", "Negativicoccus", "Dialister", "Streptococcus", "Fusobacterium", "Aerococcus", "Schaalia", "Winkia", "Bifidobacterium", "Enterococcus"),]


#Add Other as the 21st row

column_sums <- colSums(subset_pre_20[, 1:ncol(subset_pre_20)])

# Calculate the values for the new row
new_row <- data.frame(Genera = "Other",
                      `PA_015_pre` = 100 - column_sums[1],
                      `PA_057_pre` = 100 - column_sums[2],
                      `EL_070_pre` = 100 - column_sums[3],
                      `PA_016_un1` = 100 - column_sums[4],
                      `EL_075_pre` = 100 - column_sums[5],
                      `NE_010_un1` = 100 - column_sums[6],
                      `PA_014_pre` = 100 - column_sums[7],
                      `PA_076_pre` = 100 - column_sums[8], 
                      `PA_045_pre` = 100 - column_sums[9], 
                      `PA_030_pre2` = 100 - column_sums[10], 
                      `EL_036_pre` = 100 - column_sums[11], 
                      `NE_002_pre2` = 100 - column_sums[12], 
                      `NE_079_pre` = 100 - column_sums[13], 
                      `PA_061_pre2` = 100 - column_sums[14], 
                      `PA_064_pre` = 100 - column_sums[15], 
                      `EL_003_pre` = 100 - column_sums[16],
                      `EL_052_pre` = 100 - column_sums[17], 
                      `NE_069_pre` = 100 - column_sums[18], 
                      `PA_029_pre2` = 100 - column_sums[19], 
                      `EL_007_un` = 100 - column_sums[20], 
                      `PA_054_pre` = 100 - column_sums[21],
                      `PA_044_pre` = 100 - column_sums[22],
                      `PA_060_pre` = 100 - column_sums[23], 
                      `PA_049_pre` = 100 - column_sums[24],
                      `PA_046_pre` = 100 - column_sums[25], 
                      `PA_065_pre` = 100 - column_sums[26],
                      `EL_025_un1` = 100 - column_sums[27], 
                      `PA_013_pre` = 100 - column_sums[28], 
                      `NE_041_pre` = 100 - column_sums[29], 
                      `PA_055_pre` = 100 - column_sums[30],
                      `NE_035_pre` = 100 - column_sums[31],
                      `PA_080_pre` = 100 - column_sums[32],
                      `EL_078_pre` = 100 - column_sums[33],
                      `PA_053_pre` = 100 - column_sums[34],
                      `PA_008_un` = 100 - column_sums[35],
                      `PA_033_pre` = 100 - column_sums[36], 
                      `NE_067_pre` = 100 - column_sums[37],
                      `PA_001_pre` = 100 - column_sums[38], 
                      `NE_009_un` = 100 - column_sums[39],
                      `EL_059_pre` = 100 - column_sums[40],
                      `PA_042_pre` = 100 - column_sums[41], 
                      `EL_050_pre` = 100 - column_sums[42],
                      `NE_063_pre` = 100 - column_sums[43],
                      `NE_034_pre` = 100 - column_sums[44],
                      `EL_051_pre` = 100 - column_sums[45],
                      `NE_018_pre` = 100 - column_sums[46],
                      `PA_048_pre` = 100 - column_sums[47], 
                      `EL_074_pre` = 100 - column_sums[48],
                      `NE_058_pre` = 100 - column_sums[49],
                      `NE_056_pre` = 100 - column_sums[50],
                      `PA_011_un1` = 100 - column_sums[51])

rownames(new_row) <- new_row[,1]
new_row <- new_row[, -1]

# Add the new row to the bottom of the data frame
subset_pre_21 <- rbind(subset_pre_20, new_row)

#Make the genera their own column again instead of the index
subset_pre_21$Genus <- row.names(subset_pre_21)
row.names(subset_pre_21) <- NULL
subset_pre_21 <- subset_pre_21[, c(52, 1:51)]


#Convert my data into long form using number of reads

# Gather data into long format
long_pre_21 <- gather(subset_pre_21, key = "Sample", value = "Value", -Genus)

# Reorder columns
long_pre_21 <- long_pre_21[, c(2, 1, 3)]

# Rename columns
colnames(long_pre_21) <- c("Sample", "Genus", "Value")

# Reset row names
rownames(long_pre_21) <- NULL


#Make a stacked bar plot

#create your taxa colours list 
colour<-c("#0000FF" ,"#CC33FF", "#336600", "#FF99CC", "#330099",  
          "#990000", "#FF66CC", "#FFFF66",  "#33cc33", "#663300", "#FF9900", "#00CCFF", "#CCFFCC", "#CC3399", "#CC9900",  "#009999", "#FF0033",  "#99FFFF",          "#996600","#CC99CC", "#660099","#99FF66" ,"#6666CC", "#999999")

#assign colours to the taxa
names(colour) = levels(long_pre_21$Genus)

#create sample order data frame and keep only the order and UID column 
sample_order <- long_pre_21 %>%
  filter(Genus == "Hoylesella") %>%   # Filter rows where Genus is Hoylesella
  mutate(order = row_number()) %>%    # Create an order column
  dplyr::select(Sample, order)

#creating the barplot
# create a vector of genera in the order they appear in the data frame
genus_order <- unique(long_pre_21$Genus)

# create a factor with the levels in the desired order
Taxa <- factor(long_pre_21$Genus, levels = genus_order)


# create the plot for pre swabs
inner_join(long_pre_21, sample_order, by = "Sample") %>%
  mutate(Sample = factor(Sample), 
         Sample = reorder(Sample, order)) %>%
  ggplot(aes(x = Sample, y = Value, fill = Taxa)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colour) + 
  labs(x="Patients",
       y="Relative Abundance Percentage",
       title = "Relative Abundance of Taxa in the Pre Swabs") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())


# Now peds to match taxa with adults only aged 0-14
pre <- counts[, c("PA_001_pre", "FS_006_pre1", "EL_007_un", "PA_008_un", "NE_009_un", "PA_011_un1", "PA_013_pre", "PA_014_pre", "PA_015_pre", "NE_018_pre", "PA_019_pre", "NE_021_pre", "EL_022_pre", "EL_023_pre", "FS_024_pre", "EL_026_pre", "EL_027_pre2", "EL_028_pre2", "PA_029_pre2", "PA_030_pre2", "EL_031_pre", "EL_032_pre", "PA_033_pre", "NE_034_pre", "NE_035_pre", "EL_036_pre", "NE_037_pre", "NE_038_pre", "PA_039_pre", "PA_040_pre", "NE_041_pre", "PA_042_pre", "NE_043_pre", "PA_044_pre", "PA_045_pre", "EL_047_pre", "PA_048_pre", "PA_049_pre", "EL_050_pre", "EL_051_pre", "EL_052_pre", "PA_053_pre", "PA_054_pre", "PA_055_pre", "NE_056_pre", "PA_057_pre", "NE_058_pre", "EL_059_pre", "PA_060_pre", "PA_061_pre2", "PA_062_pre", "NE_063_pre", "PA_065_pre", "NE_066_pre", "EL_068_pre2", "NE_069_pre", "NE_071_pre", "NE_072_pre", "EL_073_pre", "EL_074_pre", "EL_075_pre", "PA_076_pre", "NE.077_pre", "NE_079_pre", "PA_080_pre")]


# Extract the Prevotella row
prevotella_values <- unlist(pre["Hoylesella", ])

# Order the columns by the decreasing values of Prevotella
ordered_columns <- order(prevotella_values, decreasing = FALSE)

# Reorder the counts data frame
counts_ordered <- pre[, ordered_columns]

subset_pre_20 <- counts_ordered[ c("Hoylesella", "Peptoniphilus", "Ezakiella", "Varibaculum", "Porphyromonas", "Campylobacter", "Staphylococcus", "Corynebacterium", "Prevotella", "Fenollaria", "Anaerococcus", "Finegoldia", "Actinotignum", "Mobiluncus", "Negativicoccus", "Dialister", "Streptococcus", "Fusobacterium", "Aerococcus", "Schaalia", "Winkia", "Lactobacillus", "Murdochiella", "Jonquetella", "Acinetobacter"),]


#Add Other as the 21st row

column_sums <- colSums(subset_pre_20[, 1:ncol(subset_pre_20)])

# Calculate the values for the new row
new_row <- data.frame(Genera = "Other",
                      `PA_015_pre` = 100 - column_sums[1],
                      `PA_057_pre` = 100 - column_sums[2],
                      `EL_075_pre` = 100 - column_sums[3],
                      `PA_014_pre` = 100 - column_sums[4],
                      `PA_076_pre` = 100 - column_sums[5], 
                      `PA_045_pre` = 100 - column_sums[6], 
                      `PA_030_pre2` = 100 - column_sums[7], 
                      `PA_062_pre` = 100 - column_sums[8], 
                      `FS_024_pre` = 100 - column_sums[9],                                                                                                                                           `EL_073_pre` = 100 - column_sums[10], 
                      `EL_023_pre` = 100 - column_sums[11], 
                      `EL_036_pre` = 100 - column_sums[12], 
                      `NE_079_pre` = 100 - column_sums[13], 
                      `NE_071_pre` = 100 - column_sums[14],
                      `PA_061_pre2` = 100 - column_sums[15],
                      `EL_032_pre` = 100 - column_sums[16],
                      `PA_039_pre` = 100 - column_sums[17],
                      `EL_028_pre2` = 100 - column_sums[18],
                      `EL_052_pre` = 100 - column_sums[19], 
                      `NE_069_pre` = 100 - column_sums[20], 
                      `PA_029_pre2` = 100 - column_sums[21], 
                      `EL_007_un` = 100 - column_sums[22], 
                      `PA_054_pre` = 100 - column_sums[23],
                      `PA_044_pre` = 100 - column_sums[24],
                      `PA_060_pre` = 100 - column_sums[25], 
                      `NE_037_pre` = 100 - column_sums[26],
                      `EL_026_pre` = 100 - column_sums[27],
                      `NE.077_pre` = 100 - column_sums[28],
                      `NE_072_pre` = 100 - column_sums[29],                                                                                                                                          `PA_019_pre` = 100 - column_sums[30],
                      `PA_049_pre` = 100 - column_sums[31],
                      `FS_006_pre1` = 100 - column_sums[32],
                      `PA_065_pre` = 100 - column_sums[33], 
                      `EL_027_pre2` = 100 - column_sums[34],                                                                                                                                         `NE_021_pre` = 100 - column_sums[35],
                      `PA_013_pre` = 100 - column_sums[36], 
                      `NE_041_pre` = 100 - column_sums[37],
                      `NE_038_pre` = 100 - column_sums[38], 
                      `NE_066_pre` = 100 - column_sums[39],
                      `NE_043_pre` = 100 - column_sums[40],
                      `PA_055_pre` = 100 - column_sums[41],
                      `NE_035_pre` = 100 - column_sums[42],
                      `EL_068_pre2` = 100 - column_sums[43],
                      `PA_080_pre` = 100 - column_sums[44],
                      `PA_053_pre` = 100 - column_sums[45],
                      `EL_022_pre` = 100 - column_sums[46],                                                                                                                                          `PA_040_pre` = 100 - column_sums[47],
                      `PA_008_un` = 100 - column_sums[48],
                      `PA_033_pre` = 100 - column_sums[49],
                      `EL_031_pre` = 100 - column_sums[50],
                      `PA_001_pre` = 100 - column_sums[51], 
                      `NE_009_un` = 100 - column_sums[52],
                      `EL_059_pre` = 100 - column_sums[53],
                      `PA_042_pre` = 100 - column_sums[54], 
                      `EL_050_pre` = 100 - column_sums[55],
                      `EL_047_pre` = 100 - column_sums[56],
                      `NE_063_pre` = 100 - column_sums[57],
                      `NE_034_pre` = 100 - column_sums[58],
                      `EL_051_pre` = 100 - column_sums[59],
                      `NE_018_pre` = 100 - column_sums[60],
                      `PA_048_pre` = 100 - column_sums[61], 
                      `EL_074_pre` = 100 - column_sums[62],
                      `NE_058_pre` = 100 - column_sums[63],
                      `NE_056_pre` = 100 - column_sums[64],
                      `PA_011_un1` = 100 - column_sums[65])

rownames(new_row) <- new_row[,1]
new_row <- new_row[, -1]

# Add the new row to the bottom of the data frame
subset_pre_21 <- rbind(subset_pre_20, new_row)

#Make the genera their own column again instead of the index
subset_pre_21$Genus <- row.names(subset_pre_21)
row.names(subset_pre_21) <- NULL
subset_pre_21 <- subset_pre_21[, c(66, 1:65)]


#Convert my data into long form using number of reads

# Gather data into long format
long_pre_21 <- gather(subset_pre_21, key = "Sample", value = "Value", -Genus)

# Reorder columns
long_pre_21 <- long_pre_21[, c(2, 1, 3)]

# Rename columns
colnames(long_pre_21) <- c("Sample", "Genus", "Value")

# Reset row names
rownames(long_pre_21) <- NULL


#Make a stacked bar plot

#create your taxa colours list 
colour<-c("#0000FF" ,"#CC33FF", "#336600", "#FF99CC", "#330099",  
          "#990000", "#FF66CC", "#FFFF66", "#33cc33", "#663300", "#FF9900", "#00CCFF", "#CCFFCC", "#CC3399", "#CC9900",  "#009999", "#FF0033",  "#99FFFF",          "#996600","#CC99CC", "#660099","#33CC99", "#CC6633", "#6633FF", "#FFCCCC", "#999999")

#assign colours to the taxa
names(colour) = levels(long_pre_21$Genus)

#create sample order data frame and keep only the order and UID column 
sample_order <- long_pre_21 %>%
  filter(Genus == "Hoylesella") %>%   # Filter rows where Genus is Hoylesella
  mutate(order = row_number()) %>%    # Create an order column
  dplyr::select(Sample, order)

#creating the barplot
# create a vector of genera in the order they appear in the data frame
genus_order <- unique(long_pre_21$Genus)

# create a factor with the levels in the desired order
Taxa <- factor(long_pre_21$Genus, levels = genus_order)

# create the plot for pre swabs
inner_join(long_pre_21, sample_order, by = "Sample") %>%
  mutate(Sample = factor(Sample), 
         Sample = reorder(Sample, order)) %>%
  ggplot(aes(x = Sample, y = Value, fill = Taxa)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colour) + 
  labs(x="Patients",
       y="Relative Abundance Percentage",
       title = "Relative Abundance of Taxa in the Pre Swabs") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())





#Now repeat for post samples
post <- counts[, c("PA_015_pw", "PA_057_pw", "EL_070_pw", "PA_016_un",  "EL_075_pw", "NE_010_un", "PA_014_pw", "PA_076_pw", "PA_045_pw", "PA_030_pw", "EL_036_pw", "NE_002_pw", "NE_079_pw", "PA_061_pw", "PA_064_pw", "EL_003_pw", "EL_052_pw", "NE_069_pw", "PA_029_pw", "EL_007_un1", "PA_054_pw", "PA_044_pw", "PA_060_pw", "PA_049_pw", "PA_046_pw", "PA_065_pw", "EL_025_un", "PA_013_pw", "NE_041_pw", "PA_055_pw", "EL_035_pw", "PA_080_pw", "EL_078_pw", "PA_053_pw", "PA_008_un1", "PA_033_pw", "NE_067_pw", "PA_001_pw", "NE_009_un1", "EL_059_pw", "PA_042_pw", "EL_050_pw", "NE_063_pw", "NE_034_pw", "EL_051_pw", "NE_018_pw", "PA_048_pw1", "EL_074_pw", "NE_058_pw", "NE_056_pw", "PA_011_un")]

subset_post_20 <- post[ c("Hoylesella", "Peptoniphilus", "Ezakiella", "Varibaculum", "Porphyromonas", "Campylobacter", "Staphylococcus", "Corynebacterium", "Prevotella", "Fenollaria", "Anaerococcus", "Finegoldia", "Actinotignum", "Mobiluncus", "Negativicoccus", "Dialister", "Streptococcus", "Fusobacterium", "Aerococcus", "Schaalia", "Winkia", "Bifidobacterium", "Enterococcus"),]

#Add Other as the 21st row

column_sums <- colSums(subset_post_20[, 1:ncol(subset_post_20)])

# Calculate the values for the new row
new_row <- data.frame(Genera = "Other",
                      `PA_015_pw` = 100 - column_sums[1],
                      `PA_057_pw` = 100 - column_sums[2],
                      `EL_070_pw` = 100 - column_sums[3],
                      `PA_016_un` = 100 - column_sums[4],
                      `EL_075_pw` = 100 - column_sums[5],
                      `NE_010_un` = 100 - column_sums[6],
                      `PA_014_pw` = 100 - column_sums[7],
                      `PA_076_pw` = 100 - column_sums[8], 
                      `PA_045_pw` = 100 - column_sums[9], 
                      `PA_030_pw` = 100 - column_sums[10], 
                      `EL_036_pw` = 100 - column_sums[11], 
                      `NE_002_pw` = 100 - column_sums[12], 
                      `NE_079_pw` = 100 - column_sums[13], 
                      `PA_061_pw` = 100 - column_sums[14], 
                      `PA_064_pw` = 100 - column_sums[15], 
                      `EL_003_pw` = 100 - column_sums[16],
                      `EL_052_pw` = 100 - column_sums[17], 
                      `NE_069_pw` = 100 - column_sums[18], 
                      `PA_029_pw` = 100 - column_sums[19], 
                      `EL_007_un1` = 100 - column_sums[20], 
                      `PA_054_pw` = 100 - column_sums[21],
                      `PA_044_pw` = 100 - column_sums[22],
                      `PA_060_pw` = 100 - column_sums[23], 
                      `PA_049_pw` = 100 - column_sums[24],
                      `PA_046_pw` = 100 - column_sums[25], 
                      `PA_065_pw` = 100 - column_sums[26],
                      `EL_025_un` = 100 - column_sums[27], 
                      `PA_013_pw` = 100 - column_sums[28], 
                      `NE_041_pw` = 100 - column_sums[29], 
                      `PA_055_pw` = 100 - column_sums[30],
                      `EL_035_pw` = 100 - column_sums[31],
                      `PA_080_pw` = 100 - column_sums[32],
                      `EL_078_pw` = 100 - column_sums[33],
                      `PA_053_pw` = 100 - column_sums[34],
                      `PA_008_un1` = 100 - column_sums[35],
                      `PA_033_pw` = 100 - column_sums[36], 
                      `NE_067_pw` = 100 - column_sums[37],
                      `PA_001_pw` = 100 - column_sums[38], 
                      `NE_009_un1` = 100 - column_sums[39],
                      `EL_059_pw` = 100 - column_sums[40],
                      `PA_042_pw` = 100 - column_sums[41], 
                      `EL_050_pw` = 100 - column_sums[42],
                      `NE_063_pw` = 100 - column_sums[43],
                      `NE_034_pw` = 100 - column_sums[44],
                      `EL_051_pw` = 100 - column_sums[45],
                      `NE_018_pw` = 100 - column_sums[46],
                      `PA_048_pw1` = 100 - column_sums[47], 
                      `EL_074_pw` = 100 - column_sums[48],
                      `NE_058_pw` = 100 - column_sums[49],
                      `NE_056_pw` = 100 - column_sums[50],
                      `PA_011_un` = 100 - column_sums[51])


rownames(new_row) <- new_row[,1]
new_row <- new_row[, -1]

# Add the new row to the bottom of the data frame
subset_post_21 <- rbind(subset_post_20, new_row)

#Make the genera their own column again instead of the index
subset_post_21$Genus <- row.names(subset_post_21)
row.names(subset_post_21) <- NULL
subset_post_21 <- subset_post_21[, c(52, 1:51)]


#Convert my data into long form using number of reads

# Gather data into long format
long_post_21 <- gather(subset_post_21, key = "Sample", value = "Value", -Genus)

# Reorder columns
long_post_21 <- long_post_21[, c(2, 1, 3)]

# Rename columns
colnames(long_post_21) <- c("Sample", "Genus", "Value")

# Reset row names
rownames(long_post_21) <- NULL


#Make a stacked bar plot

#create your taxa colours list 
colour<-c("#0000FF" ,"#CC33FF", "#336600", "#FF99CC", "#330099",  
          "#990000", "#FF66CC", "#FFFF66",  "#33cc33", "#663300", "#FF9900", "#00CCFF", "#CCFFCC", "#CC3399", "#CC9900",  "#009999", "#FF0033",  "#99FFFF",          "#996600","#CC99CC", "#660099","#99FF66" ,"#6666CC", "#999999")

names(colour) = levels(long_post_21$Genus)

#create sample order data frame and keep only the order and UID column 
sample_order <- long_post_21 %>%
  filter(Genus == "Hoylesella") %>%   # Filter rows where Genus is Hoylesella
  mutate(order = row_number()) %>%    # Create an order column
  dplyr::select(Sample, order)

#creating the barplot
# create a vector of genera in the order they appear in the data frame
genus_order <- unique(long_post_21$Genus)

# create a factor with the levels in the desired order
Taxa <- factor(long_post_21$Genus, levels = genus_order)

# create the plot for pre swabs
inner_join(long_post_21, sample_order, by = "Sample") %>%
  mutate(Sample = factor(Sample), 
         Sample = reorder(Sample, order)) %>%
  ggplot(aes(x = Sample, y = Value, fill = Taxa)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colour) + 
  labs(x="Patients",
       y="Relative Abundance Percentage",
       title = "Relative Abundance of Taxa in the Post Swabs") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())




#Now repeat with adults- taxa in wrong order so will not use for paper
load("abund_PC_genus.RData")
counts_ordered <- (abund_PC_genus)*100

# Get the row names (genera)
genera1 <- rownames(counts_ordered)

# Calculate the row sums (total relative abundance for each genus)
row_sums <- apply(counts_ordered[, 2:ncol(counts_ordered)], 1, sum)

# Order the genera by descending order of row sums
order_indices <- order(row_sums, decreasing = TRUE)

# Reorder the data frame by the ordered indices
subset_PC <- counts_ordered[order_indices, ]

# Set the row names back to the original values
rownames(subset_PC) <- genera1[order_indices]

# Extract the top 21 rows
subset_PC_2 <- head(subset_PC, n = 21)

# Extract the Prevotella row
prevotella_values <- as.numeric(subset_PC_2["Hoylesella", ])

# Order the columns by the decreasing values of Prevotella
ordered_columns <- order(prevotella_values, decreasing = FALSE)

# Reorder the counts data frame
subset_PC_20 <- subset_PC_2[, ordered_columns]


#Add Other as the 21st row

column_sums <- colSums(subset_PC_20[, 1:ncol(subset_PC_20)])

# Calculate the values for the new row
new_row <- data.frame(Genera = "Other",
                      `PC_018` = 100 - column_sums[1],
                      `PC_052` = 100 - column_sums[2],
                      `PC_065` = 100 - column_sums[3],
                      `PC_021` = 100 - column_sums[4],
                      `PC_026` = 100 - column_sums[5],
                      `PC_037` = 100 - column_sums[6],
                      `PC_011` = 100 - column_sums[7],
                      `PC_064` = 100 - column_sums[8], 
                      `PC_043` = 100 - column_sums[9], 
                      `PC_056` = 100 - column_sums[10], 
                      `PC_017` = 100 - column_sums[11], 
                      `PC_053` = 100 - column_sums[12], 
                      `PC_049` = 100 - column_sums[13], 
                      `PC_016` = 100 - column_sums[14],
                      `PC_109` = 100 - column_sums[15], 
                      `PC_070` = 100 - column_sums[16], 
                      `PC_029` = 100 - column_sums[17],
                      `PC_066` = 100 - column_sums[18],
                      `PC_030` = 100 - column_sums[19], 
                      `PC_007` = 100 - column_sums[20], 
                      `PC_051` = 100 - column_sums[21],
                      `PC_027` = 100 - column_sums[22],
                      `PC_023` = 100 - column_sums[23],
                      `PC_069` = 100 - column_sums[24],
                      `PC_040` = 100 - column_sums[25],
                      `PC_067` = 100 - column_sums[26],
                      `PC_004` = 100 - column_sums[27], 
                      `PC_062` = 100 - column_sums[28],
                      `PC_032` = 100 - column_sums[29],
                      `PC_003` = 100 - column_sums[30],
                      `PC_041` = 100 - column_sums[31], 
                      `PC_014` = 100 - column_sums[32],
                      `PC_068` = 100 - column_sums[33],
                      `PC_063` = 100 - column_sums[34],
                      `PC_006` = 100 - column_sums[35], 
                      `PC_058` = 100 - column_sums[36],
                      `PC_057` = 100 - column_sums[37],
                      `PC_031` = 100 - column_sums[38],
                      `PC_012` = 100 - column_sums[39],
                      `PC_054` = 100 - column_sums[40],
                      `PC_022` = 100 - column_sums[41],
                      `PC_028` = 100 - column_sums[42],
                      `PC_024` = 100 - column_sums[43],
                      `PC_036` = 100 - column_sums[44],
                      `PC_025` = 100 - column_sums[45],
                      `PC_060` = 100 - column_sums[46],
                      `PC_042` = 100 - column_sums[47],
                      `PC_038` = 100 - column_sums[48],
                      `PC_034` = 100 - column_sums[49],
                      `PC_045` = 100 - column_sums[50],                     				                `PC_108` = 100 - column_sums[51], 
                      `PC_039` = 100 - column_sums[52],
                      `PC_020` = 100 - column_sums[53],
                      `PC_033` = 100 - column_sums[54],
                      `PC_010` = 100 - column_sums[55],
                      `PC_061` = 100 - column_sums[56])


rownames(new_row) <- new_row[,1]
new_row <- new_row[, -1]

# Add the new row to the bottom of the data frame
subset_PC_21 <- rbind(subset_PC_20, new_row)

# Add the values of "Unassigned" to "Other"
subset_PC_21["Other", ] <- subset_PC_21["Other", ] + subset_PC_21["Unassigned", ]

# Remove the "Unassigned" row
subset_PC_21 <- subset_PC_21[rownames(subset_PC_21) != "Unassigned", ]

#Make the genera their own column again instead of the index
subset_PC_21$Genus <- row.names(subset_PC_21)
row.names(subset_PC_21) <- NULL
subset_PC_21 <- subset_PC_21[, c(57, 1:56)]


#Convert my data into long form using number of reads

# Gather data into long format
long_PC_21 <- gather(subset_PC_21, key = "Sample", value = "Value", -Genus)

# Reorder columns
long_PC_21 <- long_PC_21[, c(2, 1, 3)]

# Rename columns
colnames(long_PC_21) <- c("Sample", "Genus", "Value")

# Reset row names
rownames(long_PC_21) <- NULL


#Make a stacked bar plot

#create your taxa colours list 
colour<-c("#0000FF" ,"#CC33FF", "#33cc33", "#FFFF66", "#00CCFF","#330099",  "#990000",  "#663300", "#CC9900", "#FF0033", "#336600", "#CCFFCC","#FF9900", "#006666",  "#FF66CC","#CC99CC","#FF99CC", "#660066", "#FF6666", "#99CCFF", "#999999")

#assign colours to the taxa
names(colour) = levels(long_PC_21$Genus)

#create sample order data frame and keep only the order and UID column 
sample_order<-long_PC_21%>%
  filter(Genus=="Hoylesella")%>%
  mutate(order = (1:nrow(.)))%>%
  dplyr::select(Sample, order)

#creating the barplot
# create a vector of genera in the order they appear in the data frame
genus_order <- unique(long_PC_21$Genus)

# create a factor with the levels in the desired order
Taxa <- factor(long_PC_21$Genus, levels = genus_order)

# create the plot for pre swabs
inner_join(long_PC_21, sample_order, by = "Sample") %>%
  mutate(Sample = factor(Sample), 
         Sample = reorder(Sample, order)) %>%
  ggplot(aes(x = Sample, y = Value, fill = Taxa)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colour) + 
  labs(x="Patients",
       y="Relative Abundance Percentage",
       title = "Relative Abundance of Taxa in the Adult Swabs") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())






#Adults same order as uncirc peds
#Now repeat with adults
load("abund_PC_genus.RData")
counts_ordered <- (abund_PC_genus)*100

subset_PC_2 <- counts_ordered[ c("Hoylesella", "Peptoniphilus", "Ezakiella", "Varibaculum", "Porphyromonas", "Campylobacter", "Staphylococcus", "Corynebacterium", "Prevotella", "Fenollaria", "Anaerococcus", "Finegoldia", "Actinotignum", "Mobiluncus", "Negativicoccus", "Dialister", "Streptococcus", "Fusobacterium", "Aerococcus", "Schaalia", "Winkia", "Lactobacillus", "Murdochiella", "Jonquetella", "Acinetobacter"),]


# Extract the Prevotella row
prevotella_values <- as.numeric(subset_PC_2["Hoylesella", ])

# Order the columns by the decreasing values of Prevotella
ordered_columns <- order(prevotella_values, decreasing = FALSE)

# Reorder the counts data frame
subset_PC_20 <- subset_PC_2[, ordered_columns]


#Add Other as the 21st row

column_sums <- colSums(subset_PC_20[, 1:ncol(subset_PC_20)])

# Calculate the values for the new row
new_row <- data.frame(Genera = "Other",
                      `PC_018` = 100 - column_sums[1],
                      `PC_052` = 100 - column_sums[2],
                      `PC_065` = 100 - column_sums[3],
                      `PC_021` = 100 - column_sums[4],
                      `PC_026` = 100 - column_sums[5],
                      `PC_037` = 100 - column_sums[6],
                      `PC_011` = 100 - column_sums[7],
                      `PC_064` = 100 - column_sums[8], 
                      `PC_043` = 100 - column_sums[9], 
                      `PC_056` = 100 - column_sums[10], 
                      `PC_017` = 100 - column_sums[11], 
                      `PC_053` = 100 - column_sums[12], 
                      `PC_049` = 100 - column_sums[13], 
                      `PC_016` = 100 - column_sums[14],
                      `PC_109` = 100 - column_sums[15], 
                      `PC_070` = 100 - column_sums[16], 
                      `PC_029` = 100 - column_sums[17],
                      `PC_066` = 100 - column_sums[18],
                      `PC_030` = 100 - column_sums[19], 
                      `PC_007` = 100 - column_sums[20], 
                      `PC_051` = 100 - column_sums[21],
                      `PC_027` = 100 - column_sums[22],
                      `PC_023` = 100 - column_sums[23],
                      `PC_069` = 100 - column_sums[24],
                      `PC_040` = 100 - column_sums[25],
                      `PC_067` = 100 - column_sums[26],
                      `PC_004` = 100 - column_sums[27], 
                      `PC_062` = 100 - column_sums[28],
                      `PC_032` = 100 - column_sums[29],
                      `PC_003` = 100 - column_sums[30],
                      `PC_041` = 100 - column_sums[31], 
                      `PC_014` = 100 - column_sums[32],
                      `PC_068` = 100 - column_sums[33],
                      `PC_063` = 100 - column_sums[34],
                      `PC_006` = 100 - column_sums[35], 
                      `PC_058` = 100 - column_sums[36],
                      `PC_057` = 100 - column_sums[37],
                      `PC_031` = 100 - column_sums[38],
                      `PC_012` = 100 - column_sums[39],
                      `PC_054` = 100 - column_sums[40],
                      `PC_022` = 100 - column_sums[41],
                      `PC_028` = 100 - column_sums[42],
                      `PC_024` = 100 - column_sums[43],
                      `PC_036` = 100 - column_sums[44],
                      `PC_025` = 100 - column_sums[45],
                      `PC_060` = 100 - column_sums[46],
                      `PC_042` = 100 - column_sums[47],
                      `PC_038` = 100 - column_sums[48],
                      `PC_034` = 100 - column_sums[49],
                      `PC_045` = 100 - column_sums[50],                     				                  `PC_108` = 100 - column_sums[51], 
                      `PC_039` = 100 - column_sums[52],
                      `PC_020` = 100 - column_sums[53],
                      `PC_033` = 100 - column_sums[54],
                      `PC_010` = 100 - column_sums[55],
                      `PC_061` = 100 - column_sums[56])


rownames(new_row) <- new_row[,1]
new_row <- new_row[, -1]

# Add the new row to the bottom of the data frame
subset_PC_21 <- rbind(subset_PC_20, new_row)

#Make the genera their own column again instead of the index
subset_PC_21$Genus <- row.names(subset_PC_21)
row.names(subset_PC_21) <- NULL
subset_PC_21 <- subset_PC_21[, c(57, 1:56)]


#Convert my data into long form using number of reads

# Gather data into long format
long_PC_21 <- gather(subset_PC_21, key = "Sample", value = "Value", -Genus)

# Reorder columns
long_PC_21 <- long_PC_21[, c(2, 1, 3)]

# Rename columns
colnames(long_PC_21) <- c("Sample", "Genus", "Value")

# Reset row names
rownames(long_PC_21) <- NULL


#Make a stacked bar plot

#create your taxa colours list 
colour<-c("#0000FF" ,"#CC33FF", "#336600", "#FF99CC", "#330099",  
          "#990000", "#FF66CC", "#FFFF66", "#33cc33", "#663300", "#FF9900", "#00CCFF", "#CCFFCC", "#CC3399", "#CC9900",  "#009999", "#FF0033",  "#99FFFF",          "#996600","#CC99CC", "#660099","#33CC99", "#CC6633", "#6633FF", "#FFCCCC", "#999999")
#assign colours to the taxa
names(colour) = levels(long_PC_21$Genus)

#create sample order data frame and keep only the order and UID column 
sample_order<-long_PC_21%>%
  filter(Genus=="Hoylesella")%>%
  mutate(order = (1:nrow(.)))%>%
  dplyr::select(Sample, order)

#creating the barplot
# create a vector of genera in the order they appear in the data frame
genus_order <- unique(long_PC_21$Genus)

# create a factor with the levels in the desired order
Taxa <- factor(long_PC_21$Genus, levels = genus_order)

# create the plot for pre swabs
inner_join(long_PC_21, sample_order, by = "Sample") %>%
  mutate(Sample = factor(Sample), 
         Sample = reorder(Sample, order)) %>%
  ggplot(aes(x = Sample, y = Value, fill = Taxa)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colour) + 
  labs(x="Patients",
       y="Relative Abundance Percentage",
       title = "Relative Abundance of Taxa in the Adult Swabs") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())


